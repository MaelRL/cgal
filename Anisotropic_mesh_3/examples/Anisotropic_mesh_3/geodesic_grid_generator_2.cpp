// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <boost/array.hpp>
#include <boost/chrono.hpp>
#include <Eigen/Dense>
#include <omp.h>

#include <ctime>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef typename K::Segment_3                                Segment;
typedef typename K::Vector_3                                 Vector;
typedef typename K::Triangle_3                               Triangle;
typedef typename K::Line_3                                   Line;

typedef Stretched_Delaunay_3<K>                              Star;
typedef typename Star::Point_3                               Point_3;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef typename Eigen::Matrix<double, 3, 1>                 Vector3d;

typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
//typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;
typedef boost::array<std::size_t, 4>                         Tet;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_3                             EPoint;
typedef typename KExact::Segment_3                           ESegment;
typedef typename KExact::Line_3                              ELine;
typedef typename KExact::Triangle_3                          ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 2
const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 100;

// grid related stuff
Point_3 center(1., 1., 1.);
const FT grid_side = 2.;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT offset_z = center.z() - grid_side/2.;
FT n = 200.; // number of points per side of the grid
FT sq_n = n*n;
FT step = grid_side / (n-1);

// the metric field and the seeds
MF mf;
std::vector<Point_3> seeds;
std::vector<Metric> seeds_m;

//refinement
int n_refine = 0;

// debug & stuff
std::size_t known_count=0, trial_count=0, far_count=0;
std::clock_t start;

// -----------------------------------------------------------------------------
std::size_t fact(std::size_t n)
{
  return (n==0)?1:n*fact(n-1);
}

// Utility functions below...
template<std::size_t k_max>
std::vector<boost::array<std::size_t, k_max> >
combinations(const std::set<std::size_t>& set,
             std::set<std::size_t>::const_iterator it,
             const std::size_t k)
{
  // the idea is to build the combinations recursively :
  // at one level we pick one element, and we'll then build the combination
  // as this element added to the k-1 sized combinations

  CGAL_assertion(set.size() >= k_max && k > 0);
//    if(k_max == k)
//      combis.reserve(fact(n)/(fact(k)*fact(n-k))); // is there really a point?

  std::vector<boost::array<std::size_t, k_max> > combis; // current level combis

  std::set<std::size_t>::const_iterator sit = it;
  std::set<std::size_t>::const_iterator send = set.end();

  // we still need to pick k-1 elements afterwards so we need to keep some space
  std::advance(send, -(k-1));

  for(; sit!=send; ++sit)
  {
    if(k == 1)
    {
      boost::array<std::size_t, k_max> combi;

      //this is only useful to debug, but it should be bug free now!
//        for(std::size_t i=0; i<k_max; ++i)
//          combi[i] = static_cast<std::size_t>(-1);

      // initialize the last element
      combi[k_max-k] = *sit;
      combis.push_back(combi);
    }
    else
    {
      // initialize the k_max-k element
      std::vector<boost::array<std::size_t, k_max> > lower_level_combis =
                                               combinations<k_max>(set, ++sit, k-1);
      --sit;
      for(std::size_t i=0; i<lower_level_combis.size(); ++i)
        lower_level_combis[i][k_max-k] = *(sit);

      combis.insert(combis.end(), lower_level_combis.begin(),
                                  lower_level_combis.end());
    }
  }
  return combis;
}

enum FMM_state
{
  KNOWN = 0,
  TRIAL,
  FAR
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL
};

template<typename Gp>
struct Grid_point_comparer
{
  bool operator()(Gp const * const gp1, Gp const * const gp2)
  {
    return gp1->distance_to_closest_seed > gp2->distance_to_closest_seed;
  }
};

int insert_new_seed(const FT x, const FT y, const FT z)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. ||
       z < offset_z || z > center.z() + grid_side/2.)
    {
#if (verbose > 1)
      std::cout << "filtered : " << x << " " << y << " " << z << std::endl;
#endif
      return seeds.size();
    }
#endif
#if (verbose > 0)
  std::cout << "added new seed: " << x << " " << y << " " << z << std::endl;
#endif
  seeds.push_back(Point_3(x, y, z));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in("bambimboum.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y, r_z;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
#if (verbose > 10)
  std::cout << "file has : " << nv << " pts" << std::endl;
#endif
  assert(dim == 3);

  std::size_t min_nv = (std::min)(nv, vertices_nv);

  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> r_z >> useless;
    insert_new_seed(r_x, r_y, r_z);

//    insert_new_seed(center.x(), center.y(), center.z()+grid_side/2.);
//    insert_new_seed(offset_x, offset_y, offset_z);
//    insert_new_seed(offset_x+grid_side, offset_y, offset_z);
//    insert_new_seed(offset_x, offset_y+grid_side, offset_z);
//    insert_new_seed(offset_x+grid_side, offset_y+grid_side, offset_z);

//    for(int i=0; i<6; ++i)
//      for(int j=0; j<6; ++j)
//        for(int k=0; k<6; ++k)
//          insert_new_seed(center.x() - grid_side/2. + i * grid_side/6.,
//                          center.y() - grid_side/2. + j * grid_side/6.,
//                          center.z() - grid_side/2. + k * grid_side/6.);

    if(seeds.size() >= vertices_nv)
      break;
  }
#if (verbose > 0)
  std::cout << "seeds: " << seeds.size() << std::endl;
#endif
  return seeds.size();
}

struct Grid_point
{
  typedef boost::array<Grid_point*, 6> Neighbors;

  Point_3 point;
  std::size_t index;
  mutable FT distance_to_closest_seed;
  mutable std::size_t closest_seed_id;
  mutable  FMM_state state;
  Neighbors neighbors; // 6 neighbors ORDERED (above, left, back, right, front, below)
                       // if no neighbor (borders for example) then NULL
  Metric metric;
  const Grid_point* ancestor; // used to debug & draw nice things

  void change_state(FMM_state new_state)
  {
    if(new_state == state)
      std::cerr << "WARNING: useless state change..." << std::endl;

    if(state == KNOWN)
      known_count--;
    else if(state == TRIAL)
      trial_count--;
    else if(state == FAR)
      far_count--;

    if(new_state == KNOWN)
      known_count++;
    else if(new_state == TRIAL)
      trial_count++;
    else if(new_state == FAR)
      far_count++;

    state = new_state;
  }

  bool compute_closest_seed(const Grid_point* anc)
  {
    // returns true if we improved the distance
    CGAL_assertion(anc->state == KNOWN);

    const int k = 8; // depth of the ancestor edge
    FT d = FT_inf;

    // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
    // stored as 'this', ancestor1, ancestor2, etc.
    boost::array<const Grid_point*, k+1> ancestor_path;
    for(int i=1; i<k+1; ++i)
      ancestor_path[i] = NULL;
    ancestor_path[0] = this;

    const Grid_point* curr_ancestor = anc;
    for(int i=1; i<=k; ++i)
    {
      // add the new segment to the ancestor path
      ancestor_path[i] = curr_ancestor;

      Vector3d ancestor_edge;
      ancestor_edge(0) = point.x() - curr_ancestor->point.x();
      ancestor_edge(1) = point.y() - curr_ancestor->point.y();
      ancestor_edge(2) = point.z() - curr_ancestor->point.z();
      FT ancestor_edge_length = ancestor_edge.norm();
      Vector3d normalized_anc_edge = ancestor_edge/ancestor_edge_length;

      // compute the distance for the current depth (i)
      FT dist_to_ancestor = 0.;
      for(int j=0; j<i; ++j) // we add a part for each edge in the path
      {
        // get the metric for the current edge
        const Grid_point* e0 = ancestor_path[j];
        const Grid_point* e1 = ancestor_path[j+1];

        CGAL_assertion(e0 && e1);

        const Metric& m0 = e0->metric;
        const Metric& m1 = e1->metric;

        Vector3d curr_edge;
        curr_edge(0) = e0->point.x() - e1->point.x();
        curr_edge(1) = e0->point.y() - e1->point.y();
        curr_edge(2) = e0->point.z() - e1->point.z();

        // interpolate between both metric and transform the normalized edge
        // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
        Eigen::Matrix3d f = 0.5*(m0.get_transformation() + m1.get_transformation());
        Vector3d transformed_curr_edge = f*normalized_anc_edge;

        FT sp = curr_edge.dot(normalized_anc_edge);
        FT l = transformed_curr_edge.norm(); // length of the normalized anc edge in the metric

        dist_to_ancestor += sp * l;
      }
      dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

      // add ancestor edge length to the distance at that ancestor
      FT dist_at_anc = curr_ancestor->distance_to_closest_seed;
      FT new_d = dist_at_anc + dist_to_ancestor;

      if(new_d < d)
        d = new_d;

      if(!curr_ancestor->ancestor) // can't go any farther up in the ancestor tree
        break;

      curr_ancestor = curr_ancestor->ancestor;
    }

    if(d < distance_to_closest_seed)
    {
      ancestor = anc;
      distance_to_closest_seed = d;
      closest_seed_id = anc->closest_seed_id;
      return true;
    }

    return false;
  }

  PQ_state update_neighbors_distances(std::vector<Grid_point*>& trial_pq) const
  {
    // consider all the known neighbors and compute the value to that
#if (verbose > 10)
    std::cout << "update neighbors of " << index << std::endl;
#endif
    CGAL_assertion(state == KNOWN);

    PQ_state pqs_ret = NOTHING_TO_DO;
    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      Grid_point* gp = neighbors[i]; // a neighbor of 'this'
      if(!gp)
        continue;
      else if(gp->state == KNOWN)
        continue; // dual_shenanigans(gp);
      else if(gp->state == TRIAL)
      {
        // no need to insert in the trial queue since it already is in it
        if(gp->compute_closest_seed(this))
          pqs_ret = REBUILD_TRIAL;
      }
      else // gp->state == FAR
      {
        CGAL_assertion(gp->state == FAR);
        gp->compute_closest_seed(this); // always returns true here so no need to test it
        gp->state = TRIAL;
        trial_pq.push_back(gp);
        std::push_heap(trial_pq.begin(), trial_pq.end(), Grid_point_comparer<Grid_point>());
      }
    }
    return pqs_ret;
  }

  bool is_contained(const Point_3& p,
                    const Grid_point* next) const
  {
    // above : n0
    //    n2
    //    |
    // n1-gp-n3
    //    |
    //    n4
    // below : n5

    Grid_point *n0 = neighbors[0], *n1 = neighbors[1], *n2 = neighbors[2],
               *n3 = neighbors[3], *n4 = neighbors[4], *n5 = neighbors[5];

    bool is_p_geq_n0_z = (n1) ? p.z() <= n0->point.z() : false;
    bool is_p_leq_n1_x = (n1) ? p.x() <= n1->point.x() : false;
    bool is_p_geq_n2_y = (n2) ? p.y() >= n2->point.y() : false;
    bool is_p_geq_n3_x = (n3) ? p.x() >= n3->point.x() : false;
    bool is_p_leq_n4_y = (n4) ? p.y() <= n4->point.y() : false;
    bool is_p_leq_n5_z = (n5) ? p.z() <= n5->point.z() : false;

    if(is_p_geq_n0_z) { next = n0; return false; }
    if(is_p_leq_n1_x) { next = n1; return false; }
    if(is_p_geq_n2_y) { next = n2; return false; }
    if(is_p_geq_n3_x) { next = n3; return false; }
    if(is_p_leq_n4_y) { next = n4; return false; }
    if(is_p_leq_n5_z) { next = n5; return false; }
    next = this; // we're inside !
    return true;
  }

  void initialize_from_point(const Point_3& p,
                             const std::size_t seed_id) const
  {
    Vector3d v;
    v(0) = p.x() - point.x();
    v(1) = p.y() - point.y();
    v(2) = p.z() - point.z();
    const Eigen::Matrix3d& m = metric.get_mat();
    FT d = std::sqrt(v.transpose() * m * v);
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
  }

  Grid_point(const Point_3& point_, const std::size_t index_)
    : point(point_), index(index_), distance_to_closest_seed(FT_inf),
      closest_seed_id(-1), state(FAR), neighbors(),
      metric(mf->compute_metric(point)), ancestor(NULL)
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = NULL;
  }
};

struct Geo_grid
{
  typedef std::vector<Grid_point*> Grid_point_vector;
  typedef std::map<std::pair<std::size_t, std::size_t>,
                   std::pair<const Grid_point*, const Grid_point*> > Edge_bisec_map;

  std::vector<Grid_point> points;
  Grid_point_vector trial_points;

  // Duality
  mutable Edge_bisec_map edge_bisectors;
  mutable std::set<Edge> dual_edges;
  mutable std::set<Tri> dual_triangles;
  mutable std::set<Tet> dual_tetrahedras;

  // Refinement
  mutable const Grid_point* refinement_point;
  mutable FT refinement_distance;

  const Grid_point* locate_point(const Point_3& p) const
  {
    bool found = false;
    const Grid_point& next = points[0];
    while(!found)
      found = (&next)->is_contained(p, &next);
    return &next;
  }

  void locate_and_initialize(const Point_3& p,
                             const std::size_t seed_id)
  {
    // WARNING: only usable for a full grid since we assume the step size is constant
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);
    int index_z = std::floor((p.z()-offset_z)/step);

    Grid_point* gp = &(points[index_z*sq_n + index_y*n + index_x]);
    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id << std::endl;
      std::cerr << "The seed is " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    if(gp->state == TRIAL) // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");

    gp->initialize_from_point(p, seed_id);
    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

#if (verbose > 10)
    std::cout << "looking for p: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    std::cout << "found gp: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << " " << gp->point.z() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << index_z << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << offset_z << std::endl;
#endif
  }

  void initialize_geo_grid()
  {
#if (verbose > 5)
    std::cout << "grid initialization" << std::endl;
#endif

    // create the grid points (either a full grid or a smart grid (todo) )
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          Point_3 p(offset_x+i*step, offset_y+j*step, offset_z+k*step); // fill from bot left to top right
          Grid_point gp(p, i + j*n + k*sq_n);
          points.push_back(gp);
        }
      }
    }

    // assign the neighbors
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          std::size_t curr_id = i + j*n + k*sq_n;
          if(k != n-1) // there is a neighbor above
            points[curr_id].neighbors[0] = &(points[i + j*n + (k+1)*sq_n]);
          if(i != 0) // there is a neighbor left
            points[curr_id].neighbors[1] = &(points[i-1 + j*n + k*sq_n]);
          if(j != n-1) // there is a neighbor back
            points[curr_id].neighbors[2] = &(points[i + (j+1)*n + k*sq_n]);
          if(i != n-1) // there is a neighbor right
            points[curr_id].neighbors[3] = &(points[i+1 + j*n + k*sq_n]);
          if(j != 0) // there is a neighbor front
            points[curr_id].neighbors[4] = &(points[i + (j-1)*n + k*sq_n]);
          if(k != 0) // there is a neighbor below
            points[curr_id].neighbors[5] = &(points[i + j*n + (k-1)*sq_n]);
        }
      }
    }
#if (verbose > 5)
    std::cout << "neighbors assigned" << std::endl;
#endif

    // seeds belong to a cube, find it
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_3& p = seeds[i];
//      const Grid_point* gp = locate_point(p); // todo
//      gp->initialize_from_point(p, i);
      locate_and_initialize(p, i);
    }

#if (verbose > 5)
    std::cout << "grid initialized" << std::endl;
#endif
  }

  void print_states() const
  {
    std::cout << "known: " << known_count;
    std::cout << " trial: " << trial_count;
    std::cout << " far: " << far_count << std::endl;
  }

  void dual_shenanigans(const Grid_point* gp) // temporary name (◕‿◕✿)
  {
    CGAL_assertion(gp->state == KNOWN);
    // check the neighbors of gp for points that have the state KNOWN
    // in these, consider those that have a closest_seed_id different than gp
    // those are points in the dual.
    // If we have never seen that dual already, then we have also found which
    // point will we draw to in the geodesic_triangulation

    std::size_t p_id = gp->closest_seed_id;
    std::set<std::size_t> simplex;
    simplex.insert(p_id);

    boost::array<Grid_point*, 6>::const_iterator it = gp->neighbors.begin(),
                                                 iend = gp->neighbors.end();
    for(; it!=iend; ++it)
    {
      const Grid_point* gq = *it;
      if(!gq || gq->state != KNOWN)
        continue;

      std::size_t q_id = gq->closest_seed_id;
      if(p_id != q_id) // we're on the dual of [PQ]
      {
        simplex.insert(q_id);

        std::pair<const Grid_point*, const Grid_point*> gp_pair;
        std::pair<std::size_t, std::size_t> id_pair;
        if(p_id < q_id)
        {
          gp_pair = std::make_pair(gp, gq);
          id_pair = std::make_pair(p_id, q_id);
        }
        else
        {
          gp_pair = std::make_pair(gq, gp);
          id_pair = std::make_pair(q_id, p_id);
        }
        std::pair<std::pair<std::size_t, std::size_t>,
            std::pair<const Grid_point*, const Grid_point*> > value =
            std::make_pair(id_pair, gp_pair);

        // the insert below only succeeds if we haven't already found an edge p_id/q_id
        // therefore the map has for value the closest points on D_pq
        edge_bisectors.insert(value);
      }
    }

    // the simplex of size ? is added to the duals
    add_simplex_to_triangulation(gp, simplex);
  }

  void geo_grid_loop()
  {
#if (verbose > 5)
    std::cout << "main loop" << std::endl;
#endif
    refinement_point = NULL; // prepare for the next refinement point
    refinement_distance = 0.;

    bool is_t_empty = trial_points.empty();

    Grid_point* gp;
    while(!is_t_empty)
    {
      if(known_count%10000 == 0)
        print_states();

#if (verbose > 5)
      std::cout << "Queue sizes. Trial: " << trial_points.size();
      std::cout << " (known: " << known_count << " out of " << sq_n*n << ")" << std::endl;
#endif

#if (verbose > 15)
      std::cout << "trial heap: " << std::endl;
      for (std::vector<Grid_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif

      gp = trial_points.front();
      CGAL_assertion(gp && gp->state == TRIAL);
      std::pop_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
      trial_points.pop_back();

#if (verbose > 10)
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << " " << gp->point.y() << " " << gp->point.z() << ")";
      std::cout << " at distance : " << gp->distance_to_closest_seed << " from " << gp->closest_seed_id << std::endl;
#endif

      gp->state = KNOWN;

      known_count++; // tmp
      dual_shenanigans(gp);

      PQ_state pqs = gp->update_neighbors_distances(trial_points);

      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

      is_t_empty = trial_points.empty();
    }

    std::cerr << "End of geo_grid_loop. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;

    cosphericity_sweeper();
  }

  void debug()
  {
    CGAL_assertion(trial_points.empty() );

    std::cout << "BRUTE FORCE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point* gp = &(points[i]);
      std::cout << "point " << i << " min distance is supposedly: ";
      std::cout << gp->distance_to_closest_seed << std::endl;
      boost::array<Grid_point*, 6>::const_iterator it = gp->neighbors.begin(),
                                                   iend = gp->neighbors.end();
      for(; it!=iend; ++it)
      {
        const Grid_point* gq = *it;
        if(gq)
          CGAL_assertion(!gp->compute_closest_seed(gq));
      }
    }

    std::cerr << "End of debug. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void build_grid()
  {
    initialize_geo_grid();
    geo_grid_loop();
  }

  void refresh_grid_point_states()
  {
    CGAL_assertion(trial_points.empty());
    for(std::size_t i=0; i<points.size(); ++i)
      points[i].state = FAR;
  }

  void refine_grid(const Point_3 new_seed)
  {
    vertices_nv = insert_new_seed(new_seed.x(), new_seed.y(), new_seed.z());
    refresh_grid_point_states(); // we can't spread from the new seed if all states are 'KNOWN'
    locate_and_initialize(new_seed, seeds.size()-1);

    reset_dual();
    geo_grid_loop();
  }

  void refine_grid_with_self_computed_ref_point()
  {
    if(!refinement_point)
    {
      std::cerr << "Couldn't find a ref point, need more initial points" << std::endl;
      output_grid("failed_geo_grid");
      std::exit(EXIT_FAILURE);
    }
    refine_grid(refinement_point->point);
  }

  void cosphericity_sweeper() const
  {
    // we define a tiny (2k+1)*(2k+1) grid ('sweeper'), and we look at the colors
    // if we have more than 5 colors at a time, we have a (quasi-) cosphericity

    std::set<boost::array<std::size_t, 5> > quasi_cosphericities;

    const int h = (std::max)(1., n/100);
    const int sweeper_size = 2*h+1;
    if(sweeper_size > n)
    {
      std::cerr << "WARNING: sweeper too big for that grid !" << std::endl;
      return;
    }

    // we move the center of the sweeper in the grid
    for(int i=h; i<n-h; ++i)
    {
      for(int j=h; j<n-h; ++j)
      {
        for(int k=h; k<n-h; ++k)
        {
          //const Grid_point& sweeper_center = points[j*n + i];

          // gather the pts
          std::set<std::size_t> colors;
          for(int lx=-h; lx<=h; ++lx)
          {
            for(int ly=-h; ly<=h; ++ly)
            {
              for(int lz=-h; lz<=h; ++lz)
              {
                const Grid_point& gp = points[(k+lz)*sq_n + (j+ly)*n + (i+lx)];
                colors.insert(gp.closest_seed_id);
              }
            }
          }

          if(colors.size() > 4)
          {
            // we have at least 5 pts in a quasi cospherical configuration...
            // get them all with combinations!

            std::vector<boost::array<std::size_t, 5> > combis =
                combinations<5>(colors, colors.begin(), 5);

            quasi_cosphericities.insert(combis.begin(), combis.end());
          }
        }
      }
    }

    std::cout << "Sweeper no sweeping ! (with size: " << h << ")" << std::endl;
    std::cout << "Found " << quasi_cosphericities.size() << " (quasi-)cosphericities" << std::endl;

    std::set<boost::array<std::size_t, 5> >::const_iterator sit =
                                                   quasi_cosphericities.begin();
    std::set<boost::array<std::size_t, 5> >::const_iterator send =
                                                     quasi_cosphericities.end();
    for(; sit!=send; ++sit)
    {
      const boost::array<std::size_t, 5>& quasi_cosphericity = *sit;
      for(std::size_t i=0; i<sit->size(); ++i)
        std::cout << quasi_cosphericity[i] << " ";
      std::cout << std::endl;
    }
    std::cout << "End of quasi cosphericities" << std::endl;
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Grid_point& gp,
                            const Grid_point* n1, const Grid_point* n2, const Grid_point* n3) const
  {
    Tet t;
    if(n1 && n2 && n3)
    {
      t[0] = gp.index;
      t[1] = n1->index;
      t[2] = n2->index;
      t[3] = n3->index;
      std::sort(t.begin(), t.end()); // is that really needed ?
      simplices.insert(t);
    }
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Grid_point& gp,
                            std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return add_tet_to_simplices(simplices, gp, gp.neighbors[i1], gp.neighbors[i2], gp.neighbors[i3]);
  }

  void output_grid(const std::string str_base) const
  {
#if (verbose > 0)
    std::cout << "Output grid" << std::endl;
#endif

    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << points.size() << std::endl;
    out_bb << "3 1 " << points.size() << " 2" << std::endl;

    std::set<Tet> simplices;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << gp.point.z() << " " << i+1 << std::endl;

//      out_bb << gp.closest_seed_id << std::endl;
      out_bb << gp.distance_to_closest_seed << std::endl;

      // compute the tets...
      const boost::array<Grid_point*, 6>& ns = gp.neighbors; // just to get a shorter name...
      for(std::size_t i=0; i<ns.size(); ++i)
      {
#if 0
        // lazily : the eight corner tetrahedra (creates self intersections,
        // but it doesn't really matter except visually)
        add_tet_to_simplices(simplices, gp, 0, 1, 2);
        add_tet_to_simplices(simplices, gp, 0, 2, 3);
        add_tet_to_simplices(simplices, gp, 0, 3, 4);
        add_tet_to_simplices(simplices, gp, 0, 4, 1);
        add_tet_to_simplices(simplices, gp, 5, 1, 2);
        add_tet_to_simplices(simplices, gp, 5, 2, 3);
        add_tet_to_simplices(simplices, gp, 5, 3, 4);
        add_tet_to_simplices(simplices, gp, 5, 4, 1);
#else
        // a bit smarter : a cube is decomposed in 5 tets and we pick them
        // correctly so we have no overlap
        add_tet_to_simplices(simplices, gp, 0, 2, 3);
        add_tet_to_simplices(simplices, gp, 0, 1, 4);
        add_tet_to_simplices(simplices, gp, 3, 4, 5);
        add_tet_to_simplices(simplices, gp, 1, 2, 5);

        // the last one is a bit annoying since we have to grab neighbors of neighbors
        const Grid_point* n4 = gp.neighbors[4];
        const Grid_point* n3 = gp.neighbors[3];
        if(n4 && n3)
          add_tet_to_simplices(simplices, gp, n4->neighbors[3],
                               n4->neighbors[0], n3->neighbors[0]);
#endif
      }
    }

    out << "Tetrahedra" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      std::set<std::size_t> materials;

      for(std::size_t i=0; i<tet.size(); ++i)
      {
        out << tet[i] + 1 << " ";
        materials.insert(points[tet[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }

    std::set<Tri> triangles;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      Tri tr; // could use a permutation to make it look nicer, I suppose...
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = triangles.begin();
                                         it != triangles.end(); ++it)
    {
      const Tri& tr = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tr.size(); ++i)
      {
        out << tr[i] + 1 << " ";
        materials.insert(points[tr[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }
  }

  void add_triangles_edges_to_dual_edges(const Tri& tr) const
  {
    // no need to sort as long as we keep the indices in order
    Edge e;
    e[0] = tr[0]; e[1] = tr[1]; dual_edges.insert(e);
    e[0] = tr[0]; e[1] = tr[2]; dual_edges.insert(e);
    e[0] = tr[1]; e[1] = tr[2]; dual_edges.insert(e);
  }

  void add_tets_triangles_and_edges_to_dual_triangles(const Tet& tet) const
  {
    // no need to sort as long as we keep the indices in order
    Edge e;
    e[0] = tet[0]; e[1] = tet[1]; dual_edges.insert(e);
    e[0] = tet[0]; e[1] = tet[2]; dual_edges.insert(e);
    e[0] = tet[0]; e[1] = tet[3]; dual_edges.insert(e);
    e[0] = tet[1]; e[1] = tet[2]; dual_edges.insert(e);
    e[0] = tet[1]; e[1] = tet[3]; dual_edges.insert(e);
    e[0] = tet[2]; e[1] = tet[3]; dual_edges.insert(e);

    Tri tr;
    tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; dual_triangles.insert(tr);
    tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; dual_triangles.insert(tr);
    tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; dual_triangles.insert(tr);
    tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; dual_triangles.insert(tr);
  }

  template<typename Iterator>
  void add_tets_triangles_and_edges_to_dual_triangles(Iterator it, Iterator end) const
  {
    for(; it!=end; ++it)
    {
      const Tet& tet = *it;
      add_tets_triangles_and_edges_to_dual_triangles(tet);
    }
  }

  void add_simplex_to_triangulation(const Grid_point* gp,
                                    const std::set<std::size_t>& dual_simplex) const
  {
    if(dual_simplex.size() <= 1)
      return;

    // to get the next refinement point (farthest point in a dual...)
    if(gp->distance_to_closest_seed > refinement_distance)
    {
      refinement_distance = gp->distance_to_closest_seed;
      refinement_point = gp;
    }

    // Not sure if we should allow dual of edges/triangles to be considered...
    // Maybe only allow the dual of an edge to be considered if the dual is
    // on the border of the domain ?

    if(dual_simplex.size() == 2) // an edge!
    {
      typename Simplex::const_iterator it = dual_simplex.begin();
      Edge e; e[0] = *it; e[1] = (*++it);
      dual_edges.insert(e);
    }
    else if(dual_simplex.size() == 3) // a triangle!
    {
      typename Simplex::const_iterator it = dual_simplex.begin();
      Tri tr; tr[0] = *it; tr[1] = (*++it); tr[2] = (*++it);
      dual_triangles.insert(tr);

      add_triangles_edges_to_dual_edges(tr);
    }
    else if(dual_simplex.size() == 4)
    {
      typename Simplex::const_iterator it = dual_simplex.begin();
      Tet tet; tet[0] = *it; tet[1] = (*++it); tet[2] = (*++it); tet[3] = (*++it);
      dual_tetrahedras.insert(tet);

      add_tets_triangles_and_edges_to_dual_triangles(tet);
    }
    else // dual is made of > 4 points (that's a cosphericity)
    {
      std::cerr << "WARNING: COSPHERICITY... " << dual_simplex.size() << std::endl;

#if 0
      // insert all combinations of tetrahedras... (this will create illegal intersections)
      std::vector<boost::array<std::size_t, 4> > combis =
                         combinations<4>(dual_simplex, dual_simplex.begin(), 4);
      dual_tetrahedras.insert(combis.begin(), combis.end());
      add_tets_triangles_and_edges_to_dual_triangles(combis.begin(), combis.end());
#else
      // insert only the four closest


#endif
    }
  }

  void reset_dual() const
  {
    dual_edges.clear();
    dual_triangles.clear();
    dual_tetrahedras.clear();
  }

  void compute_dual() const
  {
#if (verbose > 0)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      Simplex dual_simplex;
      dual_simplex.insert(gp.closest_seed_id);

      boost::array<Grid_point*, 6>::const_iterator it = gp.neighbors.begin(),
                                                   iend = gp.neighbors.end();
      for(; it!=iend; ++it)
      {
        const Grid_point* gq = *it;
        if(!gq)
          continue;

        dual_simplex.insert(gq->closest_seed_id);
      }
      add_simplex_to_triangulation(&gp, dual_simplex);
    }
  }

  void output_geodesic_dual(std::ostream& out) const
  {
    // todo... when what a 3D geodesic triangle is well defined !
  }

  bool is_triangle_intersected(const Tri& tri, const std::set<Edge>& edges) const
  {
    bool is_intersected = false;

    // need to switch to epeck for the correct intersections here...
    const std::size_t i0 = tri[0];
    const std::size_t i1 = tri[1];
    const std::size_t i2 = tri[2];
    const Point_3& p0 = seeds[i0];
    const Point_3& p1 = seeds[i1];
    const Point_3& p2 = seeds[i2];
    const Triangle triangle(p0, p1, p2);

//#pragma omp parallel shared(is_intersected, p0, p1, p2) // hardly noticeable increase of speed...
{
    for(std::set<Edge>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
    {
#pragma omp single nowait // hack to parallelize the std::set
{
#pragma omp flush (is_intersected)
      if(!is_intersected) // hack because we're not allowed to 'break' (or 'continue')
                          // inside a pragma omp for
      {
        const Edge& edge = *it;
        bool local_intersected = false;

        const Segment segment(seeds[edge[0]], seeds[edge[1]]);

        ESegment esegment = to_exact(segment);
        ETriangle etriangle = to_exact(triangle);
        CGAL::cpp11::result_of<KExact::Intersect_3(ESegment, ETriangle)>::type
            result = CGAL::intersection(esegment, etriangle);

        if (result)
        {
          if (const EPoint* p = boost::get<EPoint>(&*result))
          {
            const EPoint& ep = *p;
            local_intersected = (ep!=to_exact(p0) && ep!=to_exact(p1) && ep!=to_exact(p2));
          }
          else if(const ESegment* s = boost::get<ESegment>(&*result))
          {
            const EPoint& ep0 = to_exact(p0);
            const EPoint& ep1 = to_exact(p1);
            const EPoint& ep2 = to_exact(p2);
            local_intersected = ((s->source()!=ep0 && s->source()!=ep1 && s->source()!=ep2) ||
                                 (s->target()!=ep0 && s->target()!=ep1 && s->target()!=ep2));
          }
        }
#pragma omp critical
{
        if(local_intersected)
          is_intersected = true;
}
#pragma omp flush (is_intersected)
      } // end !is_intersected
} // end single nowait
    } // end for
} // end parallel
    return is_intersected;
  }

  void output_straight_dual(const std::string str_base) const
  {
#if (verbose>0)
    std::cout << "Output dual" << std::endl;
#endif

    if(dual_edges.empty() && dual_triangles.empty() && dual_tetrahedras.empty())
      compute_dual();

#if (verbose > 5)
    std::cout << "captured: ";
    std::cout << dual_edges.size() << " edges, ";
    std::cout << dual_triangles.size() << " triangles, ";
    std::cout << dual_tetrahedras.size() << " tets" << std::endl;
#endif

    std::ofstream out((str_base + "_dual.mesh").c_str());
    std::ofstream outbb((str_base + "_dual.bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << seeds[i].z() << " " << i+1 << std::endl;

    outbb << "3 1 " << dual_tetrahedras.size() + dual_triangles.size() << " 1" << std::endl;

    // TETRAHEDRA
    out << "Tetrahedra" << std::endl;
    out << dual_tetrahedras.size() << std::endl;
    for(typename std::set<Tet>::iterator it = dual_tetrahedras.begin();
                                         it != dual_tetrahedras.end(); ++it)
    {
      const Tet& tet = *it;
      for(std::size_t i=0; i<tet.size(); ++i)
        out << tet[i] + 1 << " ";
      out << "1" << std::endl;

      // very inefficient to recompute them every time but whatever
      const Metric& m0 = mf->compute_metric(seeds[tet[0]]);
      const Metric& m1 = mf->compute_metric(seeds[tet[1]]);
      const Metric& m2 = mf->compute_metric(seeds[tet[2]]);
      const Metric& m3 = mf->compute_metric(seeds[tet[3]]);

      FT gamma01 = m0.compute_distortion(m1);
      FT gamma02 = m0.compute_distortion(m2);
      FT gamma03 = m0.compute_distortion(m3);
      FT gamma12 = m1.compute_distortion(m2);
      FT gamma13 = m1.compute_distortion(m3);
      FT gamma23 = m2.compute_distortion(m3);

      FT max_dist = (std::max)((std::max)(gamma01, gamma02), gamma03);
      max_dist = (std::max)((std::max)(max_dist, gamma12), gamma13);
      max_dist = (std::max)(max_dist, gamma23);

      outbb << max_dist << std::endl;
    }

     // TRIANGLES
    out << "Triangles" << std::endl;
    out << dual_triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = dual_triangles.begin();
                                         it != dual_triangles.end(); ++it)
    {
      const Tri& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
      out << is_triangle_intersected(tr, dual_edges) << std::endl;

      // very inefficient to recompute them every time but whatever
      const Metric& m0 = mf->compute_metric(seeds[tr[0]]);
      const Metric& m1 = mf->compute_metric(seeds[tr[1]]);
      const Metric& m2 = mf->compute_metric(seeds[tr[2]]);

      FT gamma01 = m0.compute_distortion(m1);
      FT gamma02 = m0.compute_distortion(m2);
      FT gamma12 = m1.compute_distortion(m2);

      FT max_dist = (std::max)((std::max)(gamma01, gamma02), gamma12);
      outbb << max_dist << std::endl;
    }
    outbb << "End" << std::endl;
    out << "End" << std::endl;
  }

  void output_grid_and_dual(const std::string str_base) const
  {
//    output_grid(str_base);
    output_straight_dual(str_base);
  }

  Geo_grid()
    :
      points(),
      trial_points(),
      dual_edges(),
      dual_triangles(),
      dual_tetrahedras(),
      refinement_point(NULL),
      refinement_distance(0.)
  { }
};

void initialize()
{
  mf = new Euclidean_metric_field<K>(1., 1., 1.);
//  mf = new Custom_metric_field<K>();

  vertices_nv = build_seeds();
  CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );
}

int main(int, char**)
{
  std::cout.precision(17);
//  std::freopen("geo_grid_2_log.txt", "w", stdout);

  double duration;
  start = std::clock();

  std::srand(0);
  initialize();

  Geo_grid gg;
  gg.build_grid();

  if(n_refine)
  {
    gg.output_grid_and_dual("geo_grid_2_pre");

    for(int i=0; i<n_refine; ++i)
    {
      gg.refine_grid_with_self_computed_ref_point();
      std::ostringstream out;
      out << "geo_grid_2_ref_" << i;
      gg.output_grid_and_dual(out.str());
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cerr << "End refinement: " << duration << std::endl;
  }

  gg.output_grid_and_dual("geo_grid_2");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "Total: " << duration << std::endl;
}
