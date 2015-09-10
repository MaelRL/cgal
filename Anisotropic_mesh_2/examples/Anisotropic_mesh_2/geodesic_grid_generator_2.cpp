// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>

#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <set>
#include <utility>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef typename K::Vector_2                                 Vector;

typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;
typedef typename Star::Metric                                Metric;

typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;

typedef typename Eigen::Matrix<double, 2, 1>                 Vector2d;

typedef typename CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
//typedef typename CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;

#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 2
const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 100;

// grid related stuff
Point_2 center(0.,0.);
const FT grid_side = 2.0;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT n = 1000.; // number of points per side of the grid
FT step = grid_side / (n-1);

// the metric field and the seeds
MF mf;
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

//refinement
int n_refine = 0;

//debug & info
int known_count=0, trial_count=0, far_count=0;
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

int insert_new_seed(const FT x, const FT y)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. )
    {
#if (verbose > 1)
      std::cout << "filtered : " << x << " " << y << std::endl;
#endif
      return seeds.size();
    }
#endif
#if (verbose > 0)
  std::cout << "added new seed: " << x << " " << y << std::endl;
#endif
  seeds.push_back(Point_2(x, y));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in("bambimboum.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "nv: " << nv << std::endl;
  assert(dim == 2);

  std::size_t min_nv = (std::min)(nv, vertices_nv);

  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> useless;
    insert_new_seed(r_x, r_y);

//    FT sqrt2d2 = std::sqrt(2.)/2.;
//    insert_new_seed(center.x()+1., center.y());
//    insert_new_seed(center.x()+sqrt2d2, center.y()+sqrt2d2);
//    insert_new_seed(center.x(), center.y()+1.);
//    insert_new_seed(center.x()-sqrt2d2, center.y()+sqrt2d2);
//    insert_new_seed(center.x()-1., center.y());
//    insert_new_seed(center.x()-sqrt2d2, center.y()-sqrt2d2);
//    insert_new_seed(center.x(), center.y()-1.);
//    insert_new_seed(center.x()+sqrt2d2, center.y()-sqrt2d2);

//    insert_new_seed(offset_x, offset_y);
//    insert_new_seed(center.x(), center.y()+grid_side/2.);
//    insert_new_seed(offset_x+grid_side, offset_y);
//    insert_new_seed(center.x()+grid_side/2., center.y());
//    insert_new_seed(offset_x+grid_side, offset_y);
//    insert_new_seed(center.x() + 0.5, center.y());

//    for(int i=0; i<6; ++i)
//      for(int j=0; j<6; ++j)
//        insert_new_seed(center.x()-grid_side/2.+i*grid_side/6.,
//                        center.y()-grid_side/2.+j*grid_side/6.);

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
  typedef boost::array<Grid_point*, 4> Neighbors;

  Point_2 point;
  std::size_t index;
  FT distance_to_closest_seed;
  std::size_t closest_seed_id;
  FMM_state state;
  Neighbors neighbors; // 4 neighbors ORDERED (left, up, right, down)
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

      Vector2d ancestor_edge;
      ancestor_edge(0) = point.x() - curr_ancestor->point.x();
      ancestor_edge(1) = point.y() - curr_ancestor->point.y();
      FT ancestor_edge_length = ancestor_edge.norm();
      Vector2d normalized_anc_edge = ancestor_edge/ancestor_edge_length;

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

        Vector2d curr_edge;
        curr_edge(0) = e0->point.x() - e1->point.x();
        curr_edge(1) = e0->point.y() - e1->point.y();

        // interpolate between both metric and transform the normalized edge
        // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
        Eigen::Matrix2d f = 0.5*(m0.get_transformation() + m1.get_transformation());
        Vector2d transformed_curr_edge = f*normalized_anc_edge;

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
    // consider all the neighbors of a KNOWN point and compute their distance to 'this'
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
        // note that we don't insert in trial_pq since it's already in
        if(gp->compute_closest_seed(this))
          pqs_ret = REBUILD_TRIAL;
      }
      else // gp->state == FAR
      {
        CGAL_assertion(gp->state == FAR);
        if(gp->compute_closest_seed(this))
        {
          gp->state = TRIAL;
          trial_pq.push_back(gp);
          std::push_heap(trial_pq.begin(), trial_pq.end(), Grid_point_comparer<Grid_point>());
        }
      }
    }
    return pqs_ret;
  }

  bool is_contained(const Point_2& p,
                    const Grid_point* next) const
  {
    //    n1
    //    |
    // n0-gp-n2
    //    |
    //    n3

    Grid_point *n0=neighbors[0], *n1=neighbors[1],
               *n2=neighbors[2], *n3=neighbors[3];

    bool is_p_leq_n0_x = (n0) ? p.x() <= n0->point.x() : false;
    bool is_p_geq_n1_y = (n1) ? p.y() >= n1->point.y() : false;
    bool is_p_geq_n2_x = (n2) ? p.x() >= n2->point.x() : false;
    bool is_p_leq_n3_y = (n3) ? p.y() <= n3->point.y() : false;

    if(is_p_leq_n0_x)
      next = n0;
    else if(is_p_geq_n1_y)
      next = n1;
    else if(is_p_geq_n2_x)
      next = n2;
    else if (is_p_leq_n3_y)
      next = n3;
    else
      next = this; // we're inside !
    return (!is_p_leq_n0_x && !is_p_geq_n1_y && !is_p_geq_n2_x && !is_p_leq_n3_y);
  }

  void initialize_from_point(const Point_2& p,
                             const std::size_t seed_id)
  {
    Vector2d v;
    v(0) = p.x() - point.x();
    v(1) = p.y() - point.y();
    const Eigen::Matrix2d& m = metric.get_mat();
    FT d = std::sqrt(v.transpose() * m * v);
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
    ancestor = NULL;
  }

  Grid_point(const Point_2& point_, const std::size_t index_)
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

  // For the approximate geodesic triangulation...
  // the first pair are the IDs of two Voronoi cells (which a non empty inter)
  // the second pair are the two closest points on the bisector (one in each cell)

  // WARNING: currently it's on the part of the bisector that exists, not the whole
  // bisector, which means we have a 'local' geodesic triangulation, not the real one...
  typedef std::map<std::pair<std::size_t, std::size_t>,
                   std::pair<const Grid_point*, const Grid_point*> > Edge_bisec_map;

  std::vector<Grid_point> points;
  Grid_point_vector trial_points;

  // Duality
  mutable Edge_bisec_map edge_bisectors;
  mutable std::set<Edge> dual_edges;
  mutable std::set<Tri> dual_triangles;

  // Refinement
  mutable const Grid_point* refinement_point;
  mutable FT refinement_distance;

  void clear_dual()
  {
    edge_bisectors.clear();
    dual_edges.clear();
    dual_triangles.clear();

    refinement_point = NULL;
    refinement_distance = 0.;
  }

  const Grid_point* locate_point(const Point_2& p) const
  {
    bool found = false;
    const Grid_point& next = points[0];
    while(!found)
      found = (&next)->is_contained(p, &next);
    return &next;
  }

  void locate_and_initialize(const Point_2& p,
                             const std::size_t seed_id)
  {
    // WARNING: only usable for a full grid since we assume the step size is constant
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);

    Grid_point* gp = &(points[index_y*n + index_x]);
    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id << std::endl;
    }

    if(gp->state == TRIAL) // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");

    gp->initialize_from_point(p, seed_id);
    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

#if (verbose > 10)
    std::cout << "looking for p: " << p.x() << " " << p.y() << std::endl;
    std::cout << "found gp: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << std::endl;
#endif
  }

  void initialize_geo_grid()
  {
#if (verbose > 5)
    std::cout << "grid initialization" << std::endl;
#endif

    // create the grid points (either a full grid or a smart grid (todo) )
    for(unsigned int j=0; j<n; ++j)
    {
      for(unsigned int i=0; i<n; ++i)
      {
        Point_2 p(offset_x+i*step, offset_y+j*step); // fill from bot left to top right
        Grid_point gp(p, i + j*n);
        points.push_back(gp);
      }
    }

    // assign the neighbors
    for(unsigned int j=0; j<n; ++j)
    {
      for(unsigned int i=0; i<n; ++i)
      {
        std::size_t curr_id = i + j*n;
        if(i != 0) // there is a neighbor left
          points[curr_id].neighbors[0] = &(points[i-1 + j*n]);
        if(j != n-1) // there is a neighbor above
          points[curr_id].neighbors[1] = &(points[i + (j+1)*n]);
        if(i != n-1) // there is a neighbor right
          points[curr_id].neighbors[2] = &(points[i+1 + j*n]);
        if(j != 0) // there is a neighbor below
          points[curr_id].neighbors[3] = &(points[i + (j-1)*n]);
      }
    }
#if (verbose > 5)
    std::cout << "neighbors assigned" << std::endl;
#endif

    // seeds belong to a quad, find it
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_2& p = seeds[i];
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
    // Check the neighbors of gp for points that have the state KNOWN.
    // In these, consider those that have a closest_seed_id different than gp;
    // those are points in the dual.
    // If we have never seen that dual already, then we have also found which
    // point we will draw to in the geodesic_triangulation

    std::size_t p_id = gp->closest_seed_id;
    std::set<std::size_t> simplex;
    simplex.insert(p_id);

    boost::array<Grid_point*, 4>::const_iterator it = gp->neighbors.begin(),
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
    clear_dual();
    Grid_point* gp;

    bool is_t_empty = trial_points.empty();
    while(!is_t_empty)
    {
      if(known_count%10000 == 0)
        print_states();

#if (verbose > 5)
      std::cout << "Queue sizes. Trial: " << trial_points.size() << std::endl;
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
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << ", " << gp->point.y() << ") ";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from " << gp->closest_seed_id << std::endl;
#endif

      gp->state = KNOWN;

      known_count++; // tmp
      dual_shenanigans(gp);
      boost::array<Grid_point*, 4>::const_iterator it = gp->neighbors.begin(),
                                                   iend = gp->neighbors.end();
      for(; it!=iend; ++it)
      {
        const Grid_point* gq = *it;
        if(gq && gq->state == KNOWN)
          dual_shenanigans(gq);
      }

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
      boost::array<Grid_point*, 4>::const_iterator it = gp->neighbors.begin(),
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

    known_count = 0;
    trial_count = 0;
    far_count = points.size();
  }

  void refine_grid(const Point_2 new_seed)
  {
    vertices_nv = insert_new_seed(new_seed.x(), new_seed.y());
    refresh_grid_point_states(); // we can't spread from the new seed if all states are 'KNOWN'
    locate_and_initialize(new_seed, seeds.size()-1);

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
    // we define a tiny (2k+1)*(2k+1) grid ('sweeper'), and we look at the colors within
    // if we have more than 4 colors at a time, we have a (quasi-) cosphericity

    // this is a bit expensive right now since every point is visited (2*k+1)² times
    // additionally we compute combinations, have to insert in a set, yadayada

    std::set<boost::array<std::size_t, 4> > quasi_cosphericities;

    int h = (std::min)(1., n/100);
    int sweeper_size = 2*h+1;
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
        //const Grid_point& sweeper_center = points[j*n + i];

        // gather the pts
        std::set<std::size_t> colors;
        for(int lx=-h; lx<=h; ++lx)
        {
          for(int ly=-h; ly<=h; ++ly)
          {
            const Grid_point& gp = points[(j+ly)*n + (i+lx)];
            colors.insert(gp.closest_seed_id);
          }
        }

        if(colors.size() > 3)
        {
          // we have at least 4 pts in a quasi cospherical configuration...
          // get 'em all with combinations!

          std::vector<boost::array<std::size_t, 4> > combis =
                             combinations<4>(colors, colors.begin(), 4);

          quasi_cosphericities.insert(combis.begin(), combis.end());
        }
      }
    }

    std::cout << "Sweeper no sweeping ! (with size: " << h << ")" << std::endl;
    std::cout << "Found " << quasi_cosphericities.size() << " (quasi-)cosphericities" << std::endl;

    std::set<boost::array<std::size_t, 4> >::const_iterator sit =
                                                   quasi_cosphericities.begin();
    std::set<boost::array<std::size_t, 4> >::const_iterator send =
                                                     quasi_cosphericities.end();
    for(; sit!=send; ++sit)
    {
      const boost::array<std::size_t, 4>& quasi_cosphericity = *sit;
      for(std::size_t i=0; i<sit->size(); ++i)
        std::cout << quasi_cosphericity[i] << " ";
      std::cout << std::endl;
    }
    std::cout << "End of quasi cosphericities" << std::endl;
  }

  void output_grid(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << points.size() << std::endl;
    out_bb << "2 1 " << points.size() << " 2" << std::endl;

    std::set<Tri> simplices;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << i+1 << std::endl;

//      out_bb << gp.closest_seed_id << std::endl;
      out_bb << gp.distance_to_closest_seed << std::endl;

      const Grid_point* n0 = gp.neighbors[0];
      const Grid_point* n1 = gp.neighbors[1];
      const Grid_point* n2 = gp.neighbors[2];
      const Grid_point* n3 = gp.neighbors[3];

      if(n1 && n2)
      {
        Tri simplex;
        simplex[0] = gp.index;
        simplex[1] = n1->index;
        simplex[2] = n2->index;
        simplices.insert(simplex);
      }
      if(n0 && n3)
      {
        Tri simplex;
        simplex[0] = gp.index;
        simplex[1] = n0->index;
        simplex[2] = n3->index;
        simplices.insert(simplex);
      }
    }

    out << "Triangles" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tri>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tri& t = *it;
      std::set<std::size_t> materials;

      for(std::size_t i=0; i<t.size(); ++i)
      {
        out << t[i] + 1 << " ";
        materials.insert(points[t[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
      //if only one value at all points, we're okay
      // if more than one, print a unique frontier material
    }

    // to draw an approximate geodesic triangulation... tmp
    if(!dual_edges.empty() && !edge_bisectors.empty())
    {
      std::stringstream geo_tr_out;
      output_geodesic_dual(geo_tr_out); // grab the geodesic triangulation's edges
      out << geo_tr_out.str();
    }

    out << "End" << std::endl;
  }

  void add_triangles_edges_to_dual_edges(const Tri& tr) const
  {
    // no need to sort as long as we keep the indices in order
    Edge e;
    e[0] = tr[0]; e[1] = tr[1]; dual_edges.insert(e);
    e[0] = tr[0]; e[1] = tr[2]; dual_edges.insert(e);
    e[0] = tr[1]; e[1] = tr[2]; dual_edges.insert(e);
  }

  template<typename Iterator>
  void add_triangles_edges_to_dual_edges(Iterator it, Iterator end) const
  {
    for(; it!=end; ++it)
    {
      const Tri& tr = *it;
      add_triangles_edges_to_dual_edges(tr);
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

    // Not sure if we should allow dual of edges to be considered...
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
    else // dual is made of > 3 points (that's a cosphericity)
    {
      std::cerr << "WARNING: COSPHERICITY... " << dual_simplex.size() << std::endl;

      // insert all combinations of tetrahedras... (this will create illegal intersections)
      std::vector<boost::array<std::size_t, 3> > combis =
                         combinations<3>(dual_simplex, dual_simplex.begin(), 3);
      dual_triangles.insert(combis.begin(), combis.end());
      add_triangles_edges_to_dual_edges(combis.begin(), combis.end());
    }
  }

  void compute_dual() const
  {
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      Simplex dual_simplex;
      dual_simplex.insert(gp.closest_seed_id);

      boost::array<Grid_point*, 4>::const_iterator it = gp.neighbors.begin(),
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
    // now loop on the map of dual edges...
    // each entry gives us the following : 2 cell ids and 2 points that are
    // on each side of the bisector
    // The dual edge is made of a generator to one point, then the edge between
    // both points, and the second point to the second generator
    // we navigate the ancestors to find our way to the respective generators

    std::stringstream os;
    std::size_t edge_n = 0;

    typename Edge_bisec_map::const_iterator edmit = edge_bisectors.begin();
    typename Edge_bisec_map::const_iterator edmend = edge_bisectors.end();
    for(; edmit!=edmend; ++edmit)
    {
//      const std::size_t p_id = edmit->first.first;
//      const std::size_t q_id = edmit->first.second;
      const Grid_point* gp = edmit->second.first;
      const Grid_point* gq = edmit->second.second;

      edge_n++;
      os << gp->index+1 << " " << gq->index+1 << " 6" << std::endl;

//      std::cout << "p: " << gp->index << " and his tree: " << std::endl;
      const Grid_point* p_anc = gp->ancestor;
      while(p_anc)
      {
//        std::cout << p_anc->index << " ";
        edge_n++;
        os << gp->index+1 << " " << p_anc->index+1 << " 7" << std::endl;
        gp = gp->ancestor;
        p_anc = p_anc->ancestor;
      }
//      std::cout << std::endl;

//      std::cout << "q: " << gq->index << " and his tree: " << std::endl;
      const Grid_point* q_anc = gq->ancestor;
      while(q_anc)
      {
//        std::cout << q_anc->index << " ";
        edge_n++;
        os << gq->index+1 << " " << q_anc->index+1 << " 8" << std::endl;
        gq = gq->ancestor;
        q_anc = q_anc->ancestor;
      }
//      std::cout << std::endl << std::endl;
    }

    out << "Edges" << std::endl;
    out << edge_n << std::endl;
    out << os.str();
  }

  void output_straight_dual(const std::string str_base) const
  {
    if(dual_edges.empty() && dual_triangles.empty())
      compute_dual();

#if (verbose > 5)
    std::cout << "captured: ";
    std::cout << dual_edges.size() << " edges, ";
    std::cout << dual_triangles.size() << " triangles" << std::endl;
#endif

    std::ofstream out((str_base + "_dual.mesh").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << i+1 << std::endl;

    out << "Triangles" << std::endl;
    out << dual_triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = dual_triangles.begin();
                                         it != dual_triangles.end(); ++it)
    {
      const Tri& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
      out << "1" << std::endl; // todo is_intersected like in grid_gen.cpp
    }
    out << "End" << std::endl;
  }

  void output_grid_and_dual(const std::string str_base) const
  {
    output_grid(str_base);
    output_straight_dual(str_base);
  }

  Geo_grid()
    :
      points(),
      trial_points(),
      edge_bisectors(),
      dual_edges(),
      dual_triangles(),
      refinement_point(NULL),
      refinement_distance(0.)
  { }
};

void initialize()
{
  mf = new Euclidean_metric_field<K>(1., 1.);
//  mf = new Custom_metric_field<K>();

  vertices_nv = build_seeds();
  CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );
}

int main(int, char**)
{
  std::cout.precision(17);
  //std::freopen("log.txt", "w", stdout);

  double duration;
  start = std::clock();

  std::srand(0);
  initialize();

  Geo_grid gg;
  gg.build_grid();

  if(n_refine)
  {
    gg.output_grid_and_dual("pre_ref");

    for(int i=0; i<n_refine; ++i)
    {
      gg.refine_grid_with_self_computed_ref_point();
      std::ostringstream out;
      out << "ref_" << i;
      gg.output_grid_and_dual(out.str());
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cerr << "End refinement: " << duration << std::endl;
  }

  gg.output_grid_and_dual("geo_grid_2");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
