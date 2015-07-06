#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <Eigen/Dense>
#include <omp.h>

#include <iostream>
#include <vector>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 3>                         Tri;

typedef typename Eigen::Matrix<double, 2, 1>                 Vector2d;

//typedef typename CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
typedef typename CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;

#define COMPUTE_DUAL
#define FILTER_SEEDS_OUTSIDE_GRID
#define USE_RECURSIVE_UPDATES
#define TMP_REFINEMENT_UGLY_HACK

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 3;

// grid related stuff
Point_2 center(0.,0.);
const FT grid_side = 4.0;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT n = 400.; // number of points per side of the grid
FT step = grid_side / n;

// the metric field and the seeds
MF mf;
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

// how big of an update the new value needs to be
FT recursive_tolerance = 1e-4;

#ifdef TMP_REFINEMENT_UGLY_HACK
FT best_ref_x = 1e30, best_ref_y = 1e30;
#endif

enum FMM_state
{
  CHANGED = 0, // changed implies known. It is used to not add elements that are already in the queue
  KNOWN,
  TRIAL,
  FAR
};

template<typename Gp>
struct Grid_point_comparer
{
  bool operator()(Gp const * const gp1, Gp const * const gp2)
  {
    return gp1->distance_to_closest_seed < gp2->distance_to_closest_seed;
  }
};

int insert_new_seed(const FT x, const FT y)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. )
    {
      std::cout << "filtered : " << x << " " << y << std::endl;
      return seeds.size();
    }
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

    if(seeds.size() == vertices_nv)
      break;
  }
  std::cout << "seeds: " << seeds.size() << std::endl;
  return seeds.size();
}

struct Grid_point
{
  typedef boost::array<Grid_point*, 4> Neighbors;

  Point_2 point;
  std::size_t index;
  mutable FT distance_to_closest_seed;
  mutable std::size_t closest_seed_id;
  mutable  FMM_state state;
  Neighbors neighbors; // 4 neighbors ORDERED (left, up, right, down)
                       // if no neighbor (borders for example) then NULL
  Metric metric;

  const Grid_point* ancestor; // used to debug & draw nice things

  bool compute_closest_seed_1D(const Grid_point* gp)
  {
//    std::cout << "compute closest 1D " << index << " to " << gp->index << std::endl;
//    std::cout << "gp distance to closest is : " << gp->distance_to_closest_seed << std::endl;

    CGAL_assertion(gp->state == KNOWN || gp->state == CHANGED);
    bool changed = false;

    Vector2d v;
    v(0) = gp->point.x() - point.x();
    v(1) = gp->point.y() - point.y();
    FT neighbor_d = std::sqrt(v.transpose()*metric.get_mat()*v);
    FT dcs_at_gp = gp->distance_to_closest_seed;
    FT d = dcs_at_gp + neighbor_d;
    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
//      std::cout << "new closest at : " << distance_to_closest_seed << std::endl;
    }
    return changed;
  }

  FT eval_at_p(FT p, const FT a, const FT b, const FT c, const FT d, const FT e) const
  {
    if(std::abs(c*p*p + d*p + e) < 1e-16) // avoids stuff like -1e-17
      return p*a+(1.-p)*b;
    CGAL_assertion(c*p*p + d*p + e >= 0);
    return p*a+(1.-p)*b + std::sqrt(c*p*p + d*p + e);
  }

  FT compute_min_distance_2D(const Grid_point* gp, const Grid_point* gq,
                             FT& p) const
  {
    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );

    // in 2D, we can solve the minimization analytically.
    // We're solving :
    // min_{p\in [0,1]} f(p) = F(x)*p + F(y)*(1-p) + sqrt(v(p)^t * M * v(p))
    // vp = p*gp+(1-p)*gq
    // local minimum is given by f'(p_0) = 0 _but_ we have to compare f(p_0)
    // with f(0) and f(1) as we're only interested in solutions that live in [0,1]
    // If there's no solution in [0,1], then there's no (local) minimum in [0,1]
    // and the minimum is obtained in 0 or in 1.

    const FT dcs_at_gp = gp->distance_to_closest_seed;
    const FT dcs_at_gq = gq->distance_to_closest_seed;
    Vector2d vp, vq;
    vp(0) = gp->point.x() - point.x();
    vp(1) = gp->point.y() - point.y();
    vq(0) = gq->point.x() - point.x();
    vq(1) = gq->point.y() - point.y();

    Eigen::Matrix2d m = metric.get_mat();
    const FT vpmvp = vp.transpose() * m * vp;
    const FT vpmvq = vp.transpose() * m * vq;
    const FT vqmvq = vq.transpose() * m * vq;

    // expressing f(p) as ap+(1-p)b+sqrt(cp^2+dp+e)
    // the tiny trick is to write f'(p) = 0 and then square the expression to get rid of sqrt

    const FT a = dcs_at_gp;
    const FT b = dcs_at_gq;
    const FT c = vpmvp + vqmvq - 2.*vpmvq;
    const FT d = 2.*(vpmvq - vqmvq);
    const FT e = vqmvq;

    const FT A = 4*(a-b)*(a-b)*c-4*c*c;
    const FT B = 4*(a-b)*(a-b)*d-4*c*d;
    const FT C = 4*(a-b)*(a-b)*e-d*d;
    FT delta = B*B-4*A*C;

//    std::cout << "a: " << a << std::endl;
//    std::cout << "b: " << b << std::endl;
//    std::cout << "c: " << c << std::endl;
//    std::cout << "d: " << d << std::endl;
//    std::cout << "e: " << e << std::endl;
//    std::cout << "A: " << A << std::endl;
//    std::cout << "B: " << B << std::endl;
//    std::cout << "C: " << C << std::endl;
//    std::cout << "delta: " << delta << std::endl;

    FT min_id = 0.;
    FT min = eval_at_p(0., a, b, c, d, e);

    FT d1 = eval_at_p(1., a, b, c, d, e);
    if(min > d1)
    {
      min_id = 1.;
      min = d1;
    }

    // get solutions
    CGAL_assertion(A!=0. || B!=0.);
    if(A == 0.)
    {
      p = -C/B;
      FT dp = eval_at_p(p, a, b, c, d, e);
      if(min > d && p >= 0 && p <= 1.)
      {
        min_id = p;
        min = dp;
      }
    }
    else if(delta >= 0)
    {
      FT sqrt_delta = std::sqrt(delta);
      FT denom = 1. / (2.*A);
      double p1 = (-B + sqrt_delta) * denom;
      double p2 = (-B - sqrt_delta) * denom;

      FT dp1 = eval_at_p(p1, a, b, c, d, e);
      FT dp2 = eval_at_p(p2, a, b, c, d, e);

//      std::cout << "p1 : " << p1 << " dp1: " << dp1 << std::endl;
//      std::cout << "p2 : " << p2 << " dp2: " << dp2 << std::endl;

      if(min > dp1 && p1 >= 0 && p1 <= 1.)
      {
        min_id = p1;
        min = dp1;
      }
      if(min > dp2 && p2 >= 0 && p2 <= 1.)
      {
        min_id = p2;
        min = dp2;
      }
    }

    p = min_id;
//    std::cout << "final min at : " << p << " d: " << min << std::endl;
    return min;
  }

  bool compute_closest_seed_2D(const Grid_point* gp, const Grid_point* gq)
  {
//    std::cout << "compute closest 2D " << index << " to " << gp->index << " & " << gq->index << std::endl;

    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );
    bool changed = false;

    FT p = 0.;
    const FT d = compute_min_distance_2D(gp, gq, p);

    // tolerance to ignore unsignificant changes
    if(distance_to_closest_seed - d > recursive_tolerance *d)
    {
      changed = true;
      distance_to_closest_seed = d;
      FT gp_d = gp->distance_to_closest_seed;
      FT gq_d = gq->distance_to_closest_seed;
      const Grid_point* best_p = (gp_d < gq_d) ? gp : gq; // fixme the commented line below should be the correct one...
//      const Grid_point* best_p = (p > 0.5) ? gp : gq;
      closest_seed_id = best_p->closest_seed_id;
      ancestor = best_p; // a bit ugly, create the point p ?
    }
    return changed;
  }

  bool compute_closest_seed(const Grid_point* ancestor,
                            const bool ignore_ancestor = false)
  {
//    std::cout << "compute closest high level " << index << std::endl;

    bool changed = false;
    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      const Grid_point* gp = neighbors[i];
      const Grid_point* gq = neighbors[(i+1)%(neighbors.size())];

      if(gp && (gp->state == KNOWN || gp->state == CHANGED))
      {
        if(gq && (gq->state == KNOWN || gq->state == CHANGED))
        {
          if(ignore_ancestor || gp == ancestor || gq == ancestor)
            if(compute_closest_seed_2D(gp, gq))
              changed = true;
        }
        else // we could do it++ here since we know the next 'it' has gp->state != KNOWN
        {
          if(ignore_ancestor || gp == ancestor)
            if(compute_closest_seed_1D(gp))
              changed = true;
        }
      }
    }
    return changed;
  }

  template<typename PQ>
  void update_neighbors_distances(PQ& changed_pq, PQ& trial_pq) const
  {
//    std::cout << "update neigh" << std::endl;

    // consider all the known neighbors and compute the value to that
    CGAL_assertion(state == KNOWN || state == CHANGED);

    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      Grid_point* gp = neighbors[i]; // a neighbor of 'this'
      if(!gp)
        continue;

//      std::cout << "gp: " << gp->index << " is go, status: " << gp->state << std::endl;

#ifdef USE_RECURSIVE_UPDATES
      if(gp->state == KNOWN || gp->state == CHANGED) // that's the recursive part
      {
        // recompute the distance in the direction gp-'this'
        FT mem = gp->distance_to_closest_seed;
        std::size_t mem_ancestor = (gp->ancestor)?gp->ancestor->index:-1;
        if(gp->compute_closest_seed(this))
        {
          std::cout.precision(17);
          std::cout << "new value at " << gp->index << " : " << gp->distance_to_closest_seed;
          std::cout << " with ancestor: " << gp->ancestor->index;
          std::cout << " (mem: " << mem << " ancestor: " << mem_ancestor << ") ";
          std::cout << " diff: " << gp->distance_to_closest_seed-mem << std::endl;
          if(gp->state == KNOWN)
          {
            // if it's already changed we only need to reorder the changed queue
            std::cout << "and inserting it" << std::endl;
            gp->state = CHANGED;
            changed_pq.push_back(gp);
            std::push_heap(changed_pq.begin(), changed_pq.end(), Grid_point_comparer<Grid_point>());
          }
        }
      }
      else
#endif
      if(gp->state == TRIAL)
      {
        gp->compute_closest_seed(this); // need to change the order of the element
        // note that we don't insert in trial_pq since it's already in
        // the sort after this function will change its position if needs be
      }
      else if(gp->state == FAR)
      {
        gp->state = TRIAL;
        gp->compute_closest_seed(this); // always returns true here so no need to test it
        trial_pq.push_back(gp);
        std::push_heap(trial_pq.begin(), trial_pq.end(), Grid_point_comparer<Grid_point>());
      }
    }
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
                             const std::size_t seed_id) const
  {
    Vector2d v;
    v(0) = p.x() - point.x();
    v(1) = p.y() - point.y();
    FT d = std::sqrt(v.transpose()*metric.get_mat()*v);
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
  }

  Grid_point(const Point_2& point_, const std::size_t index_)
    : point(point_), index(index_), distance_to_closest_seed(1e30),
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

  std::vector<Grid_point> points;
  Grid_point_vector changed_points;
  Grid_point_vector trial_points;

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
    // only usable for a full grid
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);

    Grid_point* gp = &(points[index_y*n + index_x]);
    if(gp->closest_seed_id != -1)
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id << std::endl;
    }

    // if gp is found the same for 2 different initial points, then your grid
    // is not dense enough...

    gp->initialize_from_point(p, seed_id);
    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    std::cout << "looking for p: " << p.x() << " " << p.y() << std::endl;
    std::cout << "found: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << "]" << std::endl;

    // if we want to initialize a bit more, we initialize all the vertices of the quad
    // that contains p
    // is this really useful... ?

/* //fixme these neighbors don't always exist
    CGAL_assertion(gp->neighbors[1] && gp->neighbors[2] && (gp->neighbors[2])->neighbors[1]);
    Grid_point* gq = gp->neighbors[2];
    Grid_point* gr = gq->neighbors[1];
    Grid_point* gs = gp->neighbors[1]; // which is also gr->nei[0]

    gq->initialize_from_point(p, seed_id);
    trial_points.push_back(gq);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    gr->initialize_from_point(p, seed_id);
    trial_points.push_back(gr);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    gs->initialize_from_point(p, seed_id);
    trial_points.push_back(gs);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    std::cout << "neighbors are : " << gq->index << " " << gq->point.x() << " " << gq->point.y() << std::endl;
    std::cout << "neighbors are : " << gr->index << " " << gr->point.x() << " " << gr->point.y() << std::endl;
    std::cout << "neighbors are : " << gs->index << " " << gs->point.x() << " " << gs->point.y() << std::endl;
*/
    std::sort_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
  }

  void initialize_geo_grid()
  {
    std::cout << "grid ini" << std::endl;

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

    // assign the neighbors (todo: find a better way...)
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
    std::cout << "assigned neighbors" << std::endl;

    // seeds belong to a quad, find it and initialize the pts of the quads
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_2& p = seeds[i];
//      const Grid_point* gp = locate_point(p);
//      gp->initialize_from_point(p, i);
      locate_and_initialize(p, i);
    }

    std::cout << "grid initialization end" << std::endl;
  }

  bool print_states() const
  {
    int changed_count=0, known_count=0, trial_count=0, far_count=0;

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      if(gp.state == CHANGED)
        changed_count++;
      if(gp.state == KNOWN)
        known_count++;
      if(gp.state == TRIAL)
        trial_count++;
      if(gp.state == FAR)
        far_count++;
    }
    std::cout << "changed: " << changed_count << " known: " << known_count;
    std::cout << " trial: " << trial_count << " far: " << far_count << std::endl;

//    if(known_count > 35000)
//      return true;
    return false;
  }

  void geo_grid_loop()
  {
    std::cout << "main loop" << std::endl;

    bool is_cp_empty = changed_points.empty();
    bool is_t_empty = trial_points.empty();

    Grid_point* gp;
    while(!is_cp_empty || !is_t_empty)
    {
//      std::cout << "changed heap: " << std::endl;
//      for (std::vector<Grid_point*>::iterator it = changed_points.begin();
//           it != changed_points.end(); ++it)
//        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
//      std::cout << std::endl;

//      std::cout << "trial heap: " << std::endl;
//      for (std::vector<Grid_point*>::iterator it = trial_points.begin();
//           it != trial_points.end(); ++it)
//        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
//      std::cout << std::endl;

      if(!is_cp_empty)
      {
        gp = changed_points.front();
        std::pop_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());
        changed_points.pop_back();
      }
      else // !is_t_empty
      {
        gp = trial_points.front();
        std::pop_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
        trial_points.pop_back();
      }

      if(print_states())
        break;
      std::cout << "Queue sizes. Trial: " << trial_points.size() << " Changed: " << changed_points.size() << std::endl;
      std::cout << "picked: " << gp->index << " | " << gp->point.x() << " " << gp->point.y() << std::endl;

      gp->state = KNOWN;
      gp->update_neighbors_distances(changed_points, trial_points);
      std::sort_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
      std::sort_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());

      is_cp_empty = changed_points.empty();
      is_t_empty = trial_points.empty();
    }

#ifdef USE_RECURSIVE_UPDATES
    std::cout << "DOUBLE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point* gp = &(points[i]);
      CGAL_assertion(!(gp->compute_closest_seed(gp, true))); // if 'true', we forgot a move
    }
#endif
  }

  void build_grid()
  {
    initialize_geo_grid();
    geo_grid_loop();
  }

  Point_2 compute_refinement_point()
  {
    // todo
  }

  void refresh_grid_after_new_seed_creation()
  {
    std::size_t seed_id = seeds.size() - 1;
    const Point_2& p = seeds.back();
    locate_and_initialize(p, seed_id);
    geo_grid_loop();
  }

  void refine_grid()
  {
#ifdef TMP_REFINEMENT_UGLY_HACK
    Point_2 p(best_ref_x, best_ref_y);
    CGAL_assertion(best_ref_x != 1e30 && best_ref_y != 1e30);
#else
    Point_2 p = compute_refinement_point();
#endif
    vertices_nv = insert_new_seed(p.x(), p.y());
    refresh_grid_after_new_seed_creation();

#ifdef TMP_REFINEMENT_UGLY_HACK
    best_ref_x = 1e30; best_ref_y = 1e30;
#endif
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

//      CGAL_assertion((!n0 || n0->point.x() == gp.point.x() - step) &&
//                     (!n1 || n1->point.y() == gp.point.y() + step) &&
//                     (!n2 || n2->point.x() == gp.point.x() + step) &&
//                     (!n3 || n3->point.y() == gp.point.y() - step) );

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

//      // we want to draw a botleft-topright diagonal so we only draw
//      // gp-n0-n1 and gp-n2-n3 (if possible obviously)
//      typename Neighbors::const_iterator it = gp.neighbors.begin(),
//                                         iend = gp.neighbors.end();
//      for(; it!=iend; ++it)
//      {
//        const Grid_point* gq = *it;
//        const Grid_point* gr = *((++it));
//        if(gq && gr)
//        {
//          Simplex simplex;
//          simplex.insert(gp.index);
//          simplex.insert(gq->index);
//          simplex.insert(gr->index);
//          simplices.insert(simplex);
//        }
//      }
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
        out << t[i]+1 << " ";
        materials.insert(points[t[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
      //if only one value at all points, we're okay. if more than one, print a unique frontier material
    }
  }

  std::set<Simplex> compute_dual() const
  {
    std::cout << "dual computations" << std::endl;

#ifdef TMP_REFINEMENT_UGLY_HACK
    FT farthest_dual_distance = 0.;
#endif
    std::set<Simplex> dual_simplices;
#pragma omp parallel for
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

      // no one cares if dual_simplex.size() == 1
      if(dual_simplex.size() == 2)
      {
        // that'd be an edge in the dual if we cared about it
      }
      else if(dual_simplex.size() >= 3)
      {
        if(dual_simplex.size() == 3)
#pragma omp critical
          dual_simplices.insert(dual_simplex); // a triangle!
        else
          std::cerr << "WARNING: COSPHERICITY..." << std::endl;
#ifdef TMP_REFINEMENT_UGLY_HACK
        if(gp.distance_to_closest_seed > farthest_dual_distance)
        {
#pragma omp critical
{
          farthest_dual_distance = gp.distance_to_closest_seed;
          best_ref_x = gp.point.x();
          best_ref_y = gp.point.y();
}
        }
#endif
      }
    }

//    // a dual exists if a simplex has 3 different values
//    for(typename std::set<Simplex>::iterator it = simplices.begin();
//                                             it != simplices.end(); ++it)
//    {
//      std::set<std::size_t> closest_seeds;
//      const Simplex& s = *it;
//      typename Simplex::iterator sit = s.begin();
//      typename Simplex::iterator siend = s.end();
//      for(; sit!=siend; ++sit)
//        closest_seeds.insert(points[*sit].closest_seed_id);
//      if(closest_seeds.size() == 3)
//      {
//        Simplex simplex;
//        for(sit=s.begin(); sit!=siend; ++sit)
//        {

//          simplex.insert(points[*sit].closest_seed_id);

//        }
//        dual_simplices.insert(simplex);
//      }
//    }

#ifdef TMP_REFINEMENT_UGLY_HACK
    if(dual_simplices.empty())
      std::cerr << "WARNING: no simplices captured, so no best_ref_x,y set" << std::endl;
#endif

    return dual_simplices;
  }

  void output_dual() const
  {
    const std::set<Simplex>& dual_simplices = compute_dual();
    std::cout << "captured: " << dual_simplices.size() << " simplices" << std::endl;

    std::ofstream out("geo_grid_dual.mesh");
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(int i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << i+1 << std::endl;

    out << "Triangles" << std::endl;
    out << dual_simplices.size() << std::endl;
    for(typename std::set<Simplex>::iterator it = dual_simplices.begin();
        it != dual_simplices.end(); ++it)
    {
      const Simplex& s = *it;
      typename Simplex::iterator sit = s.begin();
      typename Simplex::iterator siend = s.end();
      for(; sit!=siend; ++sit)
        out << *sit+1 << " ";
      out << "1" << std::endl;
    }
    out << "End" << std::endl;
  }

  Geo_grid() : points(), changed_points(), trial_points()
  {
    std::make_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());
    std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
  }
};

void initialize()
{
//  mf = new Euclidean_metric_field<K>(10,1);
  mf = new Custom_metric_field<K>();

  vertices_nv = build_seeds();
}

int main(int, char**)
{
  std::clock_t start;
  double duration;
  start = std::clock();

  std::freopen("geo_grid_log.txt", "w", stdout);
  std::srand(0);
  initialize();

  Geo_grid gg;
  gg.build_grid();
  gg.output_grid("geo_grid");

  int n_refine = 50;
  for(int i=0; i<n_refine; ++i)
  {
    gg.refine_grid();
    std::ostringstream out;
    out << "geo_grid_ref_" << i;
    gg.output_grid(out.str());
  }

  gg.output_grid("geo_grid");
  gg.output_dual();

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "dur: " << duration << std::endl;
}
