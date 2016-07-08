// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// canvas is a starset
#define ANISO_NO_CONSISTENCY

#include <CGAL/Anisotropic_mesher_2.h>
#include <CGAL/Starset.h>
#include <CGAL/IO/Star_set_output.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/External_metric_field.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

#define ANISO_USE_RECIPROCAL_NEIGHBORS
// #define COMPUTE_PRECISE_VOR_VERTICES
// #define USE_FULL_REBUILDS
#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 2

// #define COMPUTE_GEODESIC
#define COMPUTE_REAL_GEODESICS
// #define OVERWRITE_SEED_WITH_CANVAS_POINT // use that for easier geodesic drawing

namespace CGAL
{
namespace Anisotropic_mesh_2
{
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Exact_predicates_exact_constructions_kernel     KExact;

typedef K::FT                                                 FT;

typedef CGAL::Anisotropic_mesh_2::Metric_base<K>              Metric;
typedef CGAL::Anisotropic_mesh_2::Starset<K>                  Star_set;

typedef typename Star_set::Star_handle                        Star_handle;
typedef typename Star_set::Vertex_handle                      Vertex_handle;
typedef typename Star_set::Vertex_handle_handle               Vertex_handle_handle;
typedef typename Star_set::Face_handle                        Face_handle;
typedef typename Star_set::Face_handle_handle                 Face_handle_handle;

typedef K::Point_2                                            Point_2;
typedef K::Point_3                                            Point_3;
typedef std::set<std::size_t>                                 Simplex;
typedef boost::array<std::size_t, 2>                          Edge;
typedef boost::array<std::size_t, 3>                          Tri;
typedef Eigen::Matrix<double, 2, 1>                           Vector2d;

typedef K::Segment_2                                          Segment;
typedef K::Vector_2                                           Vector_2;
typedef K::Triangle_2                                         Triangle_2;
typedef K::Segment_3                                          Segment_3;
typedef K::Triangle_3                                         Triangle_3;

typedef KExact::Point_2                                       EPoint;
typedef KExact::Segment_2                                     ESegment;
typedef KExact::Triangle_2                                    ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                  To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                  Back_from_exact;

// using a lot of global variables, it's ugly
// but passing them by function is tedious

const FT FT_inf = std::numeric_limits<FT>::infinity();
To_exact to_exact;
Back_from_exact back_from_exact;

//typedef CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
typedef CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;
//typedefCGAL::Anisotropic_mesh_2::External_metric_field<K>* MF;
MF mf;

//geometry
FT a = 2.1; // half the side (x)
FT b = 2.1; // half the side (y)
Point_2 center(0., 0.);
Rectangle_domain<K> pdomain(a, b, center);

// seeds

enum Seed_status
{
  IN_DOMAIN = 0,
  ON_CONSTRAINT,
  ON_CORNER
};

FT sa = 2.; // seeds bbox
FT sb = 2.;
std::size_t vertices_nv = 1e10;
bool are_constraints_discretized = false;

std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;
std::vector<Seed_status> seeds_s;
const std::string str_seeds = "swirl_input_13641";

// base mesh
const std::string str_base_mesh = "swirl";
CGAL::Bbox_2 base_mesh_bbox;

// refinement
int n_refine = 1359;
std::size_t min_ancestor_path_length = 5;
static const int k = 8; // depth of the ancestor edge

// optimization
int max_opti_n = 0;
int max_depth = -1; // how precise are the Voronoi vertices (bigger = better)

//debug & info
int known_count=0, trial_count=0, far_count=0;
std::clock_t start;
std::vector<int> optimal_depths(k, 0);

typedef boost::unordered_map<boost::array<std::size_t, k+1>, FT> Path_map;
Path_map computed_paths;

// stuff to generate the grid using Mesh_2 since Aniso_mesh_2 is a turtle.
// Need to define our own refinement criterion based on metric shenanigans
bool ignore_children = false;

//face criteria
FT f_r0 = .1; // since geodesic takes r0=1, we need r_0 << 1 here (0.1-ish ?)
FT f_rho0 = 3.0;

//misc
FT gamma = 1.0;
FT beta = 2.5;
FT delta = 0.3;
int max_times_to_try = 60;
int nb = 20; // nb_initial_points

// max number of stars
std::size_t max_stars_n = -1;

Criteria_base<K> criteria(f_r0, f_rho0, gamma, beta, delta, nb, max_times_to_try);

void generate_starset(Star_set& ss)
{
  Anisotropic_mesher_2<K> mesher(ss, &pdomain, &criteria, mf);
  mesher.refine_mesh();

  std::ofstream out("geo_starset.mesh");
  output_medit(ss, out, false/*consistent_only*/);

  // useful to quickly recreate a starset
  std::ofstream out_dump((str_base_mesh + "_starset.dump").c_str());
  dump(ss, out_dump);
}

void read_dump(Star_set& ss,
               std::string filename)
{
  std::cout << "Reading dump..." << std::endl;
  std::ifstream in(filename.c_str());
  ss.clear();

  std::size_t stars_n, v_n, id;
  FT x, y;

  in >> stars_n;

  stars_n = (std::min)(max_stars_n, stars_n);

  for(std::size_t i=0; i<stars_n; ++i)
  {
    in >> x >> y;
    Point_2 p(x,y);

    Star_handle star = new typename Star_set::Star(&criteria, &pdomain);
    const Metric& m_p = mf->compute_metric(p);
    star->reset(p, i, m_p);
    ss.push_back(star);
  }

  for(std::size_t i=0; i<stars_n; ++i)
  {
    in >> v_n;
    Star_handle star_i = ss[i];
    for(std::size_t j=0; j<v_n; ++j)
    {
      in >> id;
      Star_handle star_j = ss[id];
      star_i->insert_to_star(star_j->center_point(), id, false);
    }
    std::cout << "built star: " << i << std::endl;
  }
  std::cout << ss.size() << " stars from dump" << std::endl;
}

void read_dump(Star_set& ss)
{
  std::string filename = str_base_mesh + "_starset.dump";
  return read_dump(ss, filename);
}

// the functions below aren't very pretty
bool is_face_outside_domain(Star_handle star, Face_handle fh)
{
  Point_2 p0 = star->metric().inverse_transform(fh->vertex(0)->point());
  Point_2 p1 = star->metric().inverse_transform(fh->vertex(1)->point());
  Point_2 p2 = star->metric().inverse_transform(fh->vertex(2)->point());
  FT third = 1./3.;
  return pdomain.side_of_constraint(CGAL::barycenter(p0, third,
                                                     p1, third,
                                                     p2, third)) == CGAL::ON_NEGATIVE_SIDE;
}

bool is_a_corner_vertex(std::size_t vi)
{
  return vi < 4;
}

bool is_a_corner_vertex(Vertex_handle vh)
{
  return vh->info() < 4;
}

bool has_a_corner_vertex(Star_handle star)
{
  Vertex_handle_handle vit = star->finite_adjacent_vertices_begin();
  Vertex_handle_handle vend = star->finite_adjacent_vertices_end();
  for(; vit != vend; vit++)
  {
    Vertex_handle vh = *vit;
    if(is_a_corner_vertex(vh))
      return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
// Geodesic stuff now !
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Utility functions below
// -----------------------------------------------------------------------------

void print_profiler_info()
{
  std::cout << "optimal depth with " << k << " possible size" << std::endl;
  for(int i=0; i<k; ++i)
    std::cout << "i: " << i << " optimal depths: " << optimal_depths[i] << std::endl;

  std::cout << computed_paths.size() << " computed paths" << std::endl;
}

inline std::ostream& operator<<(std::ostream& os, const Simplex& s)
{
  Simplex::iterator it = s.begin();
  for(; it!=s.end(); ++it)
    os << *it << " ";
  os << std::endl;

  return os;
}

std::size_t fact(std::size_t n)
{
  return (n==0)?1:n*fact(n-1);
}

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
      for(std::size_t i=0, llcs=lower_level_combis.size(); i<llcs; ++i)
        lower_level_combis[i][k_max-k] = *(sit);

      combis.insert(combis.end(), lower_level_combis.begin(),
                                  lower_level_combis.end());
    }
  }
  return combis;
}

bool is_triangle_intersected(const Tri& tri,
                             const boost::unordered_set<Edge>& edges)
{
  bool is_intersected = false;

  // need to switch to epeck for the correct intersections here...
  const std::size_t i0 = tri[0];
  const std::size_t i1 = tri[1];
  const std::size_t i2 = tri[2];
  const Point_2& p0 = seeds[i0];
  const Point_2& p1 = seeds[i1];
  const Point_2& p2 = seeds[i2];
  const Triangle_2 triangle(p0, p1, p2);

#pragma omp parallel shared(is_intersected, edges, p0, p1, p2)
  {
    for(boost::unordered_set<Edge>::const_iterator it=edges.begin();
                                                   it!=edges.end(); ++it)
    {
#pragma omp single nowait // hack to parallelize the std::set loop
      {
#pragma omp flush (is_intersected)
        if(!is_intersected) // hack because we're not allowed to 'break' (nor to 'continue')
          // inside a pragma omp for
        {
          const Edge& edge = *it;
          bool local_intersected = false;

          const Segment segment(seeds[edge[0]], seeds[edge[1]]);

          ESegment esegment = to_exact(segment);
          ETriangle etriangle = to_exact(triangle);
          CGAL::cpp11::result_of<KExact::Intersect_2(ESegment, ETriangle)>::type
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
          {      if(local_intersected)
              is_intersected = true;
          }
#pragma omp flush (is_intersected)
        } // end !is_intersected
      } // end single nowait
    } // end for
  } // end parallel
  return is_intersected;
}

FT sq_bbox_diagonal_length(const CGAL::Bbox_2& bbox)
{
  FT dx = bbox.xmax() - bbox.xmin();
  FT dy = bbox.ymax() - bbox.ymin();

  return dx*dx + dy*dy;
}

Eigen::Matrix2d get_interpolated_transformation(const Metric& m0, const Metric& m1)
{
  // note that we return the transformation and not the full metric !
  return 0.5*(m0.get_transformation() + m1.get_transformation());
}

FT quality(const Point_2& p, const Point_2& q, const Point_2& r)
{
  K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();

  FT alpha = 4.*std::sqrt(3.);
  FT A = std::abs(K::Compute_area_2()(p, q, r));
  FT a = sqd(p, q);
  FT b = sqd(p, r);
  FT c = sqd(q, r);
  FT quality = alpha*A/(a+b+c);

  return quality;
}

FT compute_quality(const Tri& tr)
{
  const Point_2& p0 = seeds[tr[0]];
  const Point_2& p1 = seeds[tr[1]];
  const Point_2& p2 = seeds[tr[2]];
  const Metric& m0 = seeds_m[tr[0]];
  const Metric& m1 = seeds_m[tr[1]];
  const Metric& m2 = seeds_m[tr[2]];

  const Point_2 tp0_0 = m0.transform(p0);
  const Point_2 tp1_0 = m0.transform(p1);
  const Point_2 tp2_0 = m0.transform(p2);

  const Point_2 tp0_1 = m1.transform(p0);
  const Point_2 tp1_1 = m1.transform(p1);
  const Point_2 tp2_1 = m1.transform(p2);

  const Point_2 tp0_2 = m2.transform(p0);
  const Point_2 tp1_2 = m2.transform(p1);
  const Point_2 tp2_2 = m2.transform(p2);

  FT qual = (std::min)((std::min)(quality(tp0_0, tp1_0, tp2_0),
                                  quality(tp0_1, tp1_1, tp2_1)),
                       quality(tp0_2, tp1_2, tp2_2));
  return qual;
}

// -----------------------------------------------------------------------------
// Geodesic classes
// -----------------------------------------------------------------------------

struct Voronoi_constraint
{
  Point_2 p1, p2;
  Vector2d dir;
  FT sq_length; // squared but it doesn't matter, only used for comparison

  Voronoi_constraint() { }
  Voronoi_constraint(const Point_2& p1_, const Point_2& p2_)
    :
      p1(p1_),
      p2(p2_),
      dir()
  {
    sq_length = CGAL::squared_distance(p1, p2);
    dir(0) = p2.x() - p1.x();
    dir(1) = p2.y() - p1.y();
    dir /= dir.norm();
  }

  bool has_point(const Point_2& p) const
  {
    typename K::Collinear_are_ordered_along_line_2  collinear_are_ordered_along_line;
    typename K::Orientation_2                       orientation;

    typename K::Orientation o = orientation(p1, p, p2);

    return (o == COLLINEAR &&
            collinear_are_ordered_along_line(p1, p, p2));
  }
};

enum FMM_state
{
  KNOWN = 0,
  TRIAL,
  FAR,
  ORPHAN
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL
};

struct Base_mesh;

struct Canvas_point
{
  typedef boost::unordered_set<std::size_t>                           Point_set;
  typedef std::list<std::size_t>                                      Point_list;

  // immuable stuff
  Base_mesh* bm;
  Point_2 point; // both this and the metric info are redundant with the starset info
  std::size_t index; //  return bm->ss[index]->metric();
  Metric metric; // but it's needed for virtual points... fixme, I suppose.
  bool is_on_domain_border;

#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
  Point_list reciprocal_neighbors;
#endif

  // stuff that depends on the seeds
  FMM_state state;
  int is_seed_holder;
  FT distance_to_closest_seed;
  std::size_t closest_seed_id;
  std::size_t ancestor;
  std::size_t depth;
  Point_set children;
  bool is_Voronoi_vertex;
  std::size_t ancestor_path_length;

  void change_state(FMM_state new_state);
  std::size_t anc_path_length() const;
  void print_children() const;
  void print_ancestors() const;
  bool find_in_ancestry(const std::size_t i) const;
  void remove_from_children(const std::size_t c);
  void mark_descendants(std::size_t& count);
  void reset_descendants();
  void deal_with_descendants();
  bool compute_closest_seed(const std::size_t n_anc, const bool overwrite = true);
  void update_neighbor_distance(Canvas_point& cp,
                                std::vector<std::size_t> & trial_pq,
                                PQ_state& pqs_ret) const;
  PQ_state update_neighbors_distances(std::vector<std::size_t> & trial_pq) const;
  FT distortion_to_seed() const;
  void reset();
  void initialize_from_point(const FT d, const std::size_t seed_id);

  bool operator==(const Canvas_point& cp) const;

  Canvas_point();
  Canvas_point(Base_mesh* bm_,
              const std::size_t index_,
              const bool is_on_domain_border_ = false);
  Canvas_point(Base_mesh * bm_,
               const Point_2& p,
               const std::size_t index_,
               const Metric& metric_,
               const bool is_on_domain_border_ = false);
  Canvas_point(const Canvas_point& cp);
};

template<typename BM>
struct Canvas_point_comparer
{
  BM const * const bm;

  bool operator()(const std::size_t i1, const std::size_t i2)
  {
    const Canvas_point& cp1 = bm->points[i1];
    const Canvas_point& cp2 = bm->points[i2];
    return cp1.distance_to_closest_seed > cp2.distance_to_closest_seed;
  }

  Canvas_point_comparer(BM const * const bm_) : bm(bm_) { }
};

inline std::ostream& operator<<(std::ostream& os,
                                const boost::tuple<Simplex, std::size_t, FT>& pqe)
{
  std::cout << "PQ entry element:" << std::endl;
  std::cout << "Simplex: " << pqe.get<0>();
  std::cout << "cp: " << pqe.get<1>();
  std::cout << "val: " << pqe.get<2>() << std::endl;

  return os;
}

template<typename PQ_entry, typename Base_mesh>
struct PQ_entry_comparator
{
  Base_mesh* bm;

  bool operator()(const PQ_entry& left, const PQ_entry& right) const
  {
    const Canvas_point& cpl = bm->points[left.get<1>()];
    const Canvas_point& cpr = bm->points[right.get<1>()];

    if(left.get<2>() == right.get<2>()) // that's the value
    {
      if(cpl.distance_to_closest_seed ==
         cpr.distance_to_closest_seed) // that's the dist to the ref point
      {
        const Simplex& s1 = left.get<0>();
        const Simplex& s2 = right.get<0>();

        if(s1.size() != s2.size())
          return s1.size() > s2.size(); // values & dist equal, prioritize the biggest simplex
        else
          return s1 > s2;
      }
      else
        return left.get<1>() > right.get<1>(); // first the farthest ref point if values are the same
    }

    return left.get<2>() > right.get<2>(); // first the biggest value
  }

  PQ_entry_comparator(Base_mesh* bm_) : bm(bm_) { }
};

struct Base_mesh
{
private:
  typedef Base_mesh                                         Self;

public:
  typedef std::vector<Canvas_point*>                           Canvas_point_vector;

  typedef std::vector<Canvas_point>                            Voronoi_vertices_container;
  typedef std::vector<Voronoi_vertices_container>            Voronoi_vertices_vector;

  // Refinement
  typedef boost::tuple<Simplex, std::size_t, FT>             PQ_entry;
  typedef PQ_entry_comparator<PQ_entry, Self>                PQ_comparer;
  typedef std::set<PQ_entry, PQ_comparer>                    PQ;

  // for fast seed locate (cheating by pretending we're in 3D)
  struct Base_mesh_primitive
  {
    typedef typename Star_set::Face_handle             Id;

    typedef K::Point_3                                Point;
    typedef K::Triangle_3                             Datum;

    Id m_it;
    Datum m_datum; // cache the datum
    Point m_p; // cache the reference point

    Base_mesh_primitive() { }
    Base_mesh_primitive(Id fh_, const Datum& datum, const Point& point)
      :
        m_it(fh_),
        m_datum(datum),
        m_p(point)
    { }

     Id id() const { return m_it; }
     Datum datum() const { return m_datum; }
     Point reference_point() const { return m_p; }
  };

  typedef Base_mesh_primitive                             Primitive;
  typedef Primitive::Id                                   Primitive_Id;
  typedef CGAL::AABB_traits<K, Primitive>                 AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits>           Tree;

  // The starset
  Star_set& ss;

  std::vector<Canvas_point> points;
  std::vector<std::size_t> trial_points;

  // Duality
  mutable boost::unordered_set<Edge> dual_edges;
  mutable boost::unordered_set<Tri> dual_triangles;

  // Constraint
  std::vector<Voronoi_constraint> Vor_cons;

  // Refinement
  PQ size_queue;
  PQ distortion_queue;
  PQ intersection_queue;
  PQ quality_queue;

  // Llyod
  Voronoi_vertices_vector Voronoi_vertices;

  // AABB tree
  Tree tree;

  void print_states() const
  {
    std::cout << "known: " << known_count;
    std::cout << " trial: " << trial_count;
    std::cout << " far: " << far_count;
    std::cout << " (trial: " << trial_points.size() << ",";
    std::cout << " seeds: " << seeds.size() << ")" << std::endl;
  }

  std::size_t add_corner_to_corners_map(const Point_2& p,
                                        std::map<Point_2, std::size_t>& corners)
  {
    // check if the corner exists
    typedef typename std::map<Point_2, std::size_t>::iterator Map_it;
    std::pair<Map_it, bool> is_insert_successful =
                                          corners.insert(std::make_pair(p, 0));

    if(is_insert_successful.second) // didn't exist in the corners map
      is_insert_successful.first->second = insert_new_seed(p.x(), p.y(), ON_CORNER);

    return is_insert_successful.first->second;
  }

  void discretize_constraint(const Voronoi_constraint& c,
                             std::map<Point_2, std::size_t>& corners)
  {
    FT sq_dist_from_p1 = 0.;
    std::size_t seeds_n_before_discretization = seeds.size();

    Point_2 current_point = c.p1;
    add_corner_to_corners_map(c.p1, corners);

    while(true)
    {
      // not the most efficient to compute so many metrics...
      const Metric& m = mf->compute_metric(current_point);

      Vector2d transf_dir = m.get_transformation() * c.dir;
      FT nm1 = 1. / transf_dir.norm();

      FT new_x = current_point.x() + c.dir(0) * nm1;
      FT new_y = current_point.y() + c.dir(1) * nm1;
      current_point = Point_2(new_x, new_y);
      sq_dist_from_p1 = CGAL::squared_distance(c.p1, current_point);

      if(sq_dist_from_p1 >= c.sq_length)
        break;

      insert_new_seed(new_x, new_y, ON_CONSTRAINT);
    }

    add_corner_to_corners_map(c.p2, corners);

    std::cout << "Constraint: " << c.p1 << " || " << c.p2 << std::endl;
    std::cout << "added: " << seeds.size() - seeds_n_before_discretization << std::endl;
  }

  void discretize_constraints()
  {
    // hardcoded constraints for now : a square that is the bbox of the seeds
    Point_2 p1(center.x() - sa, center.y() - sb);
    Point_2 p2(center.x() + sa, center.y() - sb);
    Point_2 p3(center.x() + sa, center.y() + sb);
    Point_2 p4(center.x() - sa, center.y() + sb);
    Voronoi_constraint c1(p1, p2);
    Voronoi_constraint c2(p2, p3);
    Voronoi_constraint c3(p3, p4);
    Voronoi_constraint c4(p4, p1);
    Vor_cons.push_back(c1);
    Vor_cons.push_back(c2);
    Vor_cons.push_back(c3);
    Vor_cons.push_back(c4);

    std::map<Point_2, std::size_t> corners; // point to seed id

    discretize_constraint(c1, corners);
    discretize_constraint(c2, corners);
    discretize_constraint(c3, corners);
    discretize_constraint(c4, corners);
  }

  Seed_status determine_seed_status(const FT x, const FT y)
  {
    std::size_t vcs = Vor_cons.size();
    std::size_t contained_by_n = 0.;
    for(std::size_t i=0; i<vcs; ++i)
    {
      const Voronoi_constraint& vc = Vor_cons[i];
      if(vc.has_point(Point_2(x, y)))
        contained_by_n++;
      // could put a break if 2+ I guess...
    }

    if(contained_by_n >= 2)
      return ON_CORNER;
    else if(contained_by_n >= 1)
      return ON_CONSTRAINT;
    else
      return IN_DOMAIN;
  }

  std::size_t insert_new_seed(const FT x, const FT y,
                              const Seed_status s = IN_DOMAIN)
  {
  #ifdef FILTER_SEEDS_OUTSIDE_GRID
      if(x < center.x() - sa || x > center.x() + sa ||
         y < center.y() - sb || y > center.y() + sb)
      {
  #if (verbose > 1)
        std::cout << "filtered : " << x << " " << y << std::endl;
  #endif
        return seeds.size();
      }
  #endif
  #if (verbose > 0)
    std::cout << "added new seed: " << x << " " << y;
    std::cout << " (" << seeds.size() << ")";
    std::cout << " with status: " << s << std::endl;
  #endif
    seeds.push_back(Point_2(x, y));
    seeds_m.push_back(mf->compute_metric(seeds.back()));
    seeds_s.push_back(s);
    return seeds.size();
  }

  int build_seeds()
  {
    std::cout << "build seeds" << std::endl;

    // don't want to count constraint seeds as 'real' seeds
    std::size_t constraint_seeds_n = seeds.size();

    std::ifstream in((str_seeds + ".mesh").c_str());

    if(!in)
    {
      std::cout << "couldn't open seed file" << std::endl;
      exit(0);
    }

    std::string word;
    std::size_t useless, nv, dim;
    FT r_x, r_y;

    in >> word >> useless; // MeshVersionFormatted i
    in >> word >> dim; // Dimension d
    in >> word >> nv;
    std::cout << "seeds nv: " << nv << std::endl;
    CGAL_assertion(dim == 2);

    std::size_t min_nv = (std::min)(nv, vertices_nv);

    seeds.reserve(min_nv);
    seeds_m.reserve(min_nv);

    for(std::size_t i=0; i<nv; ++i)
    {
      in >> r_x >> r_y >> useless;

      Seed_status s = determine_seed_status(r_x, r_y);
      insert_new_seed(r_x, r_y, s);

      if(seeds.size() - constraint_seeds_n >= vertices_nv)
        break;
    }
  #if (verbose > 0)
    std::cout << "seeds: " << seeds.size() << std::endl;
  #endif
    return seeds.size();
  }

  void initialize_seeds()
  {
    // initialize constraints first
    if(are_constraints_discretized)
      discretize_constraints();

    vertices_nv = build_seeds();
    CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );
  }

  void initialize_from_starset()
  {
    std::cout << "ini bm from ss " << std::endl;

    base_mesh_bbox = CGAL::Bbox_2();

    // build the points from the starset
    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];
      bool is_on_border = has_a_corner_vertex(star);
      Canvas_point cp(this, i, is_on_border);
      points.push_back(cp);

      base_mesh_bbox += points[i].point.bbox();

#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
      // build reciprocal adjacency
      Vertex_handle_handle vit =  star->finite_adjacent_vertices_begin();
      Vertex_handle_handle vend =  star->finite_adjacent_vertices_begin();
      for(; vit!=vend; ++vit)
      {
        Vertex_handle vh = *vit;
        Star_handle star_i = ss[vh->info()];
        if(!star_i->has_vertex(i))
          points[vh->info()].reciprocal_neighbors.push_back(i);
      }
#endif
    }

    std::cout << "bbox: " << base_mesh_bbox.xmin() << " " << base_mesh_bbox.xmax();
    std::cout << " " << base_mesh_bbox.ymin() << " " << base_mesh_bbox.ymax() << std::endl;

    build_aabb_tree();
  }

  void initialize_canvas_point(std::size_t id,
                               const FT dist,
                               const std::size_t seed_id)
  {
    Canvas_point& cp = points[id];
    if(cp.closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cout << "WARNING: a new seed (" << seed_id << ") is overwriting the closest seed id at ";
      std::cout << cp.index << ". Previous closest seed is: " << cp.closest_seed_id;
      std::cout << " (at: " << seeds[cp.closest_seed_id] << ")" << std::endl;
   }

    if(cp.state == TRIAL)
    {
      // it's already a trial point, we can't accept two seeds for one canvas point
      std::cout << "Warning: the canvas is not dense enough for the input..."
                << " Ignoring seed: " << seed_id << " : " << seeds[seed_id] << std::endl;
      return;
    }

    cp.initialize_from_point(dist, seed_id);
    trial_points.push_back(cp.index);
    std::push_heap(trial_points.begin(), trial_points.end(),
                   Canvas_point_comparer<Self>(this));

#ifdef OVERWRITE_SEED_WITH_CANVAS_POINT
    seeds[seed_id] = cp.point;
    cp.distance_to_closest_seed = 0;
#endif
  }

  void find_and_initialize_closest_vertex(const Point_2& s,
                                          const std::size_t seed_id)
  {
    // this is pretty expensive, but only used when the aabb tree failed to find
    // a star

    FT min_d = FT_inf;
    std::size_t min_id = -1;
    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];
      K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();
      FT d = sqd(star->center_point(), s);
      if(d < min_d)
      {
        min_id = i;
        min_d = d;
      }
    }
    initialize_canvas_point(min_id, std::sqrt(min_d), seed_id);
  }

  void locate_and_initialize(const Point_2& s,
                             const std::size_t seed_id)
  {
    // find the triangle that contains the seed and initialize the distance at
    // these three points accordingly.
    bool found = false;

    Segment_3 query(Point_3(s.x(), s.y(), -1.), Point_3(s.x(), s.y(), 1.));
    std::list<Primitive_Id> intersections;
    tree.all_intersected_primitives(query, std::back_inserter(intersections));

    if(intersections.empty())
    {
      std::cout << "the locate-ray didn't encounter any star..." << std::endl;
      std::cout << "bbox: " << base_mesh_bbox << std::endl;
      std::cout << "seed: " << seeds[seed_id] << std::endl;

      // switch to something expensive... take the closest point...
      find_and_initialize_closest_vertex(s, seed_id);
      found = true;
      return;
    }

//    if(intersections.size() > 1)
//      std::cout << "more than one intersection in the locate" << std::endl;

    std::size_t closest_canvas_point_id = -1;
    FT distance_to_closest_canvas_point = 1e30;
    std::list<Primitive_Id>::const_iterator it = intersections.begin(),
                                            end = intersections.end();
    for(; it!=end; ++it)
    {
      const typename Star_set::Face_handle fh = *it;
      const Triangle_2 tr_2(points[fh->vertex(0)->info()].point,
                            points[fh->vertex(1)->info()].point,
                            points[fh->vertex(2)->info()].point);

      if(tr_2.bounded_side(s) != CGAL::ON_UNBOUNDED_SIDE)
      {
        found = true;
#if (verbose > 5)
        std::size_t i0 = fh->vertex(0)->info();
        std::size_t i1 = fh->vertex(1)->info();
        std::size_t i2 = fh->vertex(2)->info();

        std::cout << "locating seed " << seed_id
                  << " point: " << s.x() << " " << s.y() << std::endl;
        std::cout << "found triangle: " << std::endl;
        std::cout << i0 << " [" << points[i0].point << "] " << std::endl;
        std::cout << i1 << " [" << points[i1].point << "] " << std::endl;
        std::cout << i2 << " [" << points[i2].point << "] " << std::endl;
#endif

        // we're inside! compute the distance to the vertices of the triangle
        for(int j=0; j<3; ++j)
        {
          Canvas_point& cp = points[fh->vertex(j)->info()];
          const Metric& seed_m = mf->compute_metric(s);
          const Metric& v_m = cp.metric;
          const Eigen::Matrix2d& f = get_interpolated_transformation(seed_m, v_m);

          Vector2d v;
          v(0) = s.x() - cp.point.x();
          v(1) = s.y() - cp.point.y();
          v = f * v;
          FT d = v.norm(); // d = std::sqrt(v^t * M * v) = (f*v).norm()

          if(d < distance_to_closest_canvas_point)
          {
            closest_canvas_point_id = cp.index;
            distance_to_closest_canvas_point = d;
          }
        }
      }
    }

    Canvas_point& cp = points[closest_canvas_point_id];
    initialize_canvas_point(cp.index, distance_to_closest_canvas_point, seed_id);
    CGAL_assertion(found);
  }

  void build_aabb_tree()
  {
    tree.clear();

    // create an aabb tree of triangles to fasten it
    std::size_t filtered_stars_n = 0.;

    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];

      Face_handle_handle fhit = star->finite_incident_faces_begin();
      Face_handle_handle fhend = star->finite_incident_faces_end();
      for(; fhit!=fhend; ++fhit)
      {
        typename Star_set::Face_handle fh = *fhit;

        if(is_face_outside_domain(star, fh))
        {
          filtered_stars_n++;
          continue;
        }

        const Point_2& p0 = points[fh->vertex(0)->info()].point;
        const Point_2& p1 = points[fh->vertex(1)->info()].point;
        const Point_2& p2 = points[fh->vertex(2)->info()].point;

        Point_3 p(p0.x(), p0.y(), 0.);
        Triangle_3 tr(p,
                      Point_3(p1.x(), p1.y(), 0.),
                      Point_3(p2.x(), p2.y(), 0.));
        Primitive pri(fh, tr, p);
        tree.insert(pri);
      }
    }

    std::cout << "built aabb tree w/ " << tree.size() << " entries" << std::endl;
    std::cout << filtered_stars_n << " filtered stars" << std::endl;

    // if the assertion below crashes, you probably forgot to change the domain
    // and thus the faces get filtered...
    CGAL_postcondition(tree.size() > ss.size());
  }

  void locate_and_initialize_seeds()
  {
    for(std::size_t i=0, ss=seeds.size(); i<ss; ++i)
    {
      const Point_2& p = seeds[i];
      locate_and_initialize(p, i);
    }
    Voronoi_vertices.resize(seeds.size());
  }

  void dual_shenanigans(Canvas_point& cp)
  {
    CGAL_assertion(cp.state == KNOWN);
    const Star_handle star = ss[cp.index];

    Face_handle_handle fhit = star->finite_incident_faces_begin();
    Face_handle_handle fhend = star->finite_incident_faces_end();
    for(; fhit!=fhend; ++fhit)
    {
      typename Star_set::Face_handle fh = *fhit;
      if(is_face_outside_domain(star, fh))
        continue;

      Simplex dual_simplex;
      for(std::size_t j=0; j<3; ++j)
      {
        const Canvas_point& cq = points[fh->vertex(j)->info()];
        if(cq.state != KNOWN)
          continue;

        std::size_t q_id = cq.closest_seed_id;
        CGAL_assertion(q_id != static_cast<std::size_t>(-1));
        dual_simplex.insert(q_id);
      }
      add_simplex_to_triangulation(cp, dual_simplex);
    }
  }

  bool get_next_trial_point(std::size_t& id)
  {
    while(!trial_points.empty())
    {
      id = trial_points.front();
      std::pop_heap(trial_points.begin(), trial_points.end(),
                    Canvas_point_comparer<Self>(this));
      trial_points.pop_back();

      CGAL_assertion(id < points.size());

      const Canvas_point& cp = points[id];
      if(cp.state != TRIAL)
      {
        std::cout << "WARNING : point with a state non-TRIAL in the PQ : ";
        std::cout << cp.index << "... Ignoring it!" << std::endl;
      }
      else
        return true;
    }
    return false;
  }

  void fix_distance_and_seed_id_at_corners()
  {
    // corners are points outside the domain lazily added by the starset so
    // it doesn't have to deal with the refinement of infinite cells
    // these points are the first four stars
    for(std::size_t i=0; i<4; ++i)
    {
      points[i].closest_seed_id = -1;
      points[i].distance_to_closest_seed = 0;
    }
  }

  void spread_distances(const bool use_dual_shenanigans)
  {
#if (verbose > 5)
    std::cout << "main loop" << std::endl;
#endif
    std::clock_t s_start = std::clock();

    std::size_t id = -1;
    bool is_t_empty = trial_points.empty();

    if(is_t_empty)
      std::cout << "WARNING: trial points shouldn't be empty !" << std::endl;

    while(!is_t_empty)
    {
      if(known_count%10000 == 0)
        print_states();

#if (verbose > 10)
      std::cout << "Trial queue size : " << trial_points.size() << std::endl;
#endif

#if (verbose > 30)
      std::cout << "trial heap: " << std::endl;
#if (verbose > 40)
      for (std::vector<Canvas_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif
#endif

      if(!get_next_trial_point(id))
        break;

      Canvas_point& cp = points[id];

#if (verbose > 10)
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "picked n° " << cp.index << " (" << cp.point << ") ";
      std::cout << "at distance : " << cp.distance_to_closest_seed << " from the seed "
                << cp.closest_seed_id << " ancestor : " << cp.ancestor << std::endl;
#endif

      cp.state = KNOWN;
      known_count++; // tmp --> use change_state()

      if(false && use_dual_shenanigans) // tmp
        dual_shenanigans(cp);

      PQ_state pqs = cp.update_neighbors_distances(trial_points);
      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(),
                       Canvas_point_comparer<Self>(this));
      is_t_empty = trial_points.empty();
    }

    fix_distance_and_seed_id_at_corners();
    CGAL_expensive_assertion_code(debug());

    std::cout << "End of spread_distances. time: ";
    std::cout << ( std::clock() - s_start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void spread_distances_from_one_seed(const std::size_t seed_id,
                                      std::map<std::size_t,
                                               boost::unordered_set<std::size_t> > neighbors,
                                      std::vector<std::deque<std::pair<std::size_t, FT> > >& geodesics,
                                      std::vector<Canvas_point>& canvas_point_memory,
                                      std::map<std::pair<std::size_t, std::size_t>,
                                               FT>& geodesic_lengths)
  {
    // Spread distances from the seed_id, until it has reached its neighboring seeds
    CGAL_precondition(trial_points.empty());

    boost::unordered_set<std::size_t>& seeds_to_reach = neighbors[seed_id];
    std::cout << "neighbors to reach : ";
    typename boost::unordered_set<std::size_t>::iterator it = seeds_to_reach.begin();
    typename boost::unordered_set<std::size_t>::iterator end = seeds_to_reach.end();
    for(; it!=end; ++it)
      std::cout << *it << " ";
    std::cout << std::endl;

    // (really) uglily for now, restart from the seed by searching the previous seed holder...
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Canvas_point& cp = points[i];

      if(static_cast<std::size_t>(cp.is_seed_holder) == seed_id)
      {
        CGAL_assertion(cp.closest_seed_id == seed_id);
        CGAL_assertion(cp.ancestor == static_cast<std::size_t>(-1));
        CGAL_assertion(cp.depth == 0);

        cp.state = TRIAL;
        trial_points.push_back(cp.index);
        std::push_heap(trial_points.begin(), trial_points.end(),
                       Canvas_point_comparer<Self>(this));
        break;
      }
    }

    std::size_t id = -1;
    bool is_t_empty = trial_points.empty();
    CGAL_precondition(!is_t_empty);

    FT max_distance_to_neighboring_seed = -1e30;
    known_count = 0;
    while(!is_t_empty)
    {
      if(known_count%10000 == 0)
        print_states();

#if (verbose > 10)
      std::cout << "Trial queue size : " << trial_points.size() << std::endl;
#endif

#if (verbose > 40)
      std::cout << "trial heap: " << std::endl;
      for (std::vector<Canvas_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif

      if(!get_next_trial_point(id))
        break;

      Canvas_point& cp = points[id];

#if (verbose > 10)
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "picked n° " << cp.index << " (" << cp.point.x() << ", " << cp.point.y() << ") ";
      std::cout << "at distance : " << cp.distance_to_closest_seed << " from the seed "
                << cp.closest_seed_id << " ancestor : " << cp.ancestor << std::endl;
#endif

      // Check if we have reached one of the required seeds
      if(cp.is_seed_holder >= 0 && // really ugly, fixme (see comment at member def)
         seeds_to_reach.find(cp.is_seed_holder) != seeds_to_reach.end() )
      {
        std::cout << "found the seed " << cp.is_seed_holder;
        std::cout << " at " << cp.index << " (" << cp.point << ") " << std::endl;

        // cp is the seed holder of another seed, but now it's distance field
        // is the distance to seed_id
        max_distance_to_neighboring_seed = cp.distance_to_closest_seed;

        // That's the center of a Voronoi cell (no ancestor)
        std::size_t reached_seed_id = cp.is_seed_holder;
        CGAL_assertion(cp.closest_seed_id == seed_id);

        // Geodesic are symmetrical, no need to compute it both ways.
        // This also dodges the first point of the trial queue (which will be the
        // seed_holder for the seed 'seed_id' and we don't care about it
        if(reached_seed_id < seed_id)
          continue;

        boost::unordered_set<std::size_t>::iterator sit =
                                            seeds_to_reach.find(reached_seed_id);
        if(sit != seeds_to_reach.end())
        {
          // the seed that we have reached is part of the neighbors we want to reach
#ifdef COMPUTE_GEODESIC
          extract_geodesic(cp, geodesics);
#endif
          seeds_to_reach.erase(sit);
        }
      }

      // only spread as little as possible
#ifdef COMPUTE_REAL_GEODESICS
      // alpha is a buffer so the gradient is well defined near all the seeds
      FT alpha = 1.2;
      if(seeds_to_reach.empty() &&
         cp.distance_to_closest_seed > alpha * max_distance_to_neighboring_seed)
        break;
#else
      if(seeds_to_reach.empty())
         break;
#endif

      cp.state = KNOWN;
      known_count++; // tmp --> use change_state()

      // Change the distance_to_closest_seed of adjacent points to allow
      // for seed_id's cell to spread
      Star_handle star = ss[id]; // the star corresponding to cp
      Vertex_handle_handle vit = star->finite_adjacent_vertices_begin();
      Vertex_handle_handle vend = star->finite_adjacent_vertices_end();
      for(; vit!=vend; ++vit)
      {
        Vertex_handle vh = *vit;

        if(is_a_corner_vertex(vh))
          continue;

        Canvas_point& cq = points[vh->info()];

        if(cq.state == FAR)
        {
          canvas_point_memory.push_back(cq);
          cq.distance_to_closest_seed = FT_inf;
        }
      }

      cp.update_neighbors_distances(trial_points);

      // always rebuild the priority queue
      std::make_heap(trial_points.begin(), trial_points.end(),
                     Canvas_point_comparer<Self>(this));
      is_t_empty = trial_points.empty();
    }

    std::cout << "reached all the seeds in " << known_count << " points" << std::endl;

    // make sure we've reached all the neighboring seeds
    CGAL_postcondition(seeds_to_reach.empty());
  }

  void add_to_neighboring_map(const Edge& e,
                              std::map<std::size_t,
                                       boost::unordered_set<std::size_t> >& m)
  {
    typedef std::map<std::size_t, boost::unordered_set<std::size_t> >  Map;

    std::size_t min = e[0], max = e[1];

    if(e[0] < 0 || e[0] >= seeds.size() || e[1] < 0 || e[1] >= seeds.size())
    {
      std::cout << "Warning: trying to add weird stuff in the neighboring map: "
                << e[0] << " " << e[1] << std::endl;
      return;
    }

    CGAL_precondition(min != max);
    if(min > max)
    {
      min = e[1];
      max = e[0];
    }
    CGAL_postcondition(min < max);

    std::pair<typename Map::iterator, bool> is_insert_successful;

    boost::unordered_set<std::size_t> s;
    s.insert(max);

    is_insert_successful = m.insert(std::make_pair(min, s));
    if(!is_insert_successful.second)
      (is_insert_successful.first)->second.insert(max);
  }

  void build_neighboring_info(std::map<std::size_t,
                                       boost::unordered_set<std::size_t> >& neighbors)
  {
    boost::unordered_set<Edge>::iterator it = dual_edges.begin();
    boost::unordered_set<Edge>::iterator end = dual_edges.end();

    for(; it!=end; ++it)
    {
      const Edge& e = *it;
      add_to_neighboring_map(e, neighbors);
    }
  }

  void rollback_points(const std::vector<Canvas_point>& cp_memory)
  {
    std::cout << "rollback" << std::endl;
    for(std::size_t i=0, cpms=cp_memory.size(); i<cpms; ++i)
    {
      const Canvas_point& cp_m = cp_memory[i];
      std::size_t id = cp_m.index;
      CGAL_assertion(points[id].point == cp_m.point);
      points[id] = cp_m;
    }

    trial_points.clear();
  }

  std::size_t ancestor_with_depth(const Canvas_point& cp)
  {
    std::size_t depth = cp.depth;

    if(depth == 0)
      return -1;

    CGAL_precondition(depth > 0);
    std::size_t anc = cp.ancestor;
    --depth;

    while(depth > 0)
    {
      CGAL_precondition(anc != static_cast<std::size_t>(-1));
      anc = points[anc].ancestor;
      --depth;
    }

    CGAL_postcondition(anc != static_cast<std::size_t>(-1));
    return anc;
  }

  void extract_geodesic(const Canvas_point& cp,
                        std::vector<std::deque<std::pair<std::size_t, FT> > >& geodesics)
  {
    std::deque<std::pair<std::size_t, FT> > geodesic_path;
    FT total_length = cp.distance_to_closest_seed;
    geodesic_path.push_front(std::make_pair(cp.index, 1.));
    CGAL_assertion(total_length != 0.);
    FT total_length_inv = 1./total_length;

#define EXTRACT_WITH_DEPTH
#ifdef EXTRACT_WITH_DEPTH
    std::size_t anc = ancestor_with_depth(cp);
#else
    std::size_t anc = cp.ancestor;
#endif
    CGAL_precondition(anc != static_cast<std::size_t>(-1));

    while(anc != static_cast<std::size_t>(-1))
    {
      FT t = points[anc].distance_to_closest_seed * total_length_inv;
      if(t < 0) t = 0.;

      if(t >= 1.)
      {
        t = 0.999999999;
        std::cout << "t: " << t << std::endl;
      }

      CGAL_postcondition(t >= 0. && t < 1.);

      geodesic_path.push_front(std::make_pair(anc, t));
#ifdef EXTRACT_WITH_DEPTH
      anc = ancestor_with_depth(points[anc]);
#else
      anc = points[anc].ancestor;
#endif
    }

    std::size_t s1 = cp.closest_seed_id;
    int s2 = cp.is_seed_holder;

    CGAL_precondition(s2 >= 0);
    CGAL_precondition(s1 < static_cast<std::size_t>(s2));

    std::cout << "built geodesic of " << s1 << " " << s2 << std::endl;
    geodesics.push_back(geodesic_path);
  }

  // real geodesics
  void draw_gradient(const std::vector<Eigen::Vector2d>& vertex_gradients,
                     const std::size_t seed_id)
  {
    std::ostringstream oss;
    oss << "geodesics_gradient_" << seed_id << ".mesh";
    std::ofstream out(oss.str().c_str());
    std::size_t ps = points.size();
    std::size_t vertex_counter = 1;

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';

    out << "Vertices" << '\n';
    out << 4 * ps << '\n';
    for(std::size_t i=0; i<ps; ++i)
    {
      Vector_2 v(vertex_gradients[i](0), vertex_gradients[i](1));
      Vector_2 v_orth(vertex_gradients[i](1), -vertex_gradients[i](0));

      out << points[i].point - 0.0005 * v_orth << " " << vertex_counter++ << '\n';
      out << points[i].point + 0.0005 * v_orth << " " << vertex_counter++ << '\n';

      out << points[i].point << " " << vertex_counter++ << '\n';
      out << points[i].point + 0.01 * v << " " << vertex_counter++ << '\n';
    }

    out << "Triangles" << '\n';
    out << 2 *ps << '\n';

    for(std::size_t i=0; i<ps; ++i)
    {
      out << 4*i+1 << " " << 4*i+2 << " " << 4*i+1 << " " << i << std::endl;
      out << 4*i+3 << " " << 4*(i+1) << " " << 4*i+3 << " " << i << std::endl;
    }

    out << "End" << std::endl;
  }

  void trace_geodesic(const std::vector<Eigen::Vector2d>& vertex_gradients,
                      std::map<std::pair<std::size_t, std::size_t>,
                               std::deque<Point_2> >& real_geodesics,
                      const std::size_t seed_id,
                      const std::size_t seed_to_reach,
                      const std::size_t starting_id)
  {
    std::cout << "tracing (gradient) geodesic: " << seed_id << " "
                                             << seed_to_reach << std::endl;

    // tracing the geodesic from starting_id to seed_id
    CGAL_precondition(starting_id < points.size());

    std::pair<std::size_t, std::size_t> geo_id = std::make_pair(seed_id,
                                                                seed_to_reach);
    FT tau = 0.0001; // gradient step

//    FT dist_at_seed_to_reach = points[starting_id].distance_to_seed;

    bool is_seed_reached = false;
    std::size_t iterations = 0.;
    Point_2 current_point = points[starting_id].point;

    std::deque<Point_2> geodesic;
    geodesic.push_front(current_point);

    Eigen::Vector2d grad_at_current_point = Eigen::Vector2d::Zero();

    while(!is_seed_reached)
    {
      // find all faces that contain the point
      Segment_3 query(Point_3(current_point.x(), current_point.y(), -1.),
                      Point_3(current_point.x(), current_point.y(), 1.));
      std::list<Primitive_Id> intersections;
      tree.all_intersected_primitives(query, std::back_inserter(intersections));

//      std::cout << "current point: " << current_point << std::endl;
      if(intersections.empty())
      {
        std::cout << "Warning: no intersection in trace geodesic... ?" << std::endl;
        return;
      }

      std::list<Primitive_Id>::const_iterator it = intersections.begin(),
                                              end = intersections.end();
      for(; it!=end; ++it)
      {
        const typename Star_set::Face_handle fh = *it;

        std::size_t p_id = fh->vertex(0)->info();
        std::size_t q_id = fh->vertex(1)->info();
        std::size_t r_id = fh->vertex(2)->info();
        Point_2 p = points[fh->vertex(0)->info()].point;
        Point_2 q = points[fh->vertex(1)->info()].point;
        Point_2 r = points[fh->vertex(2)->info()].point;

        const Triangle_2 tr_2(p, q, r);

        if(tr_2.bounded_side(current_point) != CGAL::ON_UNBOUNDED_SIDE)
        {
          // compute the bary weights in the face
          Eigen::Vector2d v1, v2, v3;
          v1(0) = q.x() - p.x();
          v1(1) = q.y() - p.y();
          v2(0) = r.x() - p.x();
          v2(1) = r.y() - p.y();

          FT lambda_p, lambda_q, lambda_r;
          v3(0) = current_point.x() - p.x();
          v3(1) = current_point.y() - p.y();

          FT d11 = v1.dot(v1);
          FT d12 = v1.dot(v2);
          FT d22 = v2.dot(v2);
          FT d31 = v3.dot(v1);
          FT d32 = v3.dot(v2);
          FT den = 1. / (d11 * d22 - d12 * d12);
          lambda_q = (d22 * d31 - d12 * d32) * den;
          lambda_r = (d11 * d32 - d12 * d31) * den;
          lambda_p = 1. - lambda_q - lambda_r;

//          std::cout << "p: " << p << std::endl;
//          std::cout << "q: " << q << std::endl;
//          std::cout << "r: " << r << std::endl;
//          std::cout << "lambdas: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;

          CGAL_postcondition(lambda_p >= 0);
          CGAL_postcondition(lambda_q >= 0);
          CGAL_postcondition(lambda_r >= 0);

          FT x = lambda_p * p.x() + lambda_q * q.x() + lambda_r * r.x();
          FT y = lambda_p * p.y() + lambda_q * q.y() + lambda_r * r.y();
          if(std::abs(x - current_point.x()) > 1e-8 ||
             std::abs(y - current_point.y()) > 1e-8)
          {
            std::cout << "Warning: problem in the barycentric coordinates" << '\n';
            std::cout << current_point << " vs " << std::endl << x << " " << y << std::endl;
          }

          if(std::abs(vertex_gradients[p_id].norm() - 1.) > 1e-8 ||
             std::abs(vertex_gradients[q_id].norm() - 1.) > 1e-8 ||
             std::abs(vertex_gradients[r_id].norm() - 1.) > 1e-8)
          {
            std::cout << "Warning: problem in the normals" << '\n';
            std::cout << vertex_gradients[p_id] << '\n';
            std::cout << vertex_gradients[q_id] << '\n';
            std::cout << vertex_gradients[r_id] << std::endl;
          }

          // get gradient & value at point p
//          std::cout << "adding grad: " << (lambda_p * vertex_gradients[p_id] + lambda_q * vertex_gradients[q_id] + lambda_r * vertex_gradients[r_id]).transpose() << std::endl;

          grad_at_current_point += lambda_p * vertex_gradients[p_id] +
                                   lambda_q * vertex_gradients[q_id] +
                                   lambda_r * vertex_gradients[r_id];

//    std::cout << "closest point is: " << current_point << std::endl;
//    std::cout << "grad: " << grad_at_current_point << std::endl;
        }
      }
//      grad_at_current_point /= (intersections.size());

      CGAL_postcondition(grad_at_current_point.norm() != 0);
      grad_at_current_point.normalize();
//      std::cout << "grad: " << grad_at_current_point.transpose() << " ";

      // new point
      Vector2d cur_p;
      cur_p(0) = current_point.x();
      cur_p(1) = current_point.y();

      Vector2d new_p = cur_p - tau * grad_at_current_point;
      current_point = Point_2(new_p(0), new_p(1));

      geodesic.push_front(current_point);

      FT d = CGAL::squared_distance(current_point, seeds[seed_id]);
//      std::cout << "d: " << d << std::endl;

      // following doesn't solve cycles
      if(d < 1e-5) // hardcode fixme
        break;

      iterations++;
      if(iterations > 1e5)
      {
        std::cout << "time to debug stuff -----------------------------" << std::endl;
        break;
      }
    }

    geodesic.push_front(seeds[seed_id]); // to have a nice connection

    std::cout << "traced (gradient) geodesic: " << geo_id.first << " "
                                                << geo_id.second << std::endl;
    real_geodesics[geo_id] = geodesic;
  }

  void compute_gradient_at_face(Face_handle fh,
                                std::vector<Eigen::Vector2d>& vertex_gradients)
  {
    // compute the grad on facet f = (c,sec) and update each vertex with
    // 1 / d(p_i, g) grad(f)

    // could be made more efficient without a copy of ids and points... todo

    // get vertices
    boost::array<Point_2, 3> face_points;
    boost::array<std::size_t, 3> face_ids;
    for(std::size_t i=0; i<3; ++i)
    {
      face_ids[i] = fh->vertex(i)->info();
      face_points[i] = points[face_ids[i]].point;
    }

    // there's probably prettier than going to 3D, but whatever...
    boost::array<Eigen::Vector3d, 3> edges; // e01, e12, e20

    for(std::size_t i=0; i<3; ++i)
    {
      edges[i](0) = face_points[(i+1)%3].x() - face_points[i].x();
      edges[i](1) = face_points[(i+1)%3].y() - face_points[i].y();
      edges[i](2) = 0;
    }

    // compute normal
    Eigen::Vector3d n = edges[0].cross(-edges[2]);
    n.normalize();

    // compute area of the triangle
    typename K::Compute_area_2 o;
    FT A = std::abs(o(face_points[0], face_points[1], face_points[2]));

    // compute the gradient
    Eigen::Vector2d grad = Eigen::Vector2d::Zero();

    for(std::size_t i=0; i<3; ++i)
    {
      // edge i is e_{i,(i+1)%3}
      Eigen::Vector3d cp = n.cross(edges[i]);

      // value at the opposite vertex (i+2)%3
      FT val = points[face_ids[(i+2)%3]].distance_to_closest_seed;
      grad(0) += val * cp(0);
      grad(1) += val * cp(1);
    }
    grad *= 0.5 / A;

    // center of the triangle
    FT third = 1. / 3.;
    Point_2 bg = CGAL::barycenter(face_points[0], third,
                                  face_points[1], third,
                                  face_points[2], third);

    // add it to the vertex grad parts at each vertex
    typename K::Compute_squared_distance_2 sqd;
    for(std::size_t i=0; i<3; ++i)
    {
      FT d = sqd(bg, face_points[i]);
      CGAL_assertion(d != 0.);

      std::size_t id = face_ids[i];
      vertex_gradients[id] += grad / A;
    }
  }

  void trace_geodesic_with_gradient(
          std::size_t seed_id,
          std::map<std::pair<std::size_t, std::size_t>,
                   std::deque<Point_2> >& real_geodesics,
          std::map<std::size_t/*seed_id*/,
                   boost::unordered_set<std::size_t>/*neighbor ids*/> neighbors)
  {
    std::cout << "tracing (gradient) geodesics from " << seed_id << std::endl;

    std::size_t v_n = points.size();

    // compute gradient at all triangles/vertices of the canvas
    std::vector<Eigen::Vector2d> vertex_gradients(v_n, Eigen::Vector2d::Zero());

    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];
      Face_handle_handle fit = star->finite_incident_faces_begin();
      Face_handle_handle fend = star->finite_incident_faces_end();
      for(; fit!=fend; ++fit)
      {
        Face_handle fh = *fit;
        if(is_face_outside_domain(star, fh))
          continue;

        compute_gradient_at_face(fh, vertex_gradients);
      }
    }

    // normalize gradients
    for(std::size_t i=0; i<v_n; ++i)
      vertex_gradients[i].normalize();

    std::cout << "draw gradient at seed_id: " << seed_id << std::endl;
//    draw_gradient(vertex_gradients, seed_id);

    // trace geodesic
    const boost::unordered_set<std::size_t>& seeds_to_reach = neighbors[seed_id];
    std::cout << seeds_to_reach.size() << " seeds to reach (gradient)" << std::endl;

    typename boost::unordered_set<std::size_t>::iterator it = seeds_to_reach.begin();
    typename boost::unordered_set<std::size_t>::iterator end = seeds_to_reach.end();
    for(; it!=end; ++it)
    {
      std::size_t seed_to_reach = *it;

      if(seed_to_reach < seed_id)
        continue;

      // obviously really ugly and expensive
      std::size_t starting_id = -1;
      for(std::size_t i=0, ps=points.size(); i<ps; ++i)
      {
        if(points[i].is_seed_holder == seed_to_reach)
        {
          starting_id = i;
          break;
        }
      }

      CGAL_postcondition(starting_id != static_cast<std::size_t>(-1));

      std::cout << "trying to reach the seed at : " << seeds[seed_id];
      std::cout << " from seed at " << seeds[seed_to_reach] << std::endl;

      trace_geodesic(vertex_gradients, real_geodesics,
                     seed_id, seed_to_reach,
                     starting_id);
    }
  }

  void output_geodesic_lengths(const typename std::map<std::pair<std::size_t, std::size_t>,
                                                        FT>& geo_lengths)
  {
    std::vector<FT> values;

    typename std::map<std::pair<std::size_t, std::size_t>,
                      FT>::const_iterator it = geo_lengths.begin();
    typename std::map<std::pair<std::size_t, std::size_t>,
                      FT>::const_iterator end = geo_lengths.end();

    for(; it!=end; ++it)
      values.push_back(it->second);

    output_histogram(values, str_base_mesh + "_geodesic_lengths.txt");
  }

  void output_geodesics(const std::vector<std::deque<std::pair<std::size_t,
                                                     FT> > >& geodesics,
                        const std::string str_base)
  {
    if(geodesics.empty())
    {
      std::cout << "no geodesic to output" << std::endl;
      return;
    }

    std::vector<bool> is_used(points.size(), false);
#define GEO_FILTER_UNUSED_POINTS
#ifdef GEO_FILTER_UNUSED_POINTS
    for(std::size_t i=0; i<geodesics.size(); ++i)
    {
      const std::deque<std::pair<std::size_t, FT> >& geo_path = geodesics[i];
      std::size_t geo_path_size_m1 = geo_path.size()-1;
      for(std::size_t j=0; j<geo_path_size_m1; ++j)
      {
        const std::pair<std::size_t, FT>& p = geo_path[j];
        const std::pair<std::size_t, FT>& pp1 = geo_path[j+1];

        is_used[p.first] = true;
        is_used[pp1.first] = true;
      }
    }

    std::vector<std::size_t> renumbering(points.size(), -1);
    std::size_t n=0;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      if(is_used[i])
        renumbering[i] = n++;
    }
#endif

    std::ofstream out((str_base + "_geodesics.mesh").c_str());

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';

    out << "Vertices" << '\n';
    out << points.size() << '\n';
    for(std::size_t i=0; i<points.size(); ++i)
      if(is_used[i])
        out << points[i].point << " " << i+1 << '\n';

    std::size_t edges_n = 0.;
    for(std::size_t i=0; i<geodesics.size(); ++i)
      edges_n += geodesics[i].size() - 1;

    out << "Edges" << '\n';
    out << edges_n << '\n';
    for(std::size_t i=0; i<geodesics.size(); ++i)
    {
      const std::deque<std::pair<std::size_t, FT> >& geo_path = geodesics[i];
      std::size_t geo_path_size_m1 = geo_path.size()-1;
      for(std::size_t j=0; j<geo_path_size_m1; ++j)
      {
        const std::pair<std::size_t, FT>& p = geo_path[j];
        const std::pair<std::size_t, FT>& pp1 = geo_path[j+1];

        // +1 due to medit
        out << renumbering[p.first] + 1 << " "
            << renumbering[pp1.first] + 1 << " " << i << std::endl;
      }
    }

    out << "End" << std::endl;
  }

  void output_real_geodesics_map(const std::map<std::pair<std::size_t, std::size_t>,
                                                std::deque<Point_2> >& real_geodesics,
                                 const std::string str_base)
  {
    if(real_geodesics.empty())
    {
      std::cout << "no (gradient) geodesic map to output" << std::endl;
      return;
    }

    typedef typename std::map<std::pair<std::size_t, std::size_t>,
                              std::deque<Point_2> >               Geo_map;
    typedef typename Geo_map::const_iterator                      Geo_map_it;

    std::ofstream out((str_base + "_real_geodesics_map.txt").c_str());

    out << real_geodesics.size() << std::endl;

    Geo_map_it gmit = real_geodesics.begin();
    Geo_map_it gmend = real_geodesics.end();
    for(; gmit!=gmend; ++gmit)
    {
      const std::pair<std::size_t, std::size_t>& ids = gmit->first;
      const std::deque<Point_2>& geo = gmit->second;

      out << ids.first << " " << ids.second << " " << geo.size() << std::endl;
      for(std::size_t i=0; i<geo.size(); ++i)
        out << geo[i] << '\n';

      out << '\n';
    }
    out << std::endl;
  }

  void output_real_geodesics(const std::map<std::pair<std::size_t, std::size_t>,
                                            std::deque<Point_2> >& real_geodesics,
                             const std::string str_base)
  {
    if(real_geodesics.empty())
    {
      std::cout << "no (gradient) geodesic to output" << std::endl;
      return;
    }

    std::ofstream out((str_base + "_real_geodesics.mesh").c_str());

    typedef typename std::map<std::pair<std::size_t, std::size_t>,
                              std::deque<Point_2> >               Geo_map;
    typedef typename Geo_map::const_iterator                      Geo_map_it;
    typedef typename std::map<Point_2, std::size_t>               Numbering_map;
    typedef typename Numbering_map::iterator                      Nm_it;

    Numbering_map vertex_numbering_map;
    std::size_t n = 0;

    std::ostringstream out_v, out_t;
    std::size_t nt = 0;

    Geo_map_it gmit = real_geodesics.begin();
    Geo_map_it gmend = real_geodesics.end();
    for(; gmit!=gmend; ++gmit)
    {
      const std::deque<Point_2>& geo = gmit->second;
      for(std::size_t i=0; i<geo.size(); ++i)
      {
        const Point_2& p = geo[i];
        std::pair<Nm_it, bool> is_insert_successful =
            vertex_numbering_map.insert(std::make_pair(p,0));

        if(is_insert_successful.second)
        {
          (is_insert_successful.first)->second = ++n; // ++ first for medit
           out_v << p << " " << n << '\n';
        }


        // draw the edge with a flat triangle
        if(i > 0)
        {
          std::size_t prev = vertex_numbering_map[geo[i-1]];
          out_t << prev << " " << n << " " << prev << " " << "1" << '\n';
          nt++;
        }
      }
    }

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';

    out << "Vertices" << '\n';
    out << vertex_numbering_map.size() << '\n';
    out << out_v.str().c_str();

    out << "Triangles" << '\n';
    out << nt << '\n';
    out << out_t.str().c_str();

    out << "End" << std::endl;
  }

  void compute_geodesics(const std::string str_base)
  {
    std::cout << "computing geodesics" << std::endl;

    // For each seed, we must spread to the nearest neighbors to grab the geodesic
    // between both seeds
    CGAL_precondition(!dual_edges.empty());

    std::map<std::size_t/*seed_id*/,
             boost::unordered_set<std::size_t>/*neighbor ids*/> neighbors;
    build_neighboring_info(neighbors);

    std::cout << "built neighboring info : " << neighbors.size() << std::endl;

    // fixme change geodesics to a map too...
    std::vector<std::deque<std::pair<std::size_t/*id*/,FT/*time*/> > > geodesics;
    std::map<std::pair<std::size_t, std::size_t>, std::deque<Point_2> > real_geodesics;

    std::map<std::pair<std::size_t, std::size_t>, FT> geodesic_lengths;

    for(std::size_t seed_id=0, ss=seeds.size(); seed_id<ss; ++seed_id)
    {
      std::cout << "computing geodesics to : " << seed_id << std::endl;

      // below is kind of ugly... but it's clearer than switching everything to
      // 'FAR' earlier and then having the memory overwrite 'KNOWN' with 'FAR'
      // every loop iteration
      for(std::size_t i=0, ptss=points.size(); i<ptss; ++i)
        points[i].state = FAR;

      // since we're spreading from a seed over the other cells, we must keep
      // the previous colors in memory
      std::vector<Canvas_point> canvas_point_memory;

      spread_distances_from_one_seed(seed_id, neighbors, geodesics,
                                     canvas_point_memory, geodesic_lengths);

#ifdef COMPUTE_REAL_GEODESICS
      trace_geodesic_with_gradient(seed_id, real_geodesics, neighbors);
#endif

      rollback_points(canvas_point_memory); // reset the points we have overwritten
    }

    output_geodesics(geodesics, str_base);
    output_real_geodesics(real_geodesics, str_base);
    output_real_geodesics_map(real_geodesics, str_base);
    output_geodesic_lengths(geodesic_lengths);
    compute_Riemannian_angles(real_geodesics);
  }

  FT compute_Riemannian_angle(
      std::size_t origin,
      std::size_t extr_id_1,
      std::size_t extr_id_2,
      std::map<std::pair<std::size_t, std::size_t>,
               std::deque<Point_2> >& real_geodesics)
  {
    typedef typename std::deque<Point_2>                              Geodesic;

    // computes the angle between the tangents of the geodesics
    // ori -- extr_1 and ori -- extr_2

    std::cout << "orign: " << origin << std::endl;
    std::cout << "t1: " << extr_id_1 << std::endl;
    std::cout << "t2: " << extr_id_2 << std::endl;

    bool is_geo_inverted_1 = origin > extr_id_1;
    bool is_geo_inverted_2 = origin > extr_id_2;

    const Geodesic& g1 = is_geo_inverted_1 ? real_geodesics[std::make_pair(extr_id_1, origin)]
                                           : real_geodesics[std::make_pair(origin, extr_id_1)];

    const Geodesic& g2 = is_geo_inverted_2 ? real_geodesics[std::make_pair(extr_id_2, origin)]
                                           : real_geodesics[std::make_pair(origin, extr_id_2)];

    std::size_t g1s = g1.size();
    std::size_t g2s = g2.size();

    CGAL_precondition(!g1.empty() && !g2.empty());

    // the extremities of the two tangent vectors
    Point_2 origin_1 = is_geo_inverted_1 ? g1.back() : g1.front();
    Point_2 origin_2 = is_geo_inverted_2 ? g2.back() : g2.front();

    // origins should be roughly the same
    if(CGAL::squared_distance(origin_1, origin_2) > 1e-5)
    {
      std::cout << "Warning: origins too far apart" << std::endl;
      std::cout << origin_1 << std::endl << origin_2 << std::endl;
      std::cout << "inverted: " << is_geo_inverted_1 << " " << is_geo_inverted_2 << std::endl;
      exit(0);
    }

    // rather than simply take and g[1] or g[size-1 -1], we might want to
    // take something a bit farther to average the direction of the geodesic
    // over the last few iterations of the gradient descent method
    std::size_t averaging_step = 1;

    CGAL_precondition(g1.size() > averaging_step+1 &&
                      g2.size() > averaging_step+1);

    Point_2 tangent_1 = is_geo_inverted_1 ? g1[g1s - averaging_step]
                                          : g1[averaging_step];
    Point_2 tangent_2 = is_geo_inverted_2 ? g2[g2s - averaging_step]
                                          : g2[averaging_step];

    if(tangent_1 == tangent_2)
    {
      std::cout << "Warning: degenerate case with equal tangent segments" << std::endl;
    }

    // compute the angle and make sure it's within 0 pi
    const Metric& m = mf->compute_metric(origin_1);
    Point_2 t_o1p = m.transform(origin_1);
    Point_2 t_o2p = m.transform(origin_2);
    Point_2 t_t1p = m.transform(tangent_1);
    Point_2 t_t2p = m.transform(tangent_2);

    Vector_2 v1 = t_t1p - t_o1p;
    Vector_2 v2 = t_t2p - t_o2p;

    FT angle_cos = v1 * v2 / CGAL::sqrt(v1*v1) / CGAL::sqrt(v2 * v2);
    FT angle = std::acos(angle_cos);

    // make sure it's within 0 pi
    if(angle < 0)
      angle += CGAL_PI;

    std::cout << "angle " << angle << std::endl;
    CGAL_assertion(angle >= 0 && angle <= CGAL_PI);

     return angle;
  }

  // histograms
  void output_histogram(const std::vector<int>& histogram,
                        FT min, FT max,
                        std::string filename)
  {
    std::cout << "output: " << filename << std::endl;
    std::ofstream out(filename.c_str());
    std::size_t histo_n = histogram.size();
    std::cout << "histo_n: " << histo_n << std::endl;
    for(std::size_t i=0; i<histo_n; ++i)
    {
      FT val = min + (max-min)*((FT) i)/((FT) histo_n);
      out << i << "," << val << "," << histogram[i] << '\n';
    }
    out << std::endl;
  }

  void output_histogram(std::vector<FT>& values,
                        std::string filename)
  {
    FT min_value = *(std::min_element(values.begin(), values.end()));
    FT max_value = *(std::max_element(values.begin(), values.end()));

    std::cout << "Outputing values: " << values.size() << " " << min_value << " " << max_value << std::endl;

    int histogram_size = 1000;
    std::vector<int> histogram(histogram_size, 0);
    FT limit_val = histogram_size - 1.;
    FT step_size = (max_value - min_value) / (FT) histogram_size;

    for(std::size_t i=0; i<values.size(); ++i)
      histogram[(std::min)(limit_val, std::floor((values[i]-min_value)/step_size))]++;

    output_histogram(histogram, min_value, max_value, filename);
  }

  void build_external_real_geodesics(std::map<std::pair<std::size_t, std::size_t>,
                                     std::deque<Point_2> >& real_geodesics)
  {
    CGAL_precondition(real_geodesics.empty());

    std::ifstream in("shock_real_geodesics_map.mesh");

    std::size_t geo_map_s, id1, id2, geo_size;

    in >> geo_map_s; // number of geodesics

    for(std::size_t i=0; i<geo_map_s; ++i)
    {
      in >> id1 >> id2 >> geo_size;

      CGAL_precondition(id1 < id2);
      std::deque<Point_2> geodesic;
      std::pair<std::size_t, std::size_t> ids = std::make_pair(id1, id2);

      for(std::size_t j=0; j<geo_size; ++j)
      {
        Point_2 p;
        in >> p;
        geodesic.push_back(p);
      }

      real_geodesics[ids] = geodesic;
    }

    CGAL_postcondition(real_geodesics.size() == geo_map_s);
  }

  void compute_Riemannian_angles(std::map<std::pair<std::size_t, std::size_t>,
                                          std::deque<Point_2> >& real_geodesics)
  {
    std::cout << "computing Riemannian angles" << std::endl;

    // Evaluate the angles of the Riemannian simplices

    // The angles are computed at the vertex in the metric of the vertex
    // and are defined as the angle between the tangent vectors of the
    // geodesics

    if(real_geodesics.empty())
    {
      build_external_real_geodesics(real_geodesics);
    }

    // loop all the Riemannian triangles
    std::vector<FT> values;
    for(boost::unordered_set<Tri>::iterator it = dual_triangles.begin();
                                            it != dual_triangles.end(); ++it)

    {
      const Tri& tr = *it;
      CGAL_assertion(tr.size() == 3);

      for(std::size_t i=0; i<tr.size(); ++i)
      {
        std::size_t orig = tr[i];

        std::size_t dest_1 = tr[(i+1)%3];
        std::size_t dest_2 = tr[(i+2)%3];

        if(orig < 0 || orig >= seeds.size() ||
           dest_1 < 0 || dest_1 >= seeds.size() ||
           dest_2 < 0 || dest_2 >= seeds.size() )
        {
          std::cout << "Warning: tried to compute the angles in a weird triangle: "
                    << orig << " " << dest_1 << " " << dest_2 << std::endl;
          continue;
        }

        std::cout << "in triangle: " << tr[0] << " " << tr[1] << " " << tr[2] << std::endl;

        FT angle = compute_Riemannian_angle(orig, dest_1, dest_2, real_geodesics);
        values.push_back(angle);
      }
    }
    output_histogram(values, "histogram_riemannian_face_angles.cvs");
  }

  bool debug()
  {
    std::cout << "debugging with " << seeds.size() << " seeds" << std::endl;
    std::clock_t d_start = std::clock();

    CGAL_assertion(trial_points.empty() );
    bool failed = false;

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
#if (verbose > 25)
      points[i].print_ancestors();
      points[i].print_children();
#endif
      points[i].state = KNOWN;
    }

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      const Canvas_point& cp = points[i];

      if(cp.ancestor != static_cast<std::size_t>(-1) &&
         points[cp.ancestor].closest_seed_id != cp.closest_seed_id)
      {
         failed = true;
         std::cout << "ancestor/children don't agree on the seed @ " << cp.index << std::endl;
      }

      if(cp.ancestor != static_cast<std::size_t>(-1) &&
         points[cp.ancestor].children.find(cp.index) == points[cp.ancestor].children.end())
      {
        failed = true;
        std::cout << "failure in ancestor/children relationship at " << cp.index << std::endl;
        std::cout << "ancestor is : " << cp.ancestor << " color: " << cp.closest_seed_id << std::endl;
      }

      if(cp.distance_to_closest_seed == FT_inf || cp.closest_seed_id >= seeds.size() ||
         (cp.ancestor != static_cast<std::size_t>(-1) &&
          points[cp.ancestor].closest_seed_id != cp.closest_seed_id))
      {
        failed = true;
        std::cout << "debug: " << i << " [" << cp.point << "] " << " at distance : ";
        std::cout << cp.distance_to_closest_seed << " from " << cp.closest_seed_id;
        std::cout << " ancestor: " << cp.ancestor;
        std::cout << " ancid: ";
        if(cp.ancestor==static_cast<std::size_t>(-1))
           std::cout << "-1" << std::endl;
        else
          std::cout << points[cp.ancestor].closest_seed_id << std::endl;
        std::cout << std::endl;
      }
    }
//    if(failed)
//      return false;

    std::cout << "End of debug. time: ";
    std::cout << ( std::clock() - d_start ) / (double) CLOCKS_PER_SEC << std::endl;
    return true;
  }

  void clear_dual()
  {
    dual_edges.clear();
    dual_triangles.clear();

    size_queue.clear();
    distortion_queue.clear();
    intersection_queue.clear();
    quality_queue.clear();

    Voronoi_vertices.clear();
    Voronoi_vertices.resize(seeds.size());
  }

  void reset()
  {
    std::cout << "canvas reset" << std::endl;

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      Canvas_point& cp = points[i];
      cp.reset();
    }

    known_count = 0;
    trial_count = 0;
    far_count = points.size();

    clear_dual();
  }

  void refresh_canvas_point_states()
  {
    CGAL_assertion(trial_points.empty());
    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
      points[i].state = FAR;

    known_count = 0;
    trial_count = 0;
    far_count = points.size();
  }

  bool refine_seeds(const Point_2& new_seed)
  {
    std::size_t prev_amount_of_seeds = seeds.size();
    vertices_nv = insert_new_seed(new_seed.x(), new_seed.y());

    if(prev_amount_of_seeds == vertices_nv)
    {
      std::cout << "Warning: Insertion failed in refine_seeds()" << std::endl;
      return false;
    }

#ifdef USE_FULL_REBUILDS
    reset();
    locate_and_initialize_seeds();
#else
    // we can't spread from the new seed if all states are 'KNOWN'
    refresh_canvas_point_states();
    locate_and_initialize(new_seed, seeds.size()-1);
#endif

    clear_dual();
    spread_distances(false/*use_dual_shenanigans*/);
    compute_dual();

    return true;
  }

  bool refine_seeds_with_self_computed_ref_point()
  {
    std::cout << "state of the queues: "
              << " size: " << size_queue.size()
              << " distortion: " << distortion_queue.size()
              << " intersection: " << intersection_queue.size()
              << " quality: "  << quality_queue.size() << std::endl;

    PQ_entry best_entry;

//    std::cout << "queues: " << std::endl;
//    std::cout << "size queue : " << size_queue.size() << std::endl;
//    typename PQ::const_iterator pqcit = size_queue.begin();
//    for(; pqcit!=size_queue.end(); ++pqcit)
//      std::cout << *pqcit << std::endl;

    if(!size_queue.empty())
       best_entry = *(size_queue.begin());
    else if(!intersection_queue.empty())
      best_entry = *(intersection_queue.begin());
    else if(!distortion_queue.empty())
      best_entry = *(distortion_queue.begin());
    else if(!quality_queue.empty())
      best_entry = *(quality_queue.begin());
    else
    {
      std::cout << "Empty refinement queues" << std::endl;
      return false;
    }

    Point_2 ref_point = points[best_entry.get<1>()].point;

    std::cout << "REF: Picked : " << best_entry.get<1>() << "[ ";
    std::cout << ref_point << " ]";
    std::cout << " second: " << best_entry.get<2>() << std::endl;

    return refine_seeds(ref_point);
  }

  void output_canvas(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';
    out << "Vertices" << '\n';
    out << points.size() << '\n';
    out_bb << "2 1 " << points.size() << " 2" << '\n';

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      const Canvas_point& cp = points[i];
      out << cp.point << " " << i+1 << '\n';

      if(cp.distance_to_closest_seed == FT_inf)
        out_bb << -1. << '\n';
      else
        out_bb << cp.distance_to_closest_seed << '\n';
//      out_bb << cp.closest_seed_id << '\n';
//      out_bb << (cp.is_Voronoi_vertex?"1":"0") << '\n';
    }

    // count the faces in the starset
    std::size_t n_of_faces = 0;
    std::ostringstream out_faces;
    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];
      Face_handle_handle fit = star->finite_incident_faces_begin();
      Face_handle_handle fend = star->finite_incident_faces_end();
      for(; fit!=fend; ++fit)
      {
        Face_handle fh = *fit;
        if(is_face_outside_domain(star, fh))
          continue;

        ++n_of_faces;
        boost::unordered_set<std::size_t> materials;
        for(std::size_t j=0; j<3; ++j)
        {
          materials.insert(points[fh->vertex(j)->info()].closest_seed_id);
          out_faces << fh->vertex(j)->info() + 1 << " "; // +1 for medit
        }
        std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
        out_faces << mat << '\n';
      }
    }

    out << "Triangles" << '\n';
    out << n_of_faces << '\n';
    out << out_faces.str().c_str() << '\n';

    out_bb << "End" << '\n';
    out << "End" << '\n';
  }

  void add_triangles_edges_to_dual_edges(const Tri& tr) const
  {
    // no need to sort as long as we keep the indices in order
    Edge e;
    e[0] = tr[0]; e[1] = tr[1]; dual_edges.insert(e);
    e[0] = tr[0]; e[1] = tr[2]; dual_edges.insert(e);
    e[0] = tr[1]; e[1] = tr[2]; dual_edges.insert(e);
  }

  struct PQ_finder_comparer
  {
    PQ_entry entry;

    bool operator()(const PQ_entry& right)
    {
      const Simplex& s1 = entry.get<0>();
      const Simplex& s2 = right.get<0>();

      if(s1.size() != s2.size())
        return false;

      Simplex::const_iterator it1 = s1.begin();
      Simplex::const_iterator it2 = s2.begin();

      for(; it1!=s1.end(); ++it1, ++it2)
        if(*it1 != *it2)
          return false;
      return true;
    }

    PQ_finder_comparer(const PQ_entry& pq_entry) : entry(pq_entry) { }
  };

  void insert_in_PQ(const PQ_entry& entry, PQ& queue, int i)
  {
    // find if the element is already in the queue
    PQ::iterator it = std::find_if(queue.begin(), queue.end(),
                                   PQ_finder_comparer(entry));

    if(it != queue.end())
    {
      // that dual simplex is already on track to be refined. We want to increase
      // its priority only if the new refinement point is farther away

      const Canvas_point& git = points[it->get<1>()];
      const Canvas_point& ge = points[entry.get<1>()];

      if(git.distance_to_closest_seed < ge.distance_to_closest_seed)
      {
//        std::cout << "erase: " << *it;
        queue.erase(it);
      }
      else
        return;
    }

#if (verbose > 15)
    std::cout << "insert " << entry;
    std::cout << "queue id: " << i << std::endl;
#endif
    queue.insert(entry);
  }

  void test_simplex(const Canvas_point& cp,
                    const Simplex& dual_simplex)
  {
    // test the simplex 'dual_simplex' agaisnt the criteria to decide
    // whether it needs to be refined or not

    // Not sure if we should allow dual of edges to be considered for the size
    // and distortion criteria

    // Maybe only allow the dual of an edge to be considered if the dual is
    // on the border of the domain ?

/*
    std::cout << "ref shen with : ";
    std::set<std::size_t>::const_iterator it = dual_simplex.begin();
    for(; it != dual_simplex.end(); ++it)
      std::cout << *it << " ";
    std::cout << std::endl;
*/

    FT max_size = 1.; // <= 0 means unused
    FT max_distortion = 1.; // <= 1 means unused
    bool intersection_ref = false;
    FT min_qual = 0.; // <= 0 means unsused

    // size
    if(max_size > 0.)
    {
      FT cpd = cp.distance_to_closest_seed;
      if(cpd > max_size)
      {
        PQ_entry pqe = boost::make_tuple(dual_simplex, cp.index, cpd);
        insert_in_PQ(pqe, size_queue, 0);
        size_queue.insert(pqe);
        return;
      }
    }

    // distortion
    if(max_distortion > 1.)
    {
      FT gamma = cp.distortion_to_seed();
      if(gamma > max_distortion)
      {
        PQ_entry pqe = boost::make_tuple(dual_simplex, cp.index, gamma);
        insert_in_PQ(pqe, distortion_queue, 1);
        return;
      }
    }

    if(dual_simplex.size() == 3)
    {
      Tri tr;
      std::set<std::size_t>::iterator it = dual_simplex.begin();
      for(int i=0; i<3; ++i)
        tr[i] = *it++;

      // intersection
      if(intersection_ref)
      {
        if(is_triangle_intersected(tr, dual_edges))
        {
          Triangle_2 tri(seeds[tr[0]], seeds[tr[1]], seeds[tr[2]]);
          FT area = tri.area();
          PQ_entry pqe = boost::make_tuple(dual_simplex, cp.index, area);
          insert_in_PQ(pqe, intersection_queue, 2);
          return;
        }
      }

      // quality fixme quality should be radius edge ratio
      if(min_qual > 0.)
      {
        FT qual = compute_quality(tr);
        if(qual < min_qual)
        {
          CGAL_assertion(qual > 1e-17);
          // 1./qual to refine the worst element first
          PQ_entry pqe = boost::make_tuple(dual_simplex, cp.index, 1./qual);
          insert_in_PQ(pqe, quality_queue, 3);
          return;
        }
      }
    }
  }

  void test_simplex(const Tri& canvas_tr,
                    const Simplex& dual_simplex)
  {
    Tri::const_iterator it = canvas_tr.begin();
    Tri::const_iterator end = canvas_tr.end();
    for(; it!=end; ++it)
     test_simplex(points[*it], dual_simplex);
  }

  // Optimization --------------------------------------------------------------
  FT triangle_area_in_metric(const Canvas_point& cp,
                             const Canvas_point& cq,
                             const Canvas_point& cr) const
  {
    // compute the interpolated's metric transformation
    FT third = 1./3.;

    // very unefficient fixme
    Eigen::Matrix2d f = third*(cp.metric.get_transformation() +
                               cq.metric.get_transformation() +
                               cr.metric.get_transformation());

    // transform the points
    const Point_2& p = cp.point;
    const Point_2& q = cq.point;
    const Point_2& r = cr.point;

    Eigen::Vector2d v(p.x(), p.y());
    v = f*v;
    const Point_2 tp(v(0), v(1));

    v(0) = q.x(); v(1) = q.y();
    v = f*v;
    const Point_2 tq(v(0), v(1));

    v(0) = r.x(); v(1) = r.y();
    v = f*v;
    const Point_2 tr(v(0), v(1));

    K::Compute_area_2 o;
    return std::abs(o(tp, tq, tr));
  }

  Canvas_point Vor_vertex_on_edge(const std::size_t n_q,
                                  const std::size_t n_r,
                                  int call_count = 0)
  {
    // do the complicated split, see formula on notes todo
    // taking the easy way for now

    Point_2 centroid = CGAL::barycenter(points[n_q].point, 0.5,
                                        points[n_r].point, 0.5);

    std::size_t id = points.size();
    Canvas_point c(this, centroid, id,
                   mf->compute_metric(centroid),
                   true/*on border*/);

    // need to push it to points (temporarily) for compute_closest_seed() to work
    points.push_back(c); // WARNING: MIGHT INVALIDATE ALL REF AND POINTERS

    Canvas_point& cp = points.back();
    cp.compute_closest_seed(n_q);
    cp.compute_closest_seed(n_r);
    cp.state = KNOWN;

    if(call_count > max_depth)
      return cp;

    if(cp.closest_seed_id == points[n_q].closest_seed_id)
      return Vor_vertex_on_edge(cp.index, n_r, ++call_count);
    else
    {
      CGAL_assertion(cp.closest_seed_id == points[n_r].closest_seed_id);
      return Vor_vertex_on_edge(cp.index, n_q, ++call_count);
    }
  }

  Canvas_point Vor_vertex_in_triangle(const std::size_t n_q,
                                      const std::size_t n_r,
                                      const std::size_t n_s,
                                      int call_count = 0)
  {
    CGAL_precondition(n_q < points.size() &&
                      n_r < points.size() &&
                      n_s < points.size());

    // centroid is probably not the most optimal...
    Point_2 centroid = CGAL::centroid(points[n_q].point,
                                      points[n_r].point,
                                      points[n_s].point);

    std::size_t id = points.size();
    Canvas_point c(this, centroid, id,
                   mf->compute_metric(centroid),
                   false/*border info*/);
    // it's actually possible for the centroid to be on the domain border,
    // but it doesn't matter since it's not a real grid point
    // and we'll never use this border information

    // need to push it back to points anyway for compute_closest_seed to work
    points.push_back(c);

    Canvas_point& cp = points.back();
    cp.compute_closest_seed(n_q);
    cp.compute_closest_seed(n_r);
    cp.compute_closest_seed(n_s);
    cp.state = KNOWN;

    // another potential stop is if the (max) distance between the centroid
    // and a point of the triangle is below a given distance
    if(call_count > max_depth)
      return cp;

    if(cp.closest_seed_id == points[n_q].closest_seed_id)
      return Vor_vertex_in_triangle(id, n_r, n_s, ++call_count);
    else if(cp.closest_seed_id == points[n_r].closest_seed_id)
      return Vor_vertex_in_triangle(id, n_q, n_s, ++call_count);
    else
    {
      CGAL_assertion(cp.closest_seed_id == points[n_s].closest_seed_id);
      return Vor_vertex_in_triangle(id, n_q, n_r, ++call_count);
    }
  }

  Canvas_point compute_precise_Voronoi_vertex(const Tri& tri)
  {
    std::size_t real_points_n = points.size();
    Canvas_point centroid = Vor_vertex_in_triangle(tri[0], tri[1], tri[2]);

    // clean off the virtual points created by Vor_vertex_in_triangle
    shave_off_virtual_points(real_points_n);

    return centroid;
  }

  void add_to_centroids(const Canvas_point& p,
                        const Canvas_point& q,
                        const Canvas_point& r,
                        const std::size_t seed_id,
                        std::vector<FT>& xs, std::vector<FT>& ys,
                        std::vector<FT>& total_areas)
  {
    CGAL_precondition(seed_id < seeds.size());

    Point_2 local_centroid = CGAL::centroid(p.point, q.point, r.point);
    FT area = triangle_area_in_metric(p, q, r);
    CGAL_postcondition(area != 0.);

    xs[seed_id] += area * local_centroid.x();
    ys[seed_id] += area * local_centroid.y();
    total_areas[seed_id] += area;
  }

  void compute_centroid_of_starset_face(const typename Star_set::Face_handle fh,
                                        std::vector<FT>& xs, std::vector<FT>& ys,
                                        std::vector<FT>& total_areas)
  {
    Tri triangle;
    triangle[0] = fh->vertex(0)->info();
    triangle[1] = fh->vertex(1)->info();
    triangle[2] = fh->vertex(2)->info();

    Simplex colors;
    for(int i=0; i<3; ++i)
      colors.insert(points[triangle[i]].closest_seed_id);

    if(colors.size() == 1)
    {
      const Canvas_point& cp = points[triangle[0]];

      add_to_centroids(cp,
                       points[triangle[1]],
                       points[triangle[2]],
                       cp.closest_seed_id, xs, ys, total_areas);
    }
    else if(false && colors.size() == 2)
    {
      // see what's done in the [colors.size() == 2] case in geo_gen_on_tr
      // I don't like the way how it is split there, though.
      // Thus, always splitting in three if I can't not split.
    }
    else // colors.size() == 3
    {
      // compute the on the edges
      boost::array<Canvas_point,3> mid_pts;
      for(int i=0; i<3; ++i)
        mid_pts[i] = Vor_vertex_on_edge(triangle[i], triangle[(i+1)%3]);

      // compute roughly the center point
      const Canvas_point& v = compute_precise_Voronoi_vertex(triangle);

      const Canvas_point& cp = points[triangle[0]];
      const Canvas_point& cq = points[triangle[1]];
      const Canvas_point& cr = points[triangle[2]];

      // 'triangle' is split in 6 smaller triangles.
      // p-i1-v & p-i3-v --> associated to the seed p's closest_seed_id
      // q-i1-v & q-i2-v --> associated to the seed q's closest_seed_id
      // r-i2-v & r-i3-v --> associated to the seed r's closest_seed_id

      //            p
      //           / |
      //          /  |
      //         /   |
      //       i1    |
      //       /|    |
      //      /  |   |
      //     /   v---i3
      //    /    |   |
      // q /____i2___| r

      add_to_centroids(cp, v, mid_pts[0], cp.closest_seed_id, xs, ys, total_areas);
      add_to_centroids(cp, v, mid_pts[2], cp.closest_seed_id, xs, ys, total_areas);

      add_to_centroids(cq, v, mid_pts[0], cq.closest_seed_id, xs, ys, total_areas);
      add_to_centroids(cq, v, mid_pts[1], cq.closest_seed_id, xs, ys, total_areas);

      add_to_centroids(cr, v, mid_pts[1], cr.closest_seed_id, xs, ys, total_areas);
      add_to_centroids(cr, v, mid_pts[2], cr.closest_seed_id, xs, ys, total_areas);
    }
  }

  void shave_off_virtual_points(const std::size_t real_points_n)
  {
//    std::cout << "shaving: " << points.size() << " to " << real_points_n << std::endl;
    std::vector<Canvas_point>::iterator it = points.begin();
    std::advance(it, real_points_n);
    points.erase(it, points.end());
  }

  void set_centroid(const std::size_t seed_id, const Point_2& c,
                    std::vector<Point_2>& centroids)
  {
    Seed_status s = seeds_s[seed_id];
    const Point_2& old_centroid = seeds[seed_id];

    if(s == IN_DOMAIN)
      centroids[seed_id] = c;
    else if(s == ON_CONSTRAINT)
    {
      std::size_t vcs = Vor_cons.size();
      bool found = false;
      for(std::size_t i=0; i<vcs; ++i)
      {
        const Voronoi_constraint& vc = Vor_cons[i];
        if(!vc.has_point(old_centroid))
          continue;

        found = true;

        Vector2d old_new;
        old_new(0) = c.x() - old_centroid.x();
        old_new(1) = c.y() - old_centroid.y();

        const Vector2d& dir = vc.dir;
        FT sp = old_new.dot(dir);

        FT new_x = old_centroid.x() + sp * dir(0);
        FT new_y = old_centroid.y() + sp * dir(1);

//        std::cout << "old seed: " << old_centroid << std::endl;
//        std::cout << "constraint: " << vc.p1 << " -- " << vc.p2
//                  << " has: " << old_centroid << std::endl;
//        std::cout << "old->new: " << old_new.transpose() << std::endl;
//        std::cout << "scalar prod: " << sp << std::endl;
//        std::cout << "dir: " << dir.transpose() << std::endl;
//        std::cout << "new centroid: " << new_x << " " << new_y << std::endl;

        CGAL_postcondition(vc.has_point(Point_2(new_x, new_y)));

        centroids[seed_id] = Point_2(new_x, new_y);
        return;
      }

      if(!found)
      {
        std::cout << "ON_CONSTRAINT status but on no actual constraint..." << '\n';
        std::cout << seed_id << " at " << old_centroid << std::endl;
        exit(0);
      }
    }
    else // if (s == ON_CORNER)
    {
      centroids[seed_id] = seeds[seed_id]; // don't move corners atm.
    }
  }

  void compute_all_centroids(std::vector<Point_2>& centroids)
  {
    // compute the centroid through the sum c = sum_i (ci*area_i) / sum_i (area_i)
    // with the tiny triangles of the base mesh

    // it's not exact because we only consider the triangles with all vertices
    // having seed_id as closest seed (thus ignoring border triangles that
    // only have <= 2 vertices that have seed_id as closest_seed)

    std::size_t seeds_n = centroids.size();
    std::vector<FT> total_areas(seeds_n, 0);
    std::vector<FT> xs(seeds_n, 0);
    std::vector<FT> ys(seeds_n, 0);

    std::size_t real_points_n = points.size();

    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];

      Face_handle_handle fhit = star->finite_incident_faces_begin();
      Face_handle_handle fhend = star->finite_incident_faces_end();
      for(; fhit!=fhend; ++fhit)
      {
        typename Star_set::Face_handle fh = *fhit;

        if(is_face_outside_domain(star, fh))
          continue;

        compute_centroid_of_starset_face(fh, xs, ys, total_areas);
      }
    }

    shave_off_virtual_points(real_points_n);

    for(std::size_t i=0; i<seeds_n; ++i)
    {
      if(total_areas[i] == 0.)
      {
        std::cout << "Warning: Starset has not vertex of color: " << i << std::endl;
        centroids[i] = seeds[i];
        continue;
      }

      FT d = 1. / total_areas[i];
      FT x = xs[i] * d;
      FT y = ys[i] * d;
      Point_2 c(x, y);

      set_centroid(i, c, centroids);
    }
  }

  FT optimize_seed(const std::size_t seed_id,
                   const std::vector<Point_2>& cell_centroids)
  {
    std::cout << "optimize seed: " << seed_id << std::endl;

    const Point_2 old_seed = seeds[seed_id];
    const Point_2& new_seed = cell_centroids[seed_id];

    seeds[seed_id] = new_seed;
    seeds_m[seed_id] = mf->compute_metric(new_seed);

    // squared displacement between the old and the new seed
    K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();
    FT displacement = sqd(old_seed, new_seed);
    std::cout << "seed n: " << seed_id << " ----- ";
    std::cout << "old/new: " << old_seed << " || " << new_seed << " ";
    std::cout << " displacement: " << displacement << std::endl;

    return displacement;
  }

  FT compute_CVT_energy() const
  {
    FT e = 0.;
    FT third = 1./3.;

    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];

      Face_handle_handle fhit = star->finite_incident_faces_begin();
      Face_handle_handle fhend = star->finite_incident_faces_end();
      for(; fhit!=fhend; ++fhit)
      {
        typename Star_set::Face_handle fh = *fhit;

        if(is_face_outside_domain(star, fh))
          continue;

        const Canvas_point& cp = points[fh->vertex(0)->info()];
        const Canvas_point& cq = points[fh->vertex(1)->info()];
        const Canvas_point& cr = points[fh->vertex(2)->info()];

        FT dist = third * (cp.distance_to_closest_seed +
                           cq.distance_to_closest_seed +
                           cr.distance_to_closest_seed);
        FT area = triangle_area_in_metric(cp, cq, cr);
        e += area * dist * dist;
      }
    }

    return e;
  }

  void optimize_seeds()
  {
    std::cout << "optimize seeds" << std::endl;

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
      if(points[i].state != KNOWN)
        std::cout << "Warning: point " << i << " is not KNOWN at Opti" << std::endl;

    FT sq_bbox_diag_l = sq_bbox_diagonal_length(base_mesh_bbox);
    bool is_optimized = false;
    int counter = 0;

    do
    {
      // compute the centroid of every face of the star set
      std::vector<Point_2> cell_centroids(seeds.size());
      compute_all_centroids(cell_centroids);

      // Optimize each seed
      FT cumulated_displacement = 0;
      for(std::size_t i=0, ss=seeds.size(); i<ss; ++i)
        cumulated_displacement += optimize_seed(i, cell_centroids);

      FT e = compute_CVT_energy();

      std::cout << "at : " << counter << ", cumulated displacement : "
                << cumulated_displacement << std::endl;
      std::cout << "energy at : " << e << std::endl;

      reset();
      locate_and_initialize_seeds();
      spread_distances(true/*use_dual_shenanigans*/);

      if(counter%1 == 0)
      {
        std::ostringstream opti_out;
        opti_out << "optimized_" << str_base_mesh << "_starset_" << counter << std::ends;
        output_canvas_data_and_dual(opti_out.str().c_str());
      }

      is_optimized = (++counter >= max_opti_n); // (cumulated_displacement < sq_bbox_diag_l * 1e-5);
      if(cumulated_displacement == 0.0)
        break;
    }
    while(!is_optimized);
  }

  // Dual stuff ----------------------------------------------------------------
  void add_simplex_to_triangulation(const Simplex& dual_simplex)
  {
    CGAL_assertion(dual_simplex.size() <= 3);

    Simplex::const_iterator it = dual_simplex.begin();
    if(dual_simplex.size() == 2) // an edge!
    {
      Edge e; e[0] = *it; e[1] = (*++it);
      dual_edges.insert(e);
    }
    else if(dual_simplex.size() == 3) // a triangle!
    {
      Tri tr; tr[0] = *it; tr[1] = (*++it); tr[2] = (*++it);
      dual_triangles.insert(tr);

      add_triangles_edges_to_dual_edges(tr);
    }
  }

  void add_simplex_to_triangulation(const Canvas_point& cp,
                                    const Simplex& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    std::size_t length = cp.ancestor_path_length;
    if(length < min_ancestor_path_length)
    {
//      std::cout << "the canvas is thin (" << length << ") at : " << cp.index << " ("
//                << cp.point << ") dual simplex of size: " << dual_simplex.size()
//                << " and cellid: " << cp.closest_seed_id << std::endl;
    }

    // test against criteria should be moved somewhere else
    // than during the dual computations todo
    if(dual_simplex.size() == 3
       || (dual_simplex.size() == 2 && cp.is_on_domain_border)
       )
      test_simplex(cp, dual_simplex);

    add_simplex_to_triangulation(dual_simplex);
  }

  void add_simplex_to_triangulation(const typename Star_set::Face_handle fh,
                                    const Simplex& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    int max_dist_p = -1;
    FT max = -FT_inf;
    for(int i=0; i<3; ++i)
    {
      const Canvas_point& cp = points[fh->vertex(i)->info()];
      const FT d = cp.distance_to_closest_seed;

      if(d > max)
      {
        max = d;
        max_dist_p = fh->vertex(i)->info();
      }
    }

    // keep as dual point the farthest canvas point -- at least until we have
    // precise Voronoi computations)
    return add_simplex_to_triangulation(points[max_dist_p], dual_simplex);
  }

  void compute_dual()
  {
    if(!dual_edges.empty() || !dual_triangles.empty())
      return;

    double duration;
    std::clock_t d_start = std::clock();

#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0, ts=ss.size(); i<ts; ++i)
    {
      const Star_handle star = ss[i];

      Face_handle_handle fhit = star->finite_incident_faces_begin();
      Face_handle_handle fhend = star->finite_incident_faces_end();
      for(; fhit!=fhend; ++fhit)
      {
        Face_handle fh = *fhit;
        if(is_face_outside_domain(star, fh))
          continue;

        std::size_t i0 = fh->vertex(0)->info();
        std::size_t i1 = fh->vertex(1)->info();
        std::size_t i2 = fh->vertex(2)->info();

        // filter triangles with only one color
        if(points[i0].closest_seed_id == points[i1].closest_seed_id &&
           points[i1].closest_seed_id == points[i2].closest_seed_id)
          continue;

#ifndef COMPUTE_DUAL_FOR_ALL_DIMENSIONS
        if(points[i0].closest_seed_id != points[i1].closest_seed_id &&
           points[i1].closest_seed_id != points[i2].closest_seed_id &&
           points[i0].closest_seed_id != points[i2].closest_seed_id)
        { /* don't filter the case: 3 different colors */ }
        else
        {
          // filter triangles with _exactly_ two colors and none of the vertices
          // are on the border of the domain
          if(!points[i0].is_on_domain_border && !points[i1].is_on_domain_border &&
             !points[i2].is_on_domain_border)
            continue;
        }
#endif
        Simplex dual_simplex;
        for(int j=0; j<3; ++j)
          dual_simplex.insert(points[fh->vertex(j)->info()].closest_seed_id);

        add_simplex_to_triangulation(fh, dual_simplex);
      }
    }

    duration = ( std::clock() - d_start ) / (double) CLOCKS_PER_SEC;
    std::cout << "End dual computations: " << duration << std::endl;
  }

  // output stuff --------------------------------------------------------------
  void output_straight_dual(const std::string str_base)
  {
    if(dual_edges.empty() && dual_triangles.empty())
      compute_dual();
    else
      std::cout << "dual already computed" << std::endl;

#if (verbose > 5)
    std::cout << "captured: ";
    std::cout << seeds.size() << " points, ";
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

    for(boost::unordered_set<Tri>::iterator it = dual_triangles.begin();
                                            it != dual_triangles.end(); ++it)

    {
      const Tri& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
      out << is_triangle_intersected(tr, dual_edges) << std::endl;
    }
    out << "End" << std::endl;
  }

  void output_canvas_data_and_dual(const std::string str_base)
  {
    output_canvas(str_base);
    output_straight_dual(str_base);
  }

  Base_mesh(Star_set& ss_)
    :
      ss(ss_),
      points(),
      trial_points(),
      dual_edges(),
      dual_triangles(),
      Vor_cons(),
      size_queue(this),
      distortion_queue(this),
      intersection_queue(this),
      quality_queue(this),
      Voronoi_vertices()
  {
    initialize_from_starset();
  }
};

void Canvas_point::change_state(FMM_state new_state)
{
  if(new_state == state)
    std::cout << "WARNING: useless state change..." << std::endl;

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

std::size_t Canvas_point::anc_path_length() const
{
  std::size_t i = 0;
  std::size_t n_anc = ancestor;
  while(n_anc != static_cast<std::size_t>(-1))
  {
    ++i;
    const Canvas_point& anc = bm->points[n_anc];
    n_anc = anc.ancestor;

    CGAL_assertion(i < bm->ss.size());
  }
  return i;
}

void Canvas_point::print_children() const
{
  std::cout << "point : " << index << " children: ";

  Point_set::iterator cit = children.begin();
  Point_set::iterator cend = children.end();
  for(; cit!=cend; ++cit)
  {
    CGAL_assertion(*cit < bm->points.size());
    std::cout << *cit;
  }
  std::cout << std::endl;
}

void Canvas_point::print_ancestors() const
{
  CGAL_expensive_assertion_code(anc_path_length()); // assert non circular ancestry
  std::cout << "point : " << index << " ancestors: ";

  std::size_t anc = ancestor;
  while(anc != static_cast<std::size_t>(-1))
  {
    CGAL_assertion(anc < bm->points.size());
    std::cout << anc << " ";
    anc = bm->points[anc].ancestor;
  }
  std::cout << std::endl;
}

bool Canvas_point::find_in_ancestry(const std::size_t i) const
{
  CGAL_assertion(ancestor == static_cast<std::size_t>(-1) ||
                 anc_path_length()); // assert non circular ancestry
  std::size_t anc = ancestor;
  while(anc != static_cast<std::size_t>(-1))
  {
    if(anc == i)
      return true;
    anc = bm->points[anc].ancestor;
  }
  return false;
}

void Canvas_point::remove_from_children(const std::size_t c)
{
  // 99% of the time, the child is found, but if during the initialization of
  // a new seed, you reset (at least) two points that have an ancestry relationship
  // then the child could have been reset already...
  Point_set::iterator it = children.find(c);
  if(it != children.end())
    children.quick_erase(it);
//  else
//  {
//    std::cout << "WARNING: call to remove_from_children didn't find the child (";
//    std::cout << c << " from " << index << ")" << std::endl;
//  }
}

void Canvas_point::mark_descendants(std::size_t& count)
{
//  std::cout << "marking " << index << " (a: " << ancestor << ") and its "
//            << children.size() << " children: ";
//  Point_set::iterator cit = children.begin();
//  Point_set::iterator end = children.end();
//  for(; cit!=end; ++cit)
//    std::cout << *cit << " ";
//  std::cout << std::endl;

  ++count;
  state = ORPHAN;

  Point_set::iterator cit = children.begin();
  Point_set::iterator end = children.end();
  for(; cit!=end; ++cit)
  {
    Canvas_point& cp = bm->points[*cit];
    if(cp.state == ORPHAN)
      continue;
    cp.mark_descendants(count);
  }
}

void Canvas_point::reset_descendants()
{
//  std::cout << "reset descendant at: " << index << std::endl;
  CGAL_precondition(state == ORPHAN);

  while(!children.empty())
  {
    Canvas_point& cp = bm->points[*(children.begin())];
    cp.reset_descendants();
  }

  // we might reset points that are in the priority queue here !
  // but it's okay since get_next_trial_point() will ignore these
  reset();
  state = ORPHAN;
}

void Canvas_point::deal_with_descendants()
{
//  std::cout << "lineage starting at " << index << std::endl;

  // We recompute all the values at these orphans only using neighbors that are
  // NOT orphans. We start from the youngest child and climb up in the tree.
  // We leave these new values as FAR.

  // We mark all the lineage as orphans, we don't want them to be relevant to
  // any computation anymore since they're not linked to a seed by an ancestor
  // tree
  std::size_t count = 0;
  mark_descendants(count);
  reset_descendants();

//  std::cout <<  count << " children in the lineage" << std::endl;
}

bool Canvas_point::compute_closest_seed(const std::size_t n_anc,
                                      const bool overwrite)
{
#if (verbose > 20)
  std::cout << "closest seed at " << index << " from " << n_anc;
  std::cout << " (current d/a: " << distance_to_closest_seed << " " << ancestor;
  std::cout << " seed: " << closest_seed_id << ")" << std::endl;
#endif

  CGAL_precondition(static_cast<std::size_t>(n_anc) < bm->points.size());
  const Canvas_point& anc = bm->points[n_anc];
  if(anc.state != KNOWN)
    std::cout << "WARNING: potential ancestor is not KNOWN" << std::endl;

  if(anc.find_in_ancestry(index))
    return false;

  FT d = FT_inf;

  // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
  // stored as 'this', ancestor1, ancestor2, etc.
  boost::array<std::size_t, k+1> ancestor_path;
  for(int i=1; i<k+1; ++i)
    ancestor_path[i] = -1;
  ancestor_path[0] = this->index;
  CGAL_assertion(ancestor_path[0] != static_cast<std::size_t>(-1));

  std::size_t n_curr_ancestor = n_anc;
  std::size_t best_i = -1;
  for(int i=1; i<=k; ++i)
  {
    // add the new segment to the ancestor path
    ancestor_path[i] = n_curr_ancestor;
    const Canvas_point& curr_ancestor = bm->points[n_curr_ancestor];

    Vector2d ancestor_edge;
    ancestor_edge(0) = point.x() - curr_ancestor.point.x();
    ancestor_edge(1) = point.y() - curr_ancestor.point.y();
    FT ancestor_edge_length = ancestor_edge.norm();
    Vector2d normalized_anc_edge = ancestor_edge/ancestor_edge_length;

    // compute the distance for the current depth (i)
    FT dist_to_ancestor = 0.;
    for(int j=0; j<i; ++j) // we add a part for each edge in the path
    {
      // get the metric for the current edge

      CGAL_assertion(ancestor_path[j] < bm->points.size() &&
                     ancestor_path[j+1] < bm->points.size());
      const Canvas_point& e0 = (j==0)?*this:bm->points[ancestor_path[j]];
      const Canvas_point& e1 = bm->points[ancestor_path[j+1]];

      const Metric& m0 = e0.metric;
      const Metric& m1 = e1.metric;

      Vector2d curr_edge;
      curr_edge(0) = e0.point.x() - e1.point.x();
      curr_edge(1) = e0.point.y() - e1.point.y();

      // interpolate between both metric and transform the normalized edge
      // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
      const Eigen::Matrix2d& f = get_interpolated_transformation(m0, m1);
      Vector2d transformed_curr_edge = f*normalized_anc_edge;

      FT sp = curr_edge.dot(normalized_anc_edge);
      FT l = transformed_curr_edge.norm(); // length of the normalized anc edge in the metric

      dist_to_ancestor += sp * l;
    }
    dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

    // add ancestor edge length to the distance at that ancestor
    FT dist_at_anc = curr_ancestor.distance_to_closest_seed;
    FT new_d = dist_at_anc + dist_to_ancestor;

    if(new_d < d)
    {
      best_i = i;
      d = new_d;
    }

    // check if we can go any farther up in the ancestor tree
    if(curr_ancestor.ancestor == static_cast<std::size_t>(-1))
      break;

    n_curr_ancestor = curr_ancestor.ancestor;
  }

  ++(optimal_depths[best_i]);

  if(d < distance_to_closest_seed)
  {
#if (verbose > 20)
    std::cout << "improving " << distance_to_closest_seed << " from " << ancestor
              << " (" << closest_seed_id << ")";
    std::cout << " to " << d << " from " << n_anc << " (" << anc.closest_seed_id
              << ")" << std::endl;
#endif
    if(overwrite)
    {
      // remove 'index' from the previous ancestor's children (if needed)
      if(ancestor != static_cast<std::size_t>(-1))
        bm->points[ancestor].remove_from_children(index);

      ancestor = n_anc;
      depth = best_i;
      distance_to_closest_seed = d;
      closest_seed_id = anc.closest_seed_id;
      ancestor_path_length = anc.ancestor_path_length + 1;

      // add 'index' to the new ancestor's children
      bm->points[ancestor].children.insert(index);
//      std::cout << "insert " << index << " in " << ancestor << " 's children" << std::endl;
//      std::cout << ancestor << " children size : " << bm->points[ancestor].children.size() << std::endl;

      CGAL_postcondition(anc_path_length() == ancestor_path_length); // checks for circular ancestry

#ifndef USE_FULL_REBUILDS // if we're using full rebuilds, below doesn't matter

      // If we're refining, the point 'this' might have children.
      // If 'this' changes, then any value at the children becomes meaningless.
      // If we just ignore them the values are usually overwritten when we call
      // this->update_neighbors()... BUT in very rare cases, they do NOT and
      // we get orphans (and we don't like orphans), so when the parent changes,
      // all the children must be reset to be sure that they'll always get a
      // valid parent.

      // Last but not least, this new parent is not necessarily 'this', so we must
      // consider all potential new parents (who might even potentially have a
      // different color than 'this')

//      if(!children.empty())
//      {
//        std::cout << "dealing with the " << children.size()
//                  << " descendant(s) of " << index << std::endl;
//        std::cout << "position : " << point << std::endl;
//      }

      if(ignore_children)
        return true;

      while(!children.empty())
      {
        Canvas_point& cp = bm->points[*(children.begin())];
        cp.deal_with_descendants();
        CGAL_postcondition(cp.anc_path_length() == cp.ancestor_path_length); // checks for circular ancestry
      }
#endif
    }
    return true;
  }
  return false;
}

void Canvas_point::update_neighbor_distance(Canvas_point& cp,
                                          std::vector<std::size_t>& trial_pq,
                                          PQ_state& pqs_ret) const
{
#if (verbose > 12)
  std::cout << "neighbor: " << cp.index << " has state: " << cp.state << std::endl;
#endif

  if(cp.state == KNOWN)
    return;
  else if(cp.state == TRIAL)
  {
    // note that we don't insert in trial_pq since it's already in
    if(cp.compute_closest_seed(index))
      pqs_ret = REBUILD_TRIAL;
  }
  else
  {
    CGAL_assertion(cp.state == FAR || cp.state == ORPHAN);
    if(cp.compute_closest_seed(index))
    {
      cp.state = TRIAL;
      trial_pq.push_back(cp.index);
      std::push_heap(trial_pq.begin(), trial_pq.end(),
                     Canvas_point_comparer<Base_mesh>(bm));
    }
  }
}

PQ_state Canvas_point::update_neighbors_distances(std::vector<std::size_t>& trial_pq) const
{
  // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (verbose > 10)
  std::cout << "update neighbors of " << index << std::endl;
#endif
  CGAL_assertion(state == KNOWN);

  PQ_state pqs_ret = NOTHING_TO_DO;

  Star_handle star = bm->ss[index];
  Vertex_handle_handle vit = star->finite_adjacent_vertices_begin();
  Vertex_handle_handle vend = star->finite_adjacent_vertices_end();
  for(; vit!=vend; ++vit)
  {
    Vertex_handle vh = *vit;

    if(is_a_corner_vertex(vh))
      continue;

    Canvas_point& cp = bm->points[vh->info()];
    update_neighbor_distance(cp, trial_pq, pqs_ret);
  }

#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
  Point_list::const_iterator lit = reciprocal_neighbors.begin();
  Point_list::const_iterator lend = reciprocal_neighbors.end();
  for(; lit!=lend; ++lit)
  {
    if(is_a_corner_vertex(*lit))
      continue;

    Canvas_point& cp = bm->points[*lit];
    update_neighbor_distance(cp, trial_pq, pqs_ret);
  }
#endif

  return pqs_ret;
}

FT Canvas_point::distortion_to_seed() const
{
  FT gamma = 1.;
  //    std::cout << "init gamma: " << gamma << " " << index << std::endl;
  const Canvas_point* curr = this;
  std::size_t n_anc = ancestor;

  while(n_anc != static_cast<std::size_t>(-1))
  {
    const Canvas_point& anc = bm->points[n_anc];
    const Metric& m1 = anc.metric;
    const Metric& m2 = curr->metric;
    FT loc_gamma = m1.compute_distortion(m2);
    //      std::cout << "loc: " << loc_gamma << " " << anc->index << std::endl;
#if 1
    gamma *= loc_gamma;
#else
    gamma = (std::max)(loc_gamma, gamma);
#endif
    //      std::cout << "gamma: " << gamma << std::endl;
    n_anc = anc.ancestor;
    curr = &anc;
  }

  gamma = (std::min)(25., gamma);

  //    std::cout << "final gamma:" << gamma << std::endl;
  //    exit(0);
  return gamma;
}

void Canvas_point::reset()
{
  // this only resets 'color'-related fields
#if (verbose > 30)
  std::cout << "reset() at " << index << std::endl;
#endif

  if(ancestor != static_cast<std::size_t>(-1))
    bm->points[ancestor].remove_from_children(index);

  Point_set::iterator cit = children.begin();
  Point_set::iterator end = children.end();
  for(; cit!=end; ++cit)
    bm->points[*cit].ancestor = -1;

  state = FAR;
  is_seed_holder = -1;
  distance_to_closest_seed = FT_inf;
  closest_seed_id = -1;
  ancestor = -1;
  depth = 0;
  children.clear();
  is_Voronoi_vertex = false;
  ancestor_path_length = 0;
}

void Canvas_point::initialize_from_point(const FT d,
                                       const std::size_t seed_id)
{
  reset();

  closest_seed_id = seed_id;
  is_seed_holder = seed_id;
  distance_to_closest_seed = d;
  state = TRIAL;
  ancestor = -1;
}
bool Canvas_point::operator==(const Canvas_point& cp) const
{
  // the only thing that really matters is the point
  return (point == cp.point);
}

Canvas_point::Canvas_point()
  :
    bm(NULL),
    point(),
    index(static_cast<std::size_t>(-1)),
    metric(),
    is_on_domain_border(false),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(),
#endif
    state(FAR),
    is_seed_holder(-1),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    depth(0),
    children(),
    is_Voronoi_vertex(false),
    ancestor_path_length(0)
{ }

Canvas_point::Canvas_point(Base_mesh * bm_,
                           const std::size_t index_,
                           const bool is_on_domain_border_)
  :
    bm(bm_),
    point(bm->ss[index_]->center_point()),
    index(index_),
    metric(bm->ss[index_]->metric()),
    is_on_domain_border(is_on_domain_border_),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(),
#endif
    state(FAR),
    is_seed_holder(-1),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    depth(0),
    children(),
    is_Voronoi_vertex(false),
    ancestor_path_length(0)
{

}

Canvas_point::Canvas_point(Base_mesh * bm_,
                           const Point_2& p,
                           const std::size_t index_,
                           const Metric& metric_,
                           const bool is_on_domain_border_)
  :
    bm(bm_),
    point(p),
    index(index_),
    metric(metric_),
    is_on_domain_border(is_on_domain_border_),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(),
#endif
    state(FAR),
    is_seed_holder(-1),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    depth(0),
    children(),
    is_Voronoi_vertex(false),
    ancestor_path_length(0)
{ }

Canvas_point::Canvas_point(const Canvas_point& cp)
  :
    bm(cp.bm),
    point(cp.point),
    index(cp.index),
    metric(cp.metric),
    is_on_domain_border(cp.is_on_domain_border),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(cp.reciprocal_neighbors),
#endif
    state(cp.state),
    is_seed_holder(cp.is_seed_holder),
    distance_to_closest_seed(cp.distance_to_closest_seed),
    closest_seed_id(cp.closest_seed_id),
    ancestor(cp.ancestor),
    depth(cp.depth),
    children(cp.children),
    is_Voronoi_vertex(cp.is_Voronoi_vertex),
    ancestor_path_length(cp.ancestor_path_length)
{ }

struct Point_extracter
    : public std::unary_function<Canvas_point, Point_2>
{
  const Point_2 operator()(const Canvas_point& cp) const
  {
    return cp.point;
  }
};

void draw_metric_field(const MF mf, const Base_mesh& bm)
{
  typedef boost::transform_iterator<Point_extracter,
      std::vector<Canvas_point>::const_iterator> Extracted_iterator;
  mf->draw(Extracted_iterator(bm.points.begin(), Point_extracter()),
           Extracted_iterator(bm.points.end(), Point_extracter()));
}

} // aniso_mesh_2
} // cgal

using namespace CGAL::Anisotropic_mesh_2;

int main(int, char**)
{
  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);
  double duration;
  start = std::clock();
  std::srand(0);

//  mf = new Euclidean_metric_field<K>(0.5, 2.5);
  mf = new Custom_metric_field<K>();
//  mf = new External_metric_field<K>("freefem.mesh", "freefem.sol");

  Star_set starset;
//  generate_starset(starset);
  read_dump(starset);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "Generated starset: " << duration << std::endl;

  Base_mesh bm(starset);
  bm.output_canvas("base_starset");

//  draw_metric_field(mf, bm);
//  exit(0);

  bm.initialize_seeds();
  bm.locate_and_initialize_seeds();
  bm.spread_distances(true/*use_dual_shenanigans*/);

  if(n_refine > 0)
  {
    bm.output_canvas_data_and_dual(str_base_mesh + "_pre_ref");

    for(int i=0; i<n_refine; ++i)
    {
      // can't compute the dual while spreading if we're refining since we have already a layer
      // of paint laying on the canvas...
      bool successful_insert = bm.refine_seeds_with_self_computed_ref_point();

      if(i%10 == 0)
      {
        std::ostringstream out;
        out << "ref_" << seeds.size();
        bm.output_straight_dual(out.str());
        if(i%100 == 0)
          bm.output_canvas(out.str());
      }

      if(!successful_insert)
        break;
    }

    for(std::size_t i=0, ps=bm.points.size(); i<ps; ++i)
      bm.points[i].state = KNOWN;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "End refinement: " << duration << std::endl;
  }

  bm.output_canvas_data_and_dual(str_base_mesh + "_starset");

  // optimize stuff
  if(max_opti_n > 0)
  {
    bm.optimize_seeds();
    bm.output_canvas_data_and_dual("optimized_" + str_base_mesh + "_starset");
  }

  ignore_children = true;
//  bm.compute_geodesics(str_base_mesh);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "duration: " << duration << std::endl;
}
