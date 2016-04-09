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

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include <iostream>
#include <functional>
#include <limits>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <utility>

#define ANISO_USE_RECIPROCAL_NEIGHBORS
// #define COMPUTE_PRECISE_VOR_VERTICES
// #define USE_FULL_REBUILDS
#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 6

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
typedef K::Vector_2                                           Vector;
typedef K::Triangle_2                                         Triangle_2;
typedef K::Segment_3                                          Segment_3;
typedef K::Triangle_3                                         Triangle_3;

typedef KExact::Point_2                                       EPoint;
typedef KExact::Segment_2                                     ESegment;
typedef KExact::Triangle_2                                    ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                  To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                  Back_from_exact;

const FT FT_inf = std::numeric_limits<FT>::infinity();
To_exact to_exact;
Back_from_exact back_from_exact;

//typedef CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
typedef CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;
//typedef CGAL::Anisotropic_mesh_2::External_metric_field<K>* MF;
MF mf;

//geometry
FT a = 2; // half the side (x)
FT b = 2; // half the side (y)
Point_2 center(0., 0.);
Rectangle_domain<K> pdomain(a, b, center);

//face criteria
FT f_r0 = 0.1; // since geodesic takes r0=1, we need r_0 << 1 here (0.1-ish ?)
FT f_rho0 = 3.0;

//misc
FT gamma = 1.0;
FT beta = 2.5;
FT delta = 0.3;
int max_times_to_try = 60;
int nb = 20; // nb_initial_points

Criteria_base<K> criteria(f_r0, f_rho0, gamma, beta, delta, nb, max_times_to_try);

void generate_starset(Star_set& ss)
{
  Anisotropic_mesher_2<K> mesher(ss, &pdomain, &criteria, mf);
  mesher.refine_mesh();

  std::ofstream out("geo_starset.mesh");
  output_medit(ss, out, false/*consistent_only*/);
  std::ofstream out_dump("dump.mesh"); // useful to quickly recreate a starset
  dump(ss, out_dump);
}

void read_dump(Star_set& ss)
{
  std::cout << "Reading dump..." << std::endl;
  std::ifstream in("dump.mesh");
  ss.clear();

  std::size_t stars_n, v_n, id;
  FT x, y;

  in >> stars_n;

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

void refine_starset(Star_set& ss)
{

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

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 5;
static const int k = 8; // depth of the ancestor edge

// the metric field and the seeds
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

// seeds
const std::string str_seeds = "ref_3052_dual";

// base mesh
const std::string str_base_mesh = "geo_starset_dump";
CGAL::Bbox_2 base_mesh_bbox;

// refinement
int n_refine = 1000;
std::size_t min_ancestor_path_length = 5;

//debug & info
int known_count=0, trial_count=0, far_count=0;
std::clock_t start;
std::vector<int> optimal_depths(k, 0);

typedef boost::unordered_map<boost::array<std::size_t, k+1>, FT> Path_map;
Path_map computed_paths;

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
    for(boost::unordered_set<Edge>::const_iterator it = edges.begin();
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

FT triangle_area(const Point_2& p, const Point_2& q, const Point_2& r)
{
   K::Compute_area_2 o;
  return std::abs(o(p, q, r));

  FT pq_x = q.x() - p.x();
  FT pq_y = q.y() - p.y();
  FT sq_pq = pq_x*pq_x + pq_y*pq_y;

  FT pr_x = r.x() - p.x();
  FT pr_y = r.y() - p.y();
  FT sq_pr = pr_x*pr_x + pr_y*pr_y;

  FT qr_x = q.x() - r.x();
  FT qr_y = q.y() - r.y();
  FT sq_qr = qr_x*qr_x + qr_y*qr_y;

  FT p1 = sq_pq + sq_pr + sq_qr;
  FT sq_p1 = p1*p1;

  FT qu_pq = sq_pq*sq_pq;
  FT qu_pr = sq_pr*sq_pr;
  FT qu_qr = sq_qr*sq_qr;

  FT p2 = qu_pq + qu_pr + qu_qr;

  FT diff = sq_p1 - 2*p2;
  CGAL_assertion(diff > -1e-16);
  diff = (std::max)(0., diff);

  FT area = 0.25 * std::sqrt(diff);
  return area;
}

FT triangle_area_in_metric(const Point_2& p, const Point_2& q, const Point_2& r)
{
  // compute the interpolated's metric transformation
  FT third = 1./3.;

  // very unefficient fixme
  Eigen::Matrix2d f = third*(mf->compute_metric(p).get_transformation() +
                             mf->compute_metric(q).get_transformation() +
                             mf->compute_metric(r).get_transformation());

  // transform the points
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

void compute_bary_weights(const Point_2&p , const Point_2& a, const Point_2& b, const Point_2& c,
                          FT& lambda_a, FT& lambda_b, FT& lambda_c)
{
  CGAL_assertion(!(Triangle_2(a,b,c)).is_degenerate());
  Vector2d v_ab, v_ac, v_ap;
  v_ab(0) = b.x() - a.x();
  v_ab(1) = b.y() - a.y();

  v_ac(0) = c.x() - a.x();
  v_ac(1) = c.y() - a.y();

  v_ap(0) = p.x() - a.x();
  v_ap(1) = p.y() - a.y();

  FT d11 = v_ab.dot(v_ab);
  FT d12 = v_ab.dot(v_ac);
  FT d22 = v_ac.dot(v_ac);
  FT d31 = v_ap.dot(v_ab);
  FT d32 = v_ap.dot(v_ac);
  CGAL_assertion(d11 * d22 - d12 * d12 != 0.);
  FT den = 1. / (d11 * d22 - d12 * d12);
  lambda_b = (d22 * d31 - d12 * d32) * den;
  lambda_c = (d11 * d32 - d12 * d31) * den;
  lambda_a = 1. - lambda_b - lambda_c;

//  Triangle t(a,b,c);
//  Triangle ta(b,c,p), tb(a,c,p), tc(a,b,p);
//  FT asd = ta.area() / t.area();
//  FT bsd = tb.area() / t.area();
//  FT csd = tc.area() / t.area();

//  FT check_x = lambda_a*a.x() + lambda_b*b.x() + lambda_c*c.x();
//  FT check_y = lambda_a*a.y() + lambda_b*b.y() + lambda_c*c.y();

//  CGAL_assertion(std::abs(lambda_a*a.x() + lambda_b*b.x() + lambda_c*c.x() - p.x()) < 1e-10 &&
//                 std::abs(lambda_a*a.y() + lambda_b*b.y() + lambda_c*c.y() - p.y()) < 1e-10);
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

struct Grid_point
{
  typedef boost::unordered_set<std::size_t>                           Point_set;
  typedef std::list<std::size_t>                                      Point_list;

  // immuable stuff
  Base_mesh* bm;
  std::size_t index;
  bool is_on_domain_border;

#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
  Point_list reciprocal_neighbors;
#endif

  // stuff that depends on the seeds
  FMM_state state;
  FT distance_to_closest_seed;
  std::size_t closest_seed_id;
  std::size_t ancestor;
  Point_set children;
  bool is_Voronoi_vertex;
  std::size_t ancestor_path_length;

  const Point_2& point() const;
  const Metric& metric() const;
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
  void update_neighbor_distance(Grid_point& gp,
                                std::vector<std::size_t> & trial_pq,
                                PQ_state& pqs_ret) const;
  PQ_state update_neighbors_distances(std::vector<std::size_t> & trial_pq) const;
  FT distortion_to_seed() const;
  void reset();
  void initialize_from_point(const FT d, const std::size_t seed_id);

  bool operator==(const Grid_point& gp) const;

  Grid_point();
  Grid_point(Base_mesh* bm_,
             const std::size_t index_,
             const bool is_on_domain_border_ = false);
  Grid_point(const Grid_point& gp);
};

template<typename BM>
struct Grid_point_comparer
{
  BM const * const bm;

  bool operator()(const std::size_t i1, const std::size_t i2)
  {
    const Grid_point& gp1 = bm->points[i1];
    const Grid_point& gp2 = bm->points[i2];
    return gp1.distance_to_closest_seed > gp2.distance_to_closest_seed;
  }

  Grid_point_comparer(BM const * const bm_) : bm(bm_) { }
};

int insert_new_seed(const FT x, const FT y)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < base_mesh_bbox.xmin() || x > base_mesh_bbox.xmax() ||
       y < base_mesh_bbox.ymin() || y > base_mesh_bbox.ymax())
    {
#if (verbose > 1)
      std::cout << "filtered : " << x << " " << y << std::endl;
#endif
      return seeds.size();
    }
#endif
#if (verbose > 0)
  std::cout << "added new seed: " << x << " " << y;
  std::cout << " (" << seeds.size() << ")" << std::endl;
#endif
  seeds.push_back(Point_2(x, y));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in((str_seeds + ".mesh").c_str());
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "seeds nv: " << nv << std::endl;
  CGAL_assertion(dim == 2);

  std::size_t min_nv = (std::min)(nv, vertices_nv);

  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> useless;
    insert_new_seed(r_x, r_y);

//    insert_new_seed(center.x()-a, center.x()-b);
//    insert_new_seed(center.x()+a, center.x()-b);
//    insert_new_seed(center.x()-a, center.x()+b);
//    insert_new_seed(center.x()+a, center.x()+b);

    if(seeds.size() >= vertices_nv)
      break;
  }
#if (verbose > 0)
  std::cout << "seeds: " << seeds.size() << std::endl;
#endif
  return seeds.size();
}

inline std::ostream& operator<<(std::ostream& os,
                                const boost::tuple<Simplex, std::size_t, FT>& pqe)
{
  std::cout << "PQ entry element:" << std::endl;
  std::cout << "Simplex: " << pqe.get<0>();
  std::cout << "gp: " << pqe.get<1>();
  std::cout << "val: " << pqe.get<2>() << std::endl;

  return os;
}

template<typename PQ_entry, typename Base_mesh>
struct PQ_entry_comparator
{
  Base_mesh* bm;

  bool operator()(const PQ_entry& left, const PQ_entry& right) const
  {
    const Grid_point& gpl = bm->points[left.get<1>()];
    const Grid_point& gpr = bm->points[right.get<1>()];

    if(left.get<2>() == right.get<2>()) // that's the value
    {
      if(gpl.distance_to_closest_seed ==
         gpr.distance_to_closest_seed) // that's the dist to the ref point
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
  typedef std::vector<Grid_point*>                           Grid_point_vector;

  typedef std::vector<Grid_point>                            Voronoi_vertices_container;
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

  std::vector<Grid_point> points;
  std::vector<std::size_t> trial_points;

  // Duality
  mutable boost::unordered_set<Edge> dual_edges;
  mutable boost::unordered_set<Tri> dual_triangles;

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
    std::cout << " far: " << far_count << std::endl;
  }

  void initialize_from_starset()
  {
    base_mesh_bbox = CGAL::Bbox_2();

    // build the points from the starset
    for(std::size_t i=0, sss=ss.size(); i<sss; ++i)
    {
      Star_handle star = ss[i];
      bool is_on_border = has_a_corner_vertex(star);
      Grid_point gp(this, i, is_on_border);
      points.push_back(gp);

      base_mesh_bbox += points[i].point().bbox();

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

    build_aabb_tree();
  }

  void initialize_grid_point(std::size_t id,
                             const FT dist,
                             const std::size_t seed_id)
  {
    Grid_point& gp = points[id];
    if(gp.closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cout << "WARNING: a new seed is overwriting the closest seed id at ";
      std::cout << gp.index << ". Previous index is: " << gp.closest_seed_id;
      std::cout << " (" << seeds.size() << " seeds)" << std::endl;
    }

    if(gp.state == TRIAL)
    {
      // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");
    }

    gp.initialize_from_point(dist, seed_id);
    trial_points.push_back(gp.index);
    std::push_heap(trial_points.begin(), trial_points.end(),
                   Grid_point_comparer<Self>(this));
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
    initialize_grid_point(min_id, std::sqrt(min_d), seed_id);
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
      // switch to something expensive... take the closest point...
      find_and_initialize_closest_vertex(s, seed_id);
      found = true;
    }

    if(intersections.size() > 1)
      std::cout << "more than one intersection in the locate" << std::endl;

    std::list<Primitive_Id>::const_iterator it = intersections.begin(),
                                            end = intersections.end();
    for(; it!=end; ++it)
    {
      const typename Star_set::Face_handle fh = *it;
      const Triangle_2 tr_2(points[fh->vertex(0)->info()].point(),
                            points[fh->vertex(1)->info()].point(),
                            points[fh->vertex(2)->info()].point());

      if(K().has_on_bounded_side_2_object()(tr_2, s) ||
         K().has_on_boundary_2_object()(tr_2, s))
      {
        found = true;
#if (verbose > 5)
        std::size_t i0 = fh->vertex(0)->info();
        std::size_t i1 = fh->vertex(1)->info();
        std::size_t i2 = fh->vertex(2)->info();

        std::cout << "locating seed " << seed_id
                  << " point: " << s.x() << " " << s.y() << std::endl;
        std::cout << "found triangle: " << std::endl;
        std::cout << i0 << " [" << points[i0].point() << "] " << std::endl;
        std::cout << i1 << " [" << points[i1].point() << "] " << std::endl;
        std::cout << i2 << " [" << points[i2].point() << "] " << std::endl;
#endif

        // we're inside! compute the distance to the vertices of the triangle
        for(int j=0; j<3; ++j)
        {
          Grid_point& gp = points[fh->vertex(j)->info()];
          const Metric& seed_m = mf->compute_metric(s);
          const Metric& v_m = gp.metric();
          const Eigen::Matrix2d& f = get_interpolated_transformation(seed_m, v_m);

          Vector2d v;
          v(0) = s.x() - gp.point().x();
          v(1) = s.y() - gp.point().y();
          v = f * v;
          FT d = v.norm(); // d = std::sqrt(v^t * M * v) = (f*v).norm()
          initialize_grid_point(gp.index, d, seed_id);
        }
        break; // no need to check the other primitives
      }
    }
    CGAL_assertion(found);
  }

  void build_aabb_tree()
  {
    tree.clear();

    // create an aabb tree of triangles to fasten it
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

        const Point_2& p0 = points[fh->vertex(0)->info()].point();
        const Point_2& p1 = points[fh->vertex(1)->info()].point();
        const Point_2& p2 = points[fh->vertex(2)->info()].point();

        Point_3 p(p0.x(), p0.y(), 0.);
        Triangle_3 tr_3(p,
                        Point_3(p1.x(), p1.y(), 0.),
                        Point_3(p2.x(), p2.y(), 0.));
        Primitive pri(fh, tr_3, p);
        tree.insert(pri);
      }
    }

    std::cout << "built aabb tree w/ " << tree.size() << " entries" << std::endl;
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

  void dual_shenanigans(Grid_point& gp)
  {
    CGAL_assertion(gp.state == KNOWN);
    const Star_handle star = ss[gp.index];

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
        const Grid_point& gq = points[fh->vertex(j)->info()];
        if(gq.state != KNOWN)
          continue;

        std::size_t q_id = gq.closest_seed_id;
        CGAL_assertion(q_id != static_cast<std::size_t>(-1));
        dual_simplex.insert(q_id);
      }
      add_simplex_to_triangulation(gp, dual_simplex);
    }
  }

  bool get_next_trial_point(std::size_t& id)
  {
    while(!trial_points.empty())
    {
      id = trial_points.front();
      std::pop_heap(trial_points.begin(), trial_points.end(),
                    Grid_point_comparer<Self>(this));
      trial_points.pop_back();

      CGAL_assertion(id < points.size());

      const Grid_point& gp = points[id];
      if(gp.state != TRIAL)
      {
        std::cout << "WARNING : point with a state non-TRIAL in the PQ : ";
        std::cout << gp.index << "... Ignoring it!" << std::endl;
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

#if (verbose > 40)
      std::cout << "trial heap: " << std::endl;
      for (std::vector<Grid_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif

      if(!get_next_trial_point(id))
        break;

      Grid_point& gp = points[id];

#if (verbose > 10)
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "picked n° " << gp.index << " (" << gp.point() << ") ";
      std::cout << "at distance : " << gp.distance_to_closest_seed << " from the seed "
                << gp.closest_seed_id << " ancestor : " << gp.ancestor << std::endl;
#endif

      gp.state = KNOWN;
      known_count++; // tmp --> use change_state()

      if(false && use_dual_shenanigans) // tmp
        dual_shenanigans(gp);

      PQ_state pqs = gp.update_neighbors_distances(trial_points);
      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(),
                       Grid_point_comparer<Self>(this));
      is_t_empty = trial_points.empty();
    }

    fix_distance_and_seed_id_at_corners();
    CGAL_expensive_assertion_code(debug());

    std::cout << "End of spread_distances. time: ";
    std::cout << ( std::clock() - s_start ) / (double) CLOCKS_PER_SEC << std::endl;
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
      const Grid_point& gp = points[i];

      if(gp.ancestor != static_cast<std::size_t>(-1) &&
         points[gp.ancestor].closest_seed_id != gp.closest_seed_id)
      {
         failed = true;
         std::cout << "ancestor/children don't agree on the seed @ " << gp.index << std::endl;
      }

      if(gp.ancestor != static_cast<std::size_t>(-1) &&
         points[gp.ancestor].children.find(gp.index) == points[gp.ancestor].children.end())
      {
        failed = true;
        std::cout << "failure in ancestor/children relationship at " << gp.index << std::endl;
        std::cout << "ancestor is : " << gp.ancestor << " color: " << gp.closest_seed_id << std::endl;
      }

      if(gp.distance_to_closest_seed == FT_inf || gp.closest_seed_id >= seeds.size() ||
         (gp.ancestor != static_cast<std::size_t>(-1) &&
          points[gp.ancestor].closest_seed_id != gp.closest_seed_id))
      {
        failed = true;
        std::cout << "debug: " << i << " [" << gp.point() << "] " << " at distance : ";
        std::cout << gp.distance_to_closest_seed << " from " << gp.closest_seed_id;
        std::cout << " ancestor: " << gp.ancestor;
        std::cout << " ancid: ";
        if(gp.ancestor==static_cast<std::size_t>(-1))
           std::cout << "-1" << std::endl;
        else
          std::cout << points[gp.ancestor].closest_seed_id << std::endl;
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
    std::cout << "grid reset" << std::endl;

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      Grid_point& gp = points[i];
      gp.reset();
    }

    known_count = 0;
    trial_count = 0;
    far_count = points.size();

    clear_dual();
  }

  void refresh_grid_point_states()
  {
    CGAL_assertion(trial_points.empty());
    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
      points[i].state = FAR;

    known_count = 0;
    trial_count = 0;
    far_count = points.size();
  }

  void refine_seeds(const Point_2 new_seed)
  {
    vertices_nv = insert_new_seed(new_seed.x(), new_seed.y());

#ifdef USE_FULL_REBUILDS
    reset();
    locate_and_initialize_seeds();
#else
    refresh_grid_point_states(); // we can't spread from the new seed if all states are 'KNOWN'
    locate_and_initialize(new_seed, seeds.size()-1);
#endif

    clear_dual();
    spread_distances(false/*use_dual_shenanigans*/);
    compute_dual();
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
      std::cout << "Couldn't find a ref point, need more initial points (or we are done)" << std::endl;
      return false;
    }

    Point_2 ref_point = points[best_entry.get<1>()].point();

    std::cout << "naturally we picked : " << best_entry.get<1>() << "[ ";
    std::cout << ref_point << " ]" << std::endl;
    std::cout << "second was: " << best_entry.get<2>() << std::endl;

    refine_seeds(ref_point);
    return true;
  }

  void output_grid(const std::string str_base) const
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
      const Grid_point& gp = points[i];
      out << gp.point() << " " << i+1 << '\n';

      out_bb << gp.distance_to_closest_seed << '\n';
//      out_bb << gp.closest_seed_id << '\n';
//      out_bb << (gp.is_Voronoi_vertex?"1":"0") << '\n';
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

      const Grid_point& git = points[it->get<1>()];
      const Grid_point& ge = points[entry.get<1>()];

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

  void test_simplex(const Grid_point& gp,
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
      FT gpd = gp.distance_to_closest_seed;
      if(gpd > max_size)
      {
        PQ_entry pqe = boost::make_tuple(dual_simplex, gp.index, gpd);
        insert_in_PQ(pqe, size_queue, 0);
        size_queue.insert(pqe);
        return;
      }
    }

    // distortion
    if(max_distortion > 1.)
    {
      FT gamma = gp.distortion_to_seed();
      if(gamma > max_distortion)
      {
        PQ_entry pqe = boost::make_tuple(dual_simplex, gp.index, gamma);
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
          PQ_entry pqe = boost::make_tuple(dual_simplex, gp.index, area);
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
          PQ_entry pqe = boost::make_tuple(dual_simplex, gp.index, 1./qual);
          insert_in_PQ(pqe, quality_queue, 3);
          return;
        }
      }
    }
  }

  void test_simplex(const Tri& grid_tri,
                    const Simplex& dual_simplex)
  {
    Tri::const_iterator it = grid_tri.begin();
    Tri::const_iterator end = grid_tri.end();
    for(; it!=end; ++it)
     test_simplex(points[*it], dual_simplex);
  }

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

  void add_simplex_to_triangulation(const Grid_point& gp,
                                    const Simplex& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    std::size_t length = gp.ancestor_path_length;
    if(length < min_ancestor_path_length)
    {
//      std::cout << "the canvas is thin (" << length << ") at : " << gp.index << " ("
//                << gp.point() << ") dual simplex of size: " << dual_simplex.size()
//                << " and cellid: " << gp.closest_seed_id << std::endl;
    }

    // test against criteria should be moved somewhere else
    // than during the dual computations todo
    if(dual_simplex.size() == 3 ||
       (dual_simplex.size() == 2 && gp.is_on_domain_border))
      test_simplex(gp, dual_simplex);

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
      const Grid_point& gp = points[fh->vertex(i)->info()];
      const FT d = gp.distance_to_closest_seed;

      if(d > max)
      {
        max = d;
        max_dist_p = fh->vertex(i)->info();
      }
    }

    // keep as dual point the farthest grid point -- at least until we have
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

  void output_grid_data_and_dual(const std::string str_base)
  {
    output_grid(str_base);
    output_straight_dual(str_base);
  }

  Base_mesh(Star_set& ss_)
    :
      ss(ss_),
      points(),
      trial_points(),
      dual_edges(),
      dual_triangles(),
      size_queue(this),
      distortion_queue(this),
      intersection_queue(this),
      quality_queue(this),
      Voronoi_vertices()
  {
    initialize_from_starset();
  }
};

const Point_2& Grid_point::point() const
{
  return bm->ss[index]->center_point();
}

const Metric& Grid_point::metric() const
{
  return bm->ss[index]->metric();
}

void Grid_point::change_state(FMM_state new_state)
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

std::size_t Grid_point::anc_path_length() const
{
  std::size_t i = 0;
  std::size_t n_anc = ancestor;
  while(n_anc != static_cast<std::size_t>(-1))
  {
    ++i;
    const Grid_point& anc = bm->points[n_anc];
    n_anc = anc.ancestor;

    CGAL_assertion(i < bm->ss.size());
  }
  return i;
}

void Grid_point::print_children() const
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

void Grid_point::print_ancestors() const
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

bool Grid_point::find_in_ancestry(const std::size_t i) const
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

void Grid_point::remove_from_children(const std::size_t c)
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

void Grid_point::mark_descendants(std::size_t& count)
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
    Grid_point& cp = bm->points[*cit];
    if(cp.state == ORPHAN)
      continue;
    cp.mark_descendants(count);
  }
}

void Grid_point::reset_descendants()
{
//  std::cout << "reset descendant at: " << index << std::endl;
  CGAL_precondition(state == ORPHAN);

  while(!children.empty())
  {
    Grid_point& cp = bm->points[*(children.begin())];
    cp.reset_descendants();
  }

  // we might reset points that are in the priority queue here !
  // but it's okay since get_next_trial_point() will ignore these
  reset();
  state = ORPHAN;
}

void Grid_point::deal_with_descendants()
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

bool Grid_point::compute_closest_seed(const std::size_t n_anc,
                                      const bool overwrite)
{
#if (verbose > 20)
  std::cout << "closest seed at " << index << " from " << n_anc;
  std::cout << " (current d/a: " << distance_to_closest_seed << " " << ancestor;
  std::cout << " seed: " << closest_seed_id << ")" << std::endl;
#endif

  CGAL_precondition(static_cast<std::size_t>(n_anc) < bm->points.size());
  const Grid_point& anc = bm->points[n_anc];
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
    const Grid_point& curr_ancestor = bm->points[n_curr_ancestor];

    Vector2d ancestor_edge;
    ancestor_edge(0) = point().x() - curr_ancestor.point().x();
    ancestor_edge(1) = point().y() - curr_ancestor.point().y();
    FT ancestor_edge_length = ancestor_edge.norm();
    Vector2d normalized_anc_edge = ancestor_edge/ancestor_edge_length;

    // compute the distance for the current depth (i)
    FT dist_to_ancestor = 0.;
    for(int j=0; j<i; ++j) // we add a part for each edge in the path
    {
      // get the metric for the current edge

      CGAL_assertion(ancestor_path[j] < bm->points.size() &&
                     ancestor_path[j+1] < bm->points.size());
      const Grid_point& e0 = (j==0)?*this:bm->points[ancestor_path[j]];
      const Grid_point& e1 = bm->points[ancestor_path[j+1]];

      const Metric& m0 = e0.metric();
      const Metric& m1 = e1.metric();

      Vector2d curr_edge;
      curr_edge(0) = e0.point().x() - e1.point().x();
      curr_edge(1) = e0.point().y() - e1.point().y();

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

      while(!children.empty())
      {
        Grid_point& gp = bm->points[*(children.begin())];
        gp.deal_with_descendants();
        CGAL_postcondition(gp.anc_path_length() == gp.ancestor_path_length); // checks for circular ancestry
      }
#endif
    }
    return true;
  }
  return false;
}

void Grid_point::update_neighbor_distance(Grid_point& gp,
                                          std::vector<std::size_t>& trial_pq,
                                          PQ_state& pqs_ret) const
{
#if (verbose > 12)
  std::cout << "neighbor: " << gp.index << " has state: " << gp.state << std::endl;
#endif

  if(gp.state == KNOWN)
    return;
  else if(gp.state == TRIAL)
  {
    // note that we don't insert in trial_pq since it's already in
    if(gp.compute_closest_seed(index))
      pqs_ret = REBUILD_TRIAL;
  }
  else
  {
    CGAL_assertion(gp.state == FAR || gp.state == ORPHAN);
    if(gp.compute_closest_seed(index))
    {
      gp.state = TRIAL;
      trial_pq.push_back(gp.index);
      std::push_heap(trial_pq.begin(), trial_pq.end(),
                     Grid_point_comparer<Base_mesh>(bm));
    }
  }
}

PQ_state Grid_point::update_neighbors_distances(std::vector<std::size_t>& trial_pq) const
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

    Grid_point& gp = bm->points[vh->info()];
    update_neighbor_distance(gp, trial_pq, pqs_ret);
  }

#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
  Point_list::const_iterator lit = reciprocal_neighbors.begin();
  Point_list::const_iterator lend = reciprocal_neighbors.end();
  for(; lit!=lend; ++lit)
  {
    if(is_a_corner_vertex(*lit))
      continue;

    Grid_point& gp = bm->points[*lit];
    update_neighbor_distance(gp, trial_pq, pqs_ret);
  }
#endif

  return pqs_ret;
}

FT Grid_point::distortion_to_seed() const
{
  FT gamma = 1.;
  //    std::cout << "init gamma: " << gamma << " " << index << std::endl;
  const Grid_point* curr = this;
  std::size_t n_anc = ancestor;

  while(n_anc != static_cast<std::size_t>(-1))
  {
    const Grid_point& anc = bm->points[n_anc];
    const Metric& m1 = anc.metric();
    const Metric& m2 = curr->metric();
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

void Grid_point::reset()
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
  distance_to_closest_seed = FT_inf;
  closest_seed_id = -1;
  ancestor = -1;
  children.clear();
  is_Voronoi_vertex = false;
  ancestor_path_length = 0;
}

void Grid_point::initialize_from_point(const FT d,
                                       const std::size_t seed_id)
{
  reset();

  closest_seed_id = seed_id;
  distance_to_closest_seed = d;
  state = TRIAL;
  ancestor = -1;
}
bool Grid_point::operator==(const Grid_point& gp) const
{
  // the only thing that really matters is the point
  return (point() == gp.point());
}

Grid_point::Grid_point()
  :
    bm(NULL),
    index(static_cast<std::size_t>(-1)),
    is_on_domain_border(false),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(),
#endif
    state(FAR),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    children(),
    is_Voronoi_vertex(false),
    ancestor_path_length(0)
{ }

Grid_point::Grid_point(Base_mesh * bm_,
                       const std::size_t index_,
                       const bool is_on_domain_border_)
  :
    bm(bm_),
    index(index_),
    is_on_domain_border(is_on_domain_border_),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(),
#endif
    state(FAR),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    children(),
    is_Voronoi_vertex(false),
    ancestor_path_length(0)
{ }

Grid_point::Grid_point(const Grid_point& gp)
  :
    bm(gp.bm),
    index(gp.index),
    is_on_domain_border(gp.is_on_domain_border),
#ifdef ANISO_USE_RECIPROCAL_NEIGHBORS
    reciprocal_neighbors(gp.reciprocal_neighbors),
#endif
    state(gp.state),
    distance_to_closest_seed(gp.distance_to_closest_seed),
    closest_seed_id(gp.closest_seed_id),
    ancestor(gp.ancestor),
    children(gp.children),
    is_Voronoi_vertex(gp.is_Voronoi_vertex),
    ancestor_path_length(gp.ancestor_path_length)
{ }

void initialize_seeds()
{
  vertices_nv = build_seeds();
  CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );
}

struct Point_extracter
    : public std::unary_function<Grid_point, Point_2>
{
  const Point_2 operator()(const Grid_point& gp) const
  {
    return gp.point();
  }
};

void draw_metric_field(const MF mf, const Base_mesh& bm)
{
  typedef boost::transform_iterator<Point_extracter,
      std::vector<Grid_point>::const_iterator> Extracted_iterator;
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

//  draw_metric_field(mf, bm);
//  exit(0);

  initialize_seeds();
  bm.locate_and_initialize_seeds();
  bm.spread_distances(true/*use_dual_shenanigans*/);

  if(n_refine > 0)
  {
    bm.output_grid_data_and_dual("pre_ref");

    for(int i=0; i<n_refine; ++i)
    {
      // can't compute the dual while spreading if we're refining since we have already a layer
      // of paint laying on the canvas...
      bool successful_insert = bm.refine_seeds_with_self_computed_ref_point();

      if(i%100 == 0)
      {
        std::ostringstream out;
        out << "ref_" << seeds.size();
        bm.output_grid_data_and_dual(out.str());
      }

      if(!successful_insert)
        break;
    }

    for(std::size_t i=0, ps=bm.points.size(); i<ps; ++i)
      bm.points[i].state = KNOWN;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "End refinement: " << duration << std::endl;
  }

  bm.output_grid_data_and_dual("starset_tr");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "duration: " << duration << std::endl;
}
