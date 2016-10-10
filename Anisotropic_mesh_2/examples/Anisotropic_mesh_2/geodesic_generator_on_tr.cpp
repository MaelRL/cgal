// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

#include <CGAL/Starset.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/External_metric_field.h>

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Mesh_2/Sizing_field_2.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Unique_hash_map.h>

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

#include <deque>
#include <iostream>
#include <functional>
#include <limits>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::FT                                                FT;

using namespace CGAL::Anisotropic_mesh_2;

typedef K::Vector_2                                         Vector;
typedef Metric_base<K>                                      Metric;

typedef K::Point_2                                          Point_2;
typedef K::Vector_2                                         Vector_2;
typedef K::Point_3                                          Point_3;
typedef std::set<std::size_t>                               Simplex;
typedef boost::array<std::size_t, 2>                        Edge;
typedef boost::array<std::size_t, 3>                        Tri;
typedef Eigen::Matrix<double, 2, 1>                         Vector2d;

typedef K::Segment_2                                        Segment;
typedef K::Triangle_2                                       Triangle_2;
typedef K::Segment_3                                        Segment_3;
typedef K::Triangle_3                                       Triangle_3;

typedef CGAL::Exact_predicates_exact_constructions_kernel   KExact;
typedef KExact::Point_2                                     EPoint;
typedef KExact::Segment_2                                   ESegment;
typedef KExact::Triangle_2                                  ETriangle;

typedef CGAL::Cartesian_converter<K, KExact>                To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                Back_from_exact;

// 'bit ugly to have the metric field running around here but we need to pass
// it around and around and around and it's faster this way !
//typedef CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
typedef CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;
//typedef CGAL::Anisotropic_mesh_2::External_metric_field<K>* MF;
MF mf;

// stuff to generate the grid using Mesh_2 since Aniso_mesh_2 is a turtle.
// Need to define our own refinement criterion based on metric shenanigans
bool ignore_children = false;

To_exact to_exact;
Back_from_exact back_from_exact;

// #define COMPUTE_PRECISE_VOR_VERTICES
#define COMPUTE_DUAL_FOR_ALL_DIMENSIONS
//#define USE_FULL_REBUILDS
//#define REFINE_GRID
#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 10

const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 10;
static const int k = 8; // depth of the ancestor edge

// the metric field and the seeds
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

// seeds
const std::string seeds_str = "starred_";

// base mesh
const std::string str_base_mesh = "starred";
CGAL::Bbox_2 base_mesh_bbox;

// refinement
int n_refine = 0;
std::size_t min_ancestor_path_length = 5;

// optimization
int max_opti_n = 0;
int max_depth = -1; // how precise are the Voronoi vertices (bigger = better)

//debug & info
int known_count=0, trial_count=0, far_count=0;
std::clock_t start;
std::vector<int> optimal_depths(k, 0);

typedef boost::unordered_map<boost::array<std::size_t, k+1>, FT> Path_map;
Path_map computed_paths;

namespace CGAL
{

template<class CDT>
class Geo_sizing_field
    : public virtual CGAL::Sizing_field_2<CDT>
{
public:
  typedef typename CDT::Geom_traits             GT;
  typedef typename GT::FT                       FT;
  typedef typename GT::Point_2                  Point_2;

  FT operator()(const Point_2& p) const
  {
    FT base = 1.0;

    const Metric& m = mf->compute_metric(p);
    const FT width = 1./std::sqrt(m.get_max_eigenvalue());

//    std::cout << "mM: " << m.get_min_eigenvalue() << " " << m.get_max_eigenvalue() << std::endl;
//    std::cout << "width: " << width << std::endl;

//    return 0.05;

    CGAL_assertion(m.get_max_eigenvalue() >= m.get_min_eigenvalue());

    FT discretization = 6.;
    FT metric_based_size = width / discretization;

//    std::cout << "sizing field: " << m.get_max_eigenvalue() << " base: " << base
//              << " mbase: " << metric_based_size << std::endl;

    return (std::min)(base, metric_based_size);
  }

  Geo_sizing_field() { }
};

template<class CDT>
class Delaunay_mesh_distortion_criteria_2 :
    public virtual CGAL::Delaunay_mesh_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits                          Geom_traits;
  typedef Geo_sizing_field<CDT>                              Sizing_field;

  FT distortionbound;
  Sizing_field sizing_field;

public:
  typedef CGAL::Delaunay_mesh_criteria_2<CDT>           Base;

  FT distortion_bound() const { return distortionbound; }
  void set_distortion_bound(const double db) { distortionbound = db; }

  struct Quality : public boost::tuple<FT, FT, FT>
  {
    typedef boost::tuple<FT, FT, FT> Base;

    const FT& distortion() const { return get<0>(); }
    const FT& size() const { return get<1>(); }
    const FT& sine() const { return get<2>(); }

    // q1 == *this, q2 == q
    // q1 < q2 means q1 is prioritised over q2
    bool operator<(const Quality& q) const
    {
      if( distortion() > 1)
        if( q.distortion() > 1 )
          return ( distortion() > q.distortion() );
        else
          return true; // *this is distorted but not q
      else
        if( q.distortion() > 1 )
          return false; // q is distorted but not *this

      if( size() > 1 )
        if( q.size() > 1 )
          return ( size() > q.size() );
        else
          return true; // *this is big but not q
      else
        if( q.size() >  1 )
          return false; // q is big but not *this

      return( sine() < q.sine() );
    }

    std::ostream& operator<<(std::ostream& out) const
    {
      return out << "(distortion=" << distortion() << ", size=" << size()
                 << ", sine=" << sine() << ")";
    }

    Quality() : Base() {}
    Quality(FT _distortion, FT _size, FT _sine)
      :
        Base(_distortion, _size, _sine)
    { }
  };

  class Is_bad : public Base::Is_bad
  {
  protected:
    FT distortion_bound;
    const Sizing_field* sizing_field;

  public:
    typedef typename Base::Is_bad::Point_2                   Point_2;

    CGAL::Mesh_2::Face_badness operator()(const Quality q) const
    {
      if( q.distortion() > 1. )
        return CGAL::Mesh_2::IMPERATIVELY_BAD;
      if( q.size() > 1. )
        return CGAL::Mesh_2::IMPERATIVELY_BAD;
      if( q.sine() < this->B )
        return CGAL::Mesh_2::BAD;
      else
        return CGAL::Mesh_2::NOT_BAD;
    }

    CGAL::Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
                                          Quality& q) const
    {
      typedef typename CDT::Geom_traits Geom_traits;

      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      q.get<0>() = 1;
      if( distortion_bound > 1)
      {
        const Metric& ma = mf->compute_metric(pa);
        const Metric& mb = mf->compute_metric(pb);
        const Metric& mc = mf->compute_metric(pc);

        FT gamma_ab = ma.compute_distortion(mb);
        FT gamma_ac = ma.compute_distortion(mc);
        FT gamma_bc = mb.compute_distortion(mc);

        FT max_gamma = (std::max)((std::max)(gamma_ab, gamma_ac), gamma_bc);
        q.get<0>() = max_gamma / distortion_bound;

//        std::cout << max_gamma << std::endl;

        if( q.distortion() > 1 )
        {
          q.get<1>() = 1.; // no need to compute the size
          q.get<2>() = 1.; // no need to compute the sine
          return CGAL::Mesh_2::IMPERATIVELY_BAD;
        }
      }

      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2
          Compute_squared_distance_2;

      Geom_traits traits;

      Compute_squared_distance_2 squared_distance =
          traits.compute_squared_distance_2_object();

      double a = CGAL::to_double(squared_distance(pb, pc));
      double b = CGAL::to_double(squared_distance(pc, pa));
      double c = CGAL::to_double(squared_distance(pa, pb));

      double max_sq_length; // squared max edge length
      double second_max_sq_length;

      if(a<b)
      {
        if(b<c) {
          max_sq_length = c;
          second_max_sq_length = b;
        }
        else { // c<=b
          max_sq_length = b;
          second_max_sq_length = ( a < c ? c : a );
        }
      }
      else // b<=a
      {
        if(a<c) {
          max_sq_length = c;
          second_max_sq_length = a;
        }
        else { // c<=a
          max_sq_length = a;
          second_max_sq_length = ( b < c ? c : b );
        }
      }

      q.get<1>() = 0;
      const Point_2& g = centroid(pa, pb, pc);
      FT squared_size_bound = sizing_field->operator()(g);
      squared_size_bound *= squared_size_bound;
      if( squared_size_bound != 0 )
      {
        q.get<1>() = max_sq_length / squared_size_bound;
        // normalized by size bound to deal with size field
        if( q.size() > 1 )
        {
          q.get<2>() = 1; // no need to compute the sine
          return CGAL::Mesh_2::IMPERATIVELY_BAD;
        }
      }

      Compute_area_2 area_2 = traits.compute_area_2_object();
      double area = 2*CGAL::to_double(area_2(pa, pb, pc));
      q.get<2>() = (area * area) / (max_sq_length * second_max_sq_length); // sine

      if( q.sine() < this->B )
        return CGAL::Mesh_2::BAD;
      else
        return CGAL::Mesh_2::NOT_BAD;
    }

    Is_bad(const FT aspect_bound,
           const Sizing_field* sizing_field_,
           const FT dist_bound,
           const Geom_traits& traits)
      :
        Base::Is_bad(aspect_bound, traits),
        distortion_bound(dist_bound),
        sizing_field(sizing_field_)
    { }
  };

  Is_bad is_bad_object() const
  {
    return Is_bad(this->bound(), &sizing_field, this->distortion_bound(),
                  this->traits /* from the bad class */);
  }

  Delaunay_mesh_distortion_criteria_2(const FT aspect_bound = 0.125,
                                      const FT distortion_bound = 0.,
                                      const Geom_traits& traits = Geom_traits())
    :
      Base(aspect_bound, traits),
      distortionbound(distortion_bound),
      sizing_field()
  { }
};
}

template<typename CDT>
bool is_vertex_on_border(const CDT& cdt,
                         typename CDT::Vertex_handle& v)
{
  typename CDT::Vertex_circulator vc = cdt.incident_vertices(v), done(vc);
  if(!vc.is_empty())
  {
    do
    {
      if(cdt.is_infinite(vc))
        return true;
    }
    while(++vc != done);
  }
  return false;
}

template<typename CDT>
void output_cdt_to_mesh(const CDT& cdt,
                        const std::string str_base)
{
  CGAL_assertion(cdt.dimension() == 2);

  std::size_t n = cdt.number_of_vertices();
  std::size_t m = cdt.number_of_faces();

  CGAL::Unique_hash_map<typename CDT::Vertex_handle, std::size_t> vertex_map;

  std::ofstream out((str_base + ".mesh").c_str());
  out << std::setprecision(17);
  out << "MeshVersionFormatted 1" << '\n';
  out << "Dimension 2" << '\n';
  out << "Vertices" << '\n';
  out << n << '\n';

  std::size_t v_counter = 1;
  for(typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                                             vit != cdt.finite_vertices_end(); ++vit)
  {
    typename CDT::Vertex_handle vh = vit;
    vertex_map[vh] = v_counter++;
    out << *vit << " " << is_vertex_on_border<CDT>(cdt, vh) << '\n';
  }

  out << "Triangles" << '\n';
  out << m << '\n';

  for(typename CDT::Face_iterator fit = cdt.finite_faces_begin();
                                  fit != cdt.finite_faces_end(); ++fit)
  {
    for(int i=0; i<=cdt.dimension(); ++i)
      out << vertex_map[fit->vertex(i)] << " ";
    out << "0" << '\n';
  }
  out << "End" << '\n';
}

void generate_grid()
{
  typedef CGAL::Triangulation_vertex_base_2<K>                    Vb;
  typedef CGAL::Delaunay_mesh_face_base_2<K>                      Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>            Tds;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>      CDT;
  typedef CGAL::Delaunay_mesh_distortion_criteria_2<CDT>          Criteria;

  typedef CDT::Vertex_handle                                      Vertex_handle;
  typedef CDT::Point                                              Point;

  CDT cdt;

  FT a = 5;
  Vertex_handle va = cdt.insert(Point(-a, -a));
  Vertex_handle vb = cdt.insert(Point(a, -a));
  Vertex_handle vc = cdt.insert(Point(a, a));
  Vertex_handle vd = cdt.insert(Point(-a, a));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  FT shape = 0.125;
  FT distortion = 1.; // > 1 to use
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(shape, distortion));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  output_cdt_to_mesh(cdt, str_base_mesh + "_base");
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

  // immuable stuff
  Base_mesh* bm;
  Point_2 point;
  std::size_t index;
  Metric metric;
  Point_set neighbors;
  std::list<std::size_t> incident_triangles;
  bool is_on_domain_border;

  // stuff that depends on the seeds
  FMM_state state;
  int is_seed_holder; // fixme this is ugly, a bimap [seeds, grid_points_id] would be best
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
  void remove_from_neighbors(const std::size_t n_p);
  void remove_from_incident_triangles(const std::size_t i);
  void mark_descendants(std::size_t& count, boost::unordered_set<std::size_t>& potential_parents);
  void reset_descendants();
  void deal_with_descendants();
  bool compute_closest_seed(const std::size_t anc_id, const bool overwrite = true);
  PQ_state update_neighbors_distances(std::vector<Grid_point*>& trial_pq) const;
  FT distortion_to_seed() const;
  void reset();
  void initialize_from_point(const FT d, const std::size_t seed_id);

  bool operator==(const Grid_point& gp) const;

  Grid_point();
  Grid_point(Base_mesh* bm_,
             const Point_2& point_,
             const std::size_t index_,
             const bool is_on_domain_border_ = false);
  Grid_point(const Grid_point& gp);
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
  std::ifstream in((seeds_str + ".mesh").c_str());
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

    if(seeds.size() >= vertices_nv)
      break;
  }
#if (verbose > 0)
  std::cout << "seeds: " << seeds.size() << std::endl;
#endif
  return seeds.size();
}

inline std::ostream& operator<<(std::ostream& os,
                                const boost::tuple<Simplex, const Grid_point*, FT>& pqe)
{
  std::cout << "PQ entry element:" << std::endl;
  std::cout << "Simplex: " << pqe.get<0>();
  std::cout << "gp: " << pqe.get<1>()->index << " [" << pqe.get<1>()->point << "] ";
  std::cout << pqe.get<1>()->distance_to_closest_seed << std::endl;
  std::cout << "val: " << pqe.get<2>() << std::endl;

  return os;
}

template<typename PQ_entry>
struct PQ_entry_comparator
{
  bool operator()(const PQ_entry& left, const PQ_entry& right) const
  {
    if(left.get<2>() == right.get<2>()) // that's the value
    {
      if(left.get<1>()->distance_to_closest_seed ==
         right.get<1>()->distance_to_closest_seed) // that's the dist to the ref point
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
};

struct Base_mesh
{
  typedef std::vector<Grid_point*>                           Grid_point_vector;

  typedef std::vector<Grid_point>                            Voronoi_vertices_container;
  typedef std::vector<Voronoi_vertices_container>            Voronoi_vertices_vector;

  // Refinement
  typedef boost::tuple<Simplex, const Grid_point*, FT>       PQ_entry;
  typedef std::set<PQ_entry, PQ_entry_comparator<PQ_entry> > PQ;

  // For the approximate geodesic triangulation...
  // the first pair are the IDs of two Voronoi cells (which a non empty inter)
  // the second pair are the two closest points on the bisector (one in each cell)
  // WARNING: currently it's on the part of the bisector that exists, not the whole
  // bisector, which means we have a 'local' geodesic triangulation, not the real one...
  typedef std::map<std::pair<std::size_t, std::size_t>,
                   std::pair<const Grid_point*, const Grid_point*> > Edge_bisec_map;

  // for fast seed locate (cheating by pretending we're in 3D)

  struct Base_mesh_primitive
  {
    typedef int                            Id;

    typedef K::Point_3            Point;
    typedef K::Triangle_3         Datum;

    Id m_it;
    Datum m_datum; // cache the datum
    Point m_p; // cache the reference point

    Base_mesh_primitive() { }
    Base_mesh_primitive(int i, const Datum& datum, const Point& point)
      :
        m_it(i),
        m_datum(datum),
        m_p(point)
    { }

     Id id() const { return m_it; }
     Datum datum() const { return m_datum; }
     Point reference_point() const { return m_p; }
  };

  typedef Base_mesh_primitive                             Primitive;
  typedef Primitive::Id                          Primitive_Id;
  typedef CGAL::AABB_traits<K, Primitive>                 AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits>           Tree;

  std::vector<Grid_point> points;
  std::vector<Tri> triangles;
  std::vector<Tri> triangles_buffer; // used while refining the grid
  Grid_point_vector trial_points;

  // extra info: map[edge e] = triangles incident to e
  std::map<std::pair<std::size_t, std::size_t>, std::list<std::size_t> > edge_incident_tris;

  // Duality
  mutable Edge_bisec_map edge_bisectors;
  mutable boost::unordered_set<Edge> dual_edges;
  mutable boost::unordered_set<Tri> dual_triangles;

  // Refinement
  PQ size_queue;
  PQ distortion_queue;
  PQ intersection_queue;
  PQ quality_queue;

  mutable const Grid_point* refinement_point;

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

  void initialize_grid_point(Grid_point* gp,
                             const FT dist,
                             const std::size_t seed_id)
  {
    CGAL_precondition(gp);
    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cout << "WARNING: a new seed is overwriting the closest seed id at ";
      std::cout << gp->index << ". Previous index is: " << gp->closest_seed_id;
      std::cout << " (" << seeds.size() << " seeds)" << std::endl;
    }

    if(gp->state == TRIAL)
    {
      // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");
    }

    std::cout << "initialize grid point: " << gp->index << " (" << gp->point << ")" << std::endl;
    std::cout << "at distance: " << dist << " from " << seed_id << std::endl;

    gp->initialize_from_point(dist, seed_id);
    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(),
                   Grid_point_comparer<Grid_point>());
  }

  void locate_and_initialize(const Point_2& s,
                             const std::size_t seed_id)
  {
    CGAL_precondition(!tree.empty());

    // find the triangle that contains the seed and initialize the distance at
    // these three points accordingly.
    bool found = false;

    Segment_3 query(Point_3(s.x(), s.y(), -1.), Point_3(s.x(), s.y(), 1.));
    std::list<Primitive_Id> intersections;
    tree.all_intersected_primitives(query, std::back_inserter(intersections));

    std::list<Primitive_Id>::const_iterator it = intersections.begin(),
                                            end = intersections.end();
    for(; it!=end; ++it)
    {
      const Tri& tr = triangles[*it];
      const Triangle_2 tr_2(points[tr[0]].point,
                            points[tr[1]].point,
                            points[tr[2]].point);

      if(K().has_on_bounded_side_2_object()(tr_2, s) ||
         K().has_on_boundary_2_object()(tr_2, s))
      {
        found = true;
#if (verbose > 10)
        std::cout << "locating seed " << seed_id
                  << " point: " << s.x() << " " << s.y() << std::endl;
        std::cout << "found triangle: " << *it << std::endl;
        std::cout << tr[0] << " [" << points[tr[0]].point << "] " << std::endl;
        std::cout << tr[1] << " [" << points[tr[1]].point << "] " << std::endl;
        std::cout << tr[2] << " [" << points[tr[2]].point << "] " << std::endl;
#endif

        // we're inside! compute the distance to the vertices of the triangle
        // only initialize the closest one
        FT min_d = FT_inf;
        std::size_t best_id = -1;
        for(int j=0; j<3; ++j)
        {
          Grid_point& gp = points[tr[j]];
          const Metric& seed_m = mf->compute_metric(s);
          const Metric& v_m = gp.metric;
          const Eigen::Matrix2d& f = get_interpolated_transformation(seed_m, v_m);

          Vector2d v;
          v(0) = s.x() - gp.point.x();
          v(1) = s.y() - gp.point.y();
          v = f * v;
          FT d = v.norm(); // d = std::sqrt(v^t * M * v) = (f*v).norm()
          if(d < min_d)
          {
            best_id = tr[j];
            min_d = d;
          }
        }

        initialize_grid_point(&(points[best_id]), min_d, seed_id);
        break; // no need to check the other primitives
      }
    }
    CGAL_postcondition(found);
  }

  void build_aabb_tree()
  {
    std::cout << "build aabb tree" << std::endl;
    tree.clear();

    // create an aabb tree of triangles to fasten it
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& t = triangles[i];
      const Point_2& p0 = points[t[0]].point;
      const Point_2& p1 = points[t[1]].point;
      const Point_2& p2 = points[t[2]].point;
      Point_3 p(p0.x(), p0.y(), 0.);
      Triangle_3 tr_3(p,
                      Point_3(p1.x(), p1.y(), 0.),
                      Point_3(p2.x(), p2.y(), 0.));
      Primitive pri(i, tr_3, p);
      tree.insert(pri);
    }
  }

  void initialize_base_mesh()
  {
#if (verbose > 5)
    std::cout << "base mesh initialization" << std::endl;
#endif

    // read and create the base mesh points
    std::ifstream in((str_base_mesh + "_base.mesh").c_str());

    if(!in)
      std::cout << "couldn't open base mesh" << std::endl;

    std::string word;
    std::size_t useless, nv, nt, dim;
    FT r_x, r_y;

    in >> word >> useless; // MeshVersionFormatted i
    in >> word >> dim; // Dimension d
    CGAL_assertion(dim == 2);

    in >> word >> nv; // Vertices nv
    for(std::size_t i=0; i<nv; ++i)
    {
      int border_info;
      in >> r_x >> r_y >> border_info;

      if(border_info != 0 && border_info != 1)
        CGAL_assertion(false && "you're not using a mesh with border info");

      Point_2 p(r_x, r_y);
      Grid_point gp(this, p, i, border_info);
      points.push_back(gp);
    }

    base_mesh_bbox = CGAL::Bbox_2();
    for(std::size_t i=0; i<nv; ++i)
      base_mesh_bbox += points[i].point.bbox();

    // read the triangles and assign the neighbors
    in >> word >> nt; // Triangles nt
    for(std::size_t i=0; i<nt; ++i)
    {
      Tri tr;
      in >> tr[0] >> tr[1] >> tr[2] >> useless;
      --tr[0]; --tr[1]; --tr[2]; // because we're reading medit data

      triangles.push_back(tr);
      std::size_t new_tri_id = triangles.size() - 1; // id of the current triangle

      // initialize the neighbors and border neighbors
      for(std::size_t i=0; i<3; ++i)
      {
        std::size_t id1 = tr[i];
        for(std::size_t j=1; j<=2; ++j)
        {
          std::size_t id2 = tr[(i+j)%3];
          points[id1].neighbors.insert(id2);
        }
        points[id1].incident_triangles.push_back(new_tri_id);
      }
    }

    build_aabb_tree();

#if (verbose > 5)
    std::cout << "base mesh initialized" << std::endl;
    std::cout << nv << " vertices" << std::endl;
    std::cout << nt << " triangles" << std::endl;
#endif
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

  void fill_edge_bisectors_map(const Grid_point* gp,
                               const Grid_point* gq)
  {
    // used to build pseudo geodesic triangulations...
    std::size_t p_id = gp->closest_seed_id;
    std::size_t q_id = gq->closest_seed_id;
    if(p_id != q_id) // we're on the bisector of [PQ]
    {
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

  void shave_off_virtual_points(const std::size_t real_points_n)
  {
//    std::cout << "shaving: " << points.size() << " to " << real_points_n << std::endl;
    std::vector<Grid_point>::iterator it = points.begin();
    std::advance(it, real_points_n);
    points.erase(it, points.end());
  }

  void mark_voronoi_vertices(const Tri& tri,
                             const Simplex& simplex)
  {
    std::size_t ds_size = simplex.size();

    // mark them as Voronoi vertices if the triangle has 3 different colors
    if(ds_size == 3)
    {
      for(std::size_t j=0; j<3; ++j)
      {
        Grid_point& gq = points[tri[j]];
        gq.is_Voronoi_vertex = true;
//        std::cout << gq.index << " (" << gq.point << ") vor vert for "
//                  << gq.closest_seed_id << std::endl;
      }
    }
  }

  void visit_neighbors_and_mark(Grid_point* gp,
                                std::vector<bool>& visited_points)
  {
    visited_points[gp->index] = true;
    std::size_t seed_p = gp->closest_seed_id;

    if(gp->index < 4) // corner hack fixme
      gp->is_Voronoi_vertex = true;

    Grid_point::Point_set::iterator it = gp->neighbors.begin();
    Grid_point::Point_set::iterator end = gp->neighbors.end();
    for(; it!=end; ++it)
    {
      Grid_point* gq = &(points[*it]);

      // below is kinda bad: imagine the hypothenuse of a corner triangle, not all
      // of these are border edges... fixme
      if(!gq->is_on_domain_border)
        continue;

      if(!visited_points[gq->index])
      {
        // check if they have different colors
        std::size_t seed_q = gq->closest_seed_id;
        if(seed_p != seed_q)
        {
          gp->is_Voronoi_vertex = true;
          gq->is_Voronoi_vertex = true;
        }

        // continue visiting from gq
        visit_neighbors_and_mark(gq, visited_points);
      }
    }
  }

  void mark_voronoi_vertices_on_border()
  {
    // Loop through the border of the domain; if a segment has two different colors
    // then we're on the intersection of a 1-bisector with the border and we add
    // each extremity of the segment as a voronoi vertex

    // Also, the corners have to be added as Voronoi vertices too

    // find a first point on the border. (points[0] is usually one)
    std::vector<bool> visited_points(points.size(), false);
    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      if(points[i].is_on_domain_border) // todo keep the border neighbor in memory
      {
        visit_neighbors_and_mark(&(points[i]), visited_points);
        break;
      }
    }
  }

  Grid_point Vor_vertex_in_triangle(const std::size_t n_q, const std::size_t n_r,
                                    const std::size_t n_s, int call_count = 0)
  {
    CGAL_precondition(n_q < points.size() && n_r < points.size() && n_s < points.size());

    // centroid is probably not the most optimal...
    Point_2 centroid = CGAL::centroid(points[n_q].point,
                                      points[n_r].point,
                                      points[n_s].point);

    std::size_t id = points.size();
    Grid_point c(this, centroid, id, false/*border info*/);
    // it's actually possible for the centroid to be on the domain border,
    // but it doesn't matter since it's not a real grid point
    // and we'll never use this border information

    // need to push it back to points anyway for compute_closest_seed to work
    points.push_back(c);

    Grid_point& gp = points.back();
    gp.compute_closest_seed(n_q);
    gp.compute_closest_seed(n_r);
    gp.compute_closest_seed(n_s);
    gp.state = KNOWN;

    // another potential stop is if the (max) distance between the centroid
    // and a point of the triangle is below a given distance
    if(call_count > max_depth)
      return gp;

    if(gp.closest_seed_id == points[n_q].closest_seed_id)
      return Vor_vertex_in_triangle(id, n_r, n_s, ++call_count);
    else if(gp.closest_seed_id == points[n_r].closest_seed_id)
      return Vor_vertex_in_triangle(id, n_q, n_s, ++call_count);
    else
    {
      CGAL_assertion(gp.closest_seed_id == points[n_s].closest_seed_id);
      return Vor_vertex_in_triangle(id, n_q, n_r, ++call_count);
    }
  }

  Grid_point compute_precise_Voronoi_vertex(const Tri& tri)
  {
    std::size_t real_points_n = points.size();
    Grid_point centroid = Vor_vertex_in_triangle(tri[0], tri[1], tri[2]);

    // clean off the virtual points created by Vor_vertex_in_triangle
    shave_off_virtual_points(real_points_n);

    return centroid;
  }

  void compute_precise_Voronoi_vertex(const Tri& tri,
                                      const Simplex& simplex)
  {
    // The function 'mark_Voronoi_vertices' only picks a point on the base mesh,
    // which is usually not the real Voronoi vertex of the Voronoi cell.
    // We need something more precise to compute the center of mass of the cell.

    // the idea is to virtually subdivide the triangle (we don't want to polute
    // the base mesh with more real triangles that probably will be useless when
    // another point is inserted and the bisectors move) with its center of mass

    // Then, compute the color at this new point, and subdivide the smaller triangle
    // that has 3 different colors. Repeat till happy.

    std::size_t ds_size = simplex.size();
    if(ds_size == 3)
    {
      const Grid_point& centroid = compute_precise_Voronoi_vertex(tri);

      for(int i=0; i<3; ++i)
      {
        const Grid_point& gp = points[tri[i]];
        std::size_t seed_id = gp.closest_seed_id;
        Voronoi_vertices[seed_id].push_back(centroid);
        Voronoi_vertices[seed_id].back().closest_seed_id = seed_id;
        Voronoi_vertices[seed_id].back().ancestor = gp.index;
      }
    }
  }

  Grid_point Vor_vertex_on_edge(const std::size_t n_q, const std::size_t n_r,
                                int call_count = 0)
  {
    // do the complicated split, see formula on notes // todo
    // taking the easy way for now

    Point_2 centroid = CGAL::barycenter(points[n_q].point, 0.5,
                                        points[n_r].point, 0.5);

    std::size_t id = points.size();
    Grid_point c(this, centroid, id, true/*on border*/);

    // need to push it to points (temporarily) for compute_closest_seed to work
    points.push_back(c);

    Grid_point& gp = points.back();
    gp.compute_closest_seed(n_q);
    gp.compute_closest_seed(n_r);
    gp.state = KNOWN;

    if(call_count > max_depth)
      return gp;

    if(gp.closest_seed_id == points[n_q].closest_seed_id)
      return Vor_vertex_on_edge(gp.index, n_r, ++call_count);
    else
    {
      CGAL_assertion(gp.closest_seed_id == points[n_r].closest_seed_id);
      return Vor_vertex_on_edge(gp.index, n_q, ++call_count);
    }
  }

  void visit_border_neighbors(const std::size_t n_p,
                             std::vector<bool>& visited_points)
  {
    const Grid_point& gp = points[n_p];

    visited_points[n_p] = true;
    std::size_t seed_p = gp.closest_seed_id;

    if(n_p < 4) // corner hack fixme
      Voronoi_vertices[seed_p].push_back(gp);

    Grid_point::Point_set::iterator it = gp.neighbors.begin();
    Grid_point::Point_set::iterator end = gp.neighbors.end();
    for(; it!=end; ++it)
    {
      const Grid_point& gq = points[*it];

      //fixme, this is not very robust, easy to find two border points whose
      // edge is not a border edge...
      if(!gq.is_on_domain_border)
        continue;

      if(!visited_points[gq.index])
      {
        // check if they have different colors
        std::size_t seed_q = gq.closest_seed_id;
        if(seed_p != seed_q)
        {
          Grid_point precise_vor_vertex = Vor_vertex_on_edge(n_p, gq.index);

          Voronoi_vertices[seed_p].push_back(precise_vor_vertex);
          Voronoi_vertices[seed_p].back().closest_seed_id = seed_p;
          Voronoi_vertices[seed_p].back().ancestor = gp.index;

          Voronoi_vertices[seed_q].push_back(precise_vor_vertex);
          Voronoi_vertices[seed_q].back().closest_seed_id = seed_q;
          Voronoi_vertices[seed_q].back().ancestor = gq.index;
        }

        // continue visiting from gq
        visit_border_neighbors(gq.index, visited_points);
      }
    }
  }

  void compute_precise_Voronoi_vertices_on_border()
  {
    // fixme : at the moment, we look for two vertices P and Q with different colors
    // (let's note them Cp and Cq). We then do some dichotomy-like search to find
    // a precise point equidistant from the seeds Sp and Sq by recursively
    // inserting virtual points on the segment [PQ]

    // BUT we force ancestry to stay on the segment, which is probably always the
    // case but we could theoretically have an ancestry that uses the third
    // point of the triangle...

    // MOREOEVER the color of that third point is unknown and might be :
    // - a third different color (in that case it should be handled by
    //   'compute_voronoi_vertex_in_triangle' anyway
    // - a color that is either Cp or Cq, and we can hope that it doesn't change
    //   too much from constraining the ancestry to stay on [PQ]

    // still need to be fixed though...

    // Related to that problem is the handling of a triangle with 3 colors with
    // a border edge... At the moment, we would compute two Voronoi vertices that
    // would be very close to each other... another fixme...

    std::size_t real_points_n = points.size();

    // find a first point on the border. (points[0] is usually one)
    std::vector<bool> visited_points(points.size(), false);
    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      if(points[i].is_on_domain_border)
      {
        visit_border_neighbors(i, visited_points);
        break;
      }
    }

    // clean off the virtual points created by Vor_vertex_on_edge
    shave_off_virtual_points(real_points_n);
  }

  void dual_shenanigans(Grid_point* gp)
  {
    CGAL_assertion(gp->state == KNOWN);

    std::list<std::size_t>::const_iterator it = gp->incident_triangles.begin(),
                                           end = gp->incident_triangles.end();
    for(; it!=end; ++it)
    {
      const Tri& tri = triangles[*it];
      Simplex dual_simplex;

      for(std::size_t j=0; j<3; ++j)
      {
        const Grid_point& gq = points[tri[j]];
        if(gq.state != KNOWN)
          continue;

        std::size_t q_id = gq.closest_seed_id;
        CGAL_assertion(q_id != static_cast<std::size_t>(-1));
        dual_simplex.insert(q_id);
//        fill_edge_bisectors_map(gp, gq);
      }
      add_simplex_to_triangulation(gp, dual_simplex);
    }
  }

  bool get_next_trial_point(Grid_point*& gp)
  {
    while(!trial_points.empty())
    {
      gp = trial_points.front();
      std::pop_heap(trial_points.begin(), trial_points.end(),
                    Grid_point_comparer<Grid_point>());
      trial_points.pop_back();

      CGAL_assertion(gp);

      if(gp->state != TRIAL)
      {
        std::cout << "WARNING : point with a state non-TRIAL in the PQ : ";
        std::cout << gp->index << "... Ignoring it!" << std::endl;
      }

      return true;
    }
    return false;
  }

  void spread_distances(const bool use_dual_shenanigans)
  {
    // fixme above : it's pointless to recompute the Voronoi vertices every time
    // but need a good criterion to decide when we need it (typically the last
    // time we spread_distance before optimizing)...

#if (verbose > 5)
    std::cout << "main loop" << std::endl;
#endif
    Grid_point* gp = NULL;
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

      if(!get_next_trial_point(gp))
        break;

#if (verbose > 10)
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << ", " << gp->point.y() << ") ";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from the seed "
                << gp->closest_seed_id << " ancestor : " << gp->ancestor << std::endl;
#endif

      gp->state = KNOWN;
      known_count++; // tmp --> use change_state()

      if(false && use_dual_shenanigans) // tmp
        dual_shenanigans(gp);

      PQ_state pqs = gp->update_neighbors_distances(trial_points);
      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
      is_t_empty = trial_points.empty();
    }

    debug();

    std::cout << "End of spread_distances. time: ";
    std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;

    // print_profiler_info();
  }

  void debug()
  {
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
        std::cout << "asd: " << points[gp.ancestor].children.size() << " " << points[gp.ancestor].closest_seed_id << std::endl;
      }

      if(gp.distance_to_closest_seed == FT_inf || gp.closest_seed_id >= seeds.size() ||
         (gp.ancestor != static_cast<std::size_t>(-1) &&
          points[gp.ancestor].closest_seed_id != gp.closest_seed_id))
      {
        failed = true;
        std::cout << "debug: " << i << " [" << gp.point << "] " << " at distance : ";
        std::cout << gp.distance_to_closest_seed << " from " << gp.closest_seed_id;
        std::cout << " ancestor: " << gp.ancestor;
        std::cout << " ancid: ";
        if(gp.ancestor==static_cast<std::size_t>(-1))
           std::cout << "-1" << std::endl;
        else
          std::cout << points[gp.ancestor].closest_seed_id << std::endl;
        std::cout << "this failing point is found in the neighbors of : ";
        for(std::size_t j=0, ps=points.size(); j<ps; ++j)
        {
          const Grid_point& gq = points[j];
          if(gq.neighbors.find(gp.index) != gq.neighbors.end())
            std::cout << j << " ";
        }
        std::cout << std::endl;
      }
    }
    if(failed)
      CGAL_assertion(false);

    std::cout << "End of debug. time: ";
    std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void add_to_neighboring_map(const Edge& e,
                              std::map<std::size_t,
                                       boost::unordered_set<std::size_t> >& m)
  {
    typedef std::map<std::size_t, boost::unordered_set<std::size_t> >       Map;

    std::size_t min = e[0], max = e[1];
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

  std::size_t ancestor_with_depth(const Grid_point* gp)
  {
    CGAL_precondition(gp);
    std::size_t depth = gp->depth;

    if(depth == 0)
      return -1;

    CGAL_precondition(depth > 0);
    std::size_t anc = gp->ancestor;
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

  std::deque<std::pair<std::size_t, FT> >
  simplify_geodesic(const std::deque<std::pair<std::size_t, FT> >& geodesic_path)
  {
    std::size_t max_point_n = 50;
    FT step = 1./static_cast<FT>(max_point_n - 1);
    FT curr_time = 0.;
    std::deque<std::pair<std::size_t, FT> > simplified_geodesic_path;

    // first point is naturally in
    simplified_geodesic_path.push_back(geodesic_path[0]);

    for(std::size_t i=0; i<geodesic_path.size(); ++i)
    {
      const std::pair<std::size_t, FT>& p = geodesic_path[i];
      const FT time = p.second;

      if(time - curr_time >= step)
      {
        // this point is sufficiently far from the previous one
        simplified_geodesic_path.push_back(p);
        curr_time = time;
      }
    }

    // need to make sure we have the final point of the geodesic in
    if(simplified_geodesic_path.back() != geodesic_path.back())
      simplified_geodesic_path.push_back(geodesic_path.back());

    std::cout << "simplified from " << geodesic_path.size()
              << " to " << simplified_geodesic_path.size() << std::endl;

    return simplified_geodesic_path;
  }

  void extract_geodesic(const Grid_point* gp,
                        std::vector<std::deque<std::pair<std::size_t, FT> > >& geodesics)
  {
    std::deque<std::pair<std::size_t, FT> > geodesic_path;
    FT total_length = gp->distance_to_closest_seed;
    geodesic_path.push_front(std::make_pair(gp->index, 1.));
    CGAL_assertion(total_length != 0.);
    FT total_length_inv = 1./total_length;

#define EXTRACT_WITH_DEPTH
#ifdef EXTRACT_WITH_DEPTH
    std::size_t anc = ancestor_with_depth(gp);
#else
    std::size_t anc = gp->ancestor;
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
      anc = ancestor_with_depth(&(points[anc]));
#else
      anc = points[anc].ancestor;
#endif
    }

    std::size_t s1 = gp->closest_seed_id;
    int s2 = gp->is_seed_holder;

    CGAL_precondition(s2 >= 0);
    CGAL_precondition(s1 < static_cast<std::size_t>(s2));

    std::cout << "built geodesic of " << s1 << " " << s2 << std::endl;
    geodesics.push_back(geodesic_path);
  }

  void spread_distances_from_one_seed(const std::size_t seed_id,
                                      std::map<std::size_t,
                                               boost::unordered_set<std::size_t> >& neighbors,
                                      std::vector<std::deque<std::pair<std::size_t, FT> > >& geodesics,
                                      std::vector<Grid_point>& grid_point_memory)
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
      Grid_point& gp = points[i];

      if(static_cast<std::size_t>(gp.is_seed_holder) == seed_id)
      {
        CGAL_assertion(gp.closest_seed_id == seed_id);
        CGAL_assertion(gp.ancestor == static_cast<std::size_t>(-1));
        CGAL_assertion(gp.depth == 0);

        gp.state = TRIAL;
        trial_points.push_back(&gp);
        std::push_heap(trial_points.begin(), trial_points.end(),
                       Grid_point_comparer<Grid_point>());
        break;
      }
    }

    bool is_t_empty = trial_points.empty();
    CGAL_precondition(!is_t_empty);

    known_count = 0;
    Grid_point* gp = NULL;
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

      if(!get_next_trial_point(gp))
        break;

#if (verbose > 10)
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << ", " << gp->point.y() << ") ";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from the seed "
                << gp->closest_seed_id << " ancestor : " << gp->ancestor << std::endl;
#endif

      // Check if we have reached one of the required seeds
      if(gp->is_seed_holder >= 0) // really ugly, fixme (see comment at member def)
      {
        std::cout << "found the seed " << gp->is_seed_holder << " at " << gp->index << " (" << gp->point << ") " << std::endl;

        // That's the center of a Voronoi cell (no ancestor)
        std::size_t reached_seed_id = gp->is_seed_holder;
        CGAL_assertion(gp->closest_seed_id == seed_id);

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
          extract_geodesic(gp, geodesics);
          seeds_to_reach.erase(sit);
        }
      }

      // only spread as little as possible
      if(seeds_to_reach.empty())
        break;

      gp->state = KNOWN;
      known_count++; // tmp --> use change_state()

      // Change the distance_to_closest_seed of points to allow for seed_id's
      // cell to spread

      typename Grid_point::Point_set::const_iterator it = gp->neighbors.begin(),
                                                     end = gp->neighbors.end();
      for(; it!=end; ++it)
      {
        Grid_point& gq = points[*it];

        if(gq.state == FAR)
        {
          grid_point_memory.push_back(gq);
          gq.distance_to_closest_seed = FT_inf;
        }
      }

      gp->update_neighbors_distances(trial_points);

      // always rebuild the priority queue
      std::make_heap(trial_points.begin(), trial_points.end(),
                     Grid_point_comparer<Grid_point>());
      is_t_empty = trial_points.empty();
    }

    std::cout << "reached all the seeds in " << known_count << " points" << std::endl;

    // make sure we've reached all the neighboring seeds
    CGAL_postcondition(seeds_to_reach.empty());
  }

  void rollback_points(const std::vector<Grid_point>& gp_memory)
  {
    std::cout << "rollback" << std::endl;
    for(std::size_t i=0, gpms=gp_memory.size(); i<gpms; ++i)
    {
      const Grid_point& gp_m = gp_memory[i];
      std::size_t id = gp_m.index;
      CGAL_assertion(points[id].point == gp_m.point);
      points[id] = gp_m;
    }

    trial_points.clear();
  }

  void output_geodesics(const std::vector<std::deque<std::pair<std::size_t, FT> > >& geodesics,
                        const std::string str_base)
  {
    if(geodesics.empty())
    {
      std::cout << "no geodesic to output" << std::endl;
      return;
    }

    std::ofstream out((str_base + "_geodesics.mesh").c_str());

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';

    out << "Vertices" << '\n';
    out << points.size() << '\n';
    for(std::size_t i=0; i<points.size(); ++i)
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
        out << p.first + 1 << " " << pp1.first + 1 << " " << i << std::endl;
      }
    }

    out << "End" << std::endl;
  }

  std::size_t combi(std::size_t n, std::size_t k)
  {
    // n choose k
    // this will (silently) overflow for large n...

    if (k > n)
      return 0;

    std::size_t r = 1;
    for(std::size_t i=1; i<=k; ++i)
    {
      r *= n--;
      r /= i;
    }

    return r;
  }

  FT evaluate_Berstein_poly(const std::size_t i, const std::size_t n,
                            const  FT t)
  {
    // computes (i n) (1-t)^{n-i} t^i

    std::size_t c = combi(n, i);
    FT a = std::pow((1-t), n-i);
    FT b = std::pow(t, i);

//    std::cout << "eval berstein poly : " << i << " " << n << " " << t << std::endl;
//    std::cout << "c: " << c << " a: " << a << " b: " << b << " cab: " << c*a*b << std::endl;

    return c * a * b;
  }

  Point_2 evaluate_Bezier_curve(const std::vector<Point_2>& control_points, const  FT t)
  {
    // using de casteljau's algorithm
    std::size_t number_of_control_points = control_points.size();
    std::vector<FT> cp_xs(number_of_control_points);
    std::vector<FT> cp_ys(number_of_control_points);

    for(std::size_t i=0; i<number_of_control_points; ++i)
    {
      cp_xs[i] = control_points[i].x();
      cp_ys[i] = control_points[i].y();
    }

    for(std::size_t i=1; i<number_of_control_points; ++i)
    {
      for(std::size_t j=0; j<number_of_control_points - i; ++j)
      {
        cp_xs[j] = (1-t) * cp_xs[j] + t * cp_xs[j+1];
        cp_ys[j] = (1-t) * cp_ys[j] + t * cp_ys[j+1];
      }
    }

    return Point_2(cp_xs[0], cp_ys[0]);
  }

  void compute_control_points(const std::deque<std::pair<std::size_t, FT> >& geodesic,
                              std::vector<Point_2>& control_points)
  {
    // the deque is made of points and a time between [0,1] obtained by knowing
    // that the speed is constant along a geodesic

    // A Bézier curve is B(t) = \sum_i (i n) (1-t)^{n-i} t^i P_i
    // and we know B(t0), \dots, B(tn). This is thus solving a system AX=B

    // the construction is definitely not the most efficient...
    std::size_t geo_path_size = geodesic.size();
    Eigen::MatrixXd A(geo_path_size, geo_path_size);
    Eigen::VectorXd B(geo_path_size);
    Eigen::VectorXd X(geo_path_size), Y(geo_path_size);

    // solve for the x coefficients
    std::size_t row_id = 0.;
    for(; row_id < geo_path_size; ++row_id)
    {
      const std::pair<std::size_t, FT>& pair = geodesic[row_id];
      const Point_2& p = points[pair.first].point;
      const FT time = pair.second;

      B(row_id) = p.x();
      std::cout << "set_b@x: " << row_id << " " << p.x() << std::endl;

      for(std::size_t col_id=0; col_id < geo_path_size; ++col_id)
      {
        FT e = evaluate_Berstein_poly(col_id, geo_path_size-1, time);
        std::cout << "set_a: " << row_id << " " << col_id << " " << e << std::endl;
        A(row_id, col_id) = e;
      }
    }
    X = A.fullPivLu().solve(B);

    // solve for the y coefficients (the matrix A stays the same)
    row_id = 0.;
    for(; row_id < geo_path_size; ++row_id)
    {
      const std::pair<std::size_t, FT>& pair = geodesic[row_id];
      const Point_2& p = points[pair.first].point;

      B(row_id) = p.y();
      std::cout << "set_b@y: " << row_id << " " << p.y() << std::endl;
    }
    Y = A.fullPivLu().solve(B);

    // combine to obtain the control points !
    CGAL_precondition(control_points.size() == geodesic.size());

    std::cout << "solutions : " << std::endl << X.transpose() << std::endl << Y.transpose() << std::endl;

    for(std::size_t i=0; i<geo_path_size; ++i)
    {
      Point_2 p(X(i), Y(i));

//      std::cout << "control point i: " << i << " is " << p << std::endl;
      control_points[i] = p;
    }

    // debug (make sure that we go through the points that we wanted to interpolate)
    for(std::size_t i=0; i<geodesic.size(); ++i)
    {
      const std::pair<std::size_t, FT>& pair = geodesic[i];
      const Point_2& p = evaluate_Bezier_curve(control_points, pair.second);

//      std::cout << "i/t: " << i << " " << pair.second << std::endl;
//      std::cout << "compare : " << p << " and " << points[pair.first].point << std::endl;
      CGAL_assertion(CGAL::squared_distance(p, points[pair.first].point) < 1e-8);
    }
  }

  void compute_Bezier_curve(const std::deque<std::pair<std::size_t, FT> >& geodesic,
                            std::vector<Point_2>& control_points)
  {
    // todo simplify some notations maybe...
    std::cout << "compute Bezier approximation of the geodesic between: ";
    std::cout << points[geodesic.front().first].closest_seed_id << " and ";
    std::cout << points[geodesic.back().first].closest_seed_id << std::endl;

    // Before computing the control points, we might want to split the geodesic
    // in small segments of n points to limit the degree of the interpolating
    // Bezier curve
    CGAL_precondition(geodesic.size() == control_points.size());
    CGAL_precondition(geodesic.size() >= 2);

    std::size_t size_of_geodesic_segment = 4;
    std::size_t size_of_geodesic_path = geodesic.size();
    std::size_t start_of_segment = 0;
    std::size_t remaining_points = size_of_geodesic_path - 1;

    while(start_of_segment < size_of_geodesic_path - 1)
    {
      CGAL_assertion(remaining_points > 0); // can't make a segment if there is no remaining point
      if(remaining_points < size_of_geodesic_segment - 1) // can't make a full segment
        size_of_geodesic_segment = 1 + remaining_points;

      CGAL_assertion(size_of_geodesic_segment > 1); // proper segment

      std::cout << "total length of the path : " << size_of_geodesic_path << std::endl;
      std::cout << "Segment from : " << start_of_segment << " and size: " << size_of_geodesic_segment << std::endl;
      CGAL_assertion(start_of_segment + size_of_geodesic_segment - 1 < size_of_geodesic_path);

      // We want to go through the first and last point of the segment, thus we need
      // to locally scale up the time associated with these points
      FT start_of_segment_t = geodesic[start_of_segment].second;
      FT end_of_segment_t = geodesic[start_of_segment + size_of_geodesic_segment - 1].second;
      FT diff_t = end_of_segment_t - start_of_segment_t;
      FT scaling = 1. / diff_t;

      // fill the segment (todo make it lighter without copy)
      std::deque<std::pair<std::size_t, FT> > geodesic_segment;


      std::cout << "geodesic segment : " << std::endl;
      for(std::size_t i=0; i<size_of_geodesic_segment; ++i)
      {
        // same point but the time is scaled
        geodesic_segment.push_back(std::make_pair(geodesic[start_of_segment + i].first,
            scaling * (geodesic[start_of_segment + i].second - start_of_segment_t)));
        std::cout << "init: " << points[geodesic[start_of_segment + i].first].point
                  << " and time: " << geodesic[start_of_segment + i].second << std::endl;
        std::cout << "local : " << points[geodesic_segment[i].first].point
                  << " and time: " << geodesic_segment[i].second << std::endl;
      }

      std::vector<Point_2> segment_control_points(size_of_geodesic_segment);
      compute_control_points(geodesic_segment, segment_control_points);

      std::cout << "geodesic control points: " << std::endl;
      for(std::size_t i=0; i<size_of_geodesic_segment; ++i)
      {
        CGAL_assertion(start_of_segment + i < control_points.size());
        control_points[start_of_segment + i] = segment_control_points[i];
        std::cout << segment_control_points[i] << std::endl;
      }

      CGAL_postcondition(CGAL::squared_distance(control_points[start_of_segment],
                                                points[geodesic[start_of_segment].first].point) < 1e-8);
      CGAL_postcondition(CGAL::squared_distance(control_points[start_of_segment + size_of_geodesic_segment - 1],
                                                points[geodesic[start_of_segment + size_of_geodesic_segment - 1].first].point) < 1e-8);

      // -1 because the next segment starts at the end of the previous one
      start_of_segment += size_of_geodesic_segment - 1;

      // -1 because it's not useful to count the next starting point as part of the remaining points
      remaining_points -= size_of_geodesic_segment - 1;
    }

    for(std::size_t i=0; i<size_of_geodesic_path; ++i)
    {
      CGAL_postcondition(control_points[i] != Point_2());
    }
  }

  void output_Bezier_curves(const std::vector<std::vector<Point_2> >& control_points,
                            std::size_t offset = 0,
                            const bool draw_control_points = false)
  {
    // fixme if the geodesics are approximated by multiple low level Bezier curves
    // then what is drawn should be each curve independantly, not treating the
    // list of control points of different Bezier curves as if they were the
    // control points of a single high degree Bezier curve

    std::ofstream out("Bezier_test.mesh");
    std::size_t number_of_sample_points = 1000;
    FT step = 1./static_cast<FT>(number_of_sample_points - 1.);
    std::size_t number_of_curves = control_points.size();

    std::size_t total_vertices_n = number_of_curves * number_of_sample_points;
    std::size_t total_edges_n = number_of_curves * (number_of_sample_points - 1);

    if(draw_control_points)
    {
      for(std::size_t i=0; i<number_of_curves; ++i)
      {
        total_vertices_n += control_points[i].size();
        total_edges_n += control_points[i].size()-1;
      }
    }

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';

    out << "Vertices" << '\n';
    out << total_vertices_n << std::endl;

    for(std::size_t i=0; i<number_of_curves; ++i)
    {
      const std::vector<Point_2>& local_control_points = control_points[i];

      if(draw_control_points)
      {
        for(std::size_t j=0; j<local_control_points.size(); ++j)
          out << local_control_points[j] << " " << i << std::endl;
      }

      for(std::size_t j=0; j<number_of_sample_points; ++j)
      {
        out << evaluate_Bezier_curve(local_control_points, step * static_cast<FT>(j))
            << " " << i << std::endl;
      }
    }

    out << "Edges" << '\n';
    out << total_edges_n << '\n';
    for(std::size_t i=0; i<number_of_curves; ++i)
    {
      if(draw_control_points)
      {
        std::size_t number_of_local_control_points = control_points[i].size();

        // edges to draw the control points
        for(std::size_t j=0; j<number_of_local_control_points-1; ++j)
          out << offset+j+1 << " " << offset+j+2 << " " << i << '\n'; // +1 for medit

        offset += number_of_local_control_points;
      }

      // edges to draw the bezier curve
      for(std::size_t j=0; j<number_of_sample_points-1; ++j)
        out << offset+j+1 << " " << offset+j+2 << " " << i << '\n'; // +1 for medit

      offset += number_of_sample_points;
    }

    out << "End" << std::endl;
  }

  void approximate_geodesics_with_Bezier(std::vector<std::deque<std::pair<std::size_t, FT> > >& geodesics)
  {
    std::vector<std::vector<Point_2> > control_points;

    for(std::size_t i=0, gs=geodesics.size(); i<gs; ++i)
    {
      const std::deque<std::pair<std::size_t, FT> >& geodesic = geodesics[i];
      const std::size_t geo_size = geodesic.size();

      std::vector<Point_2> single_geodesic_control_points(geo_size);
      compute_Bezier_curve(geodesic, single_geodesic_control_points);

      std::cout << "full geodesic :" << std::endl;
      for(std::size_t j=0; j<geo_size; ++j)
      {
        const Point_2& p = points[geodesic[j].first].point;
        std::cout << "{" << p.x() << ", " << p.y() << "}, " <<  std::endl;
      }

      std::cout << "full control points: " << std::endl;
      for(std::size_t j=0; j<geo_size; ++j)
      {
        const Point_2& p = single_geodesic_control_points[j];
        std::cout << "{" << p.x() << ", " << p.y() << "}, " <<  std::endl;
      }

      control_points.push_back(single_geodesic_control_points);
    }

    CGAL_postcondition(geodesics.size() == control_points.size());

    output_Bezier_curves(control_points,
                         0, //points.size() /*offset*/,
                         false/*draw control points*/);
  }

  FT compute_Riemannian_angle(
      std::size_t origin,
      std::size_t extr_id_1,
      std::size_t extr_id_2,
      const std::vector<std::deque<std::pair<std::size_t/*id*/,
                                             FT/*time*/> > >& geodesics)
  {
    // computes the angle between the tangents of the geodesics
    // ori -- extr_1 and ori -- extr_2

    std::cout << "orign: " << origin << std::endl;
    std::cout << "t1: " << extr_id_1 << std::endl;
    std::cout << "t2: " << extr_id_2 << std::endl;

    // first, find the geodesics in the geodesics map
    // geodesics were only draw from id1 to id2 with id1 < id2 so we need to know
    // if it's at the end of a geodesic
    std::size_t look_up_orig_1 = origin, look_up_orig_2 = origin;
    std::size_t look_up_end_1 = extr_id_1, look_up_end_2 = extr_id_2;

    // now find the corresponding geodesics
    // and extract the corresponding segments and points
    std::size_t ori_1 = -1, ori_2 = -1; // second one not really needed
    std::size_t tangent_extr_1 = -1, tangent_extr_2 = -1;

    // warning: disgustingly unefficient
    // proper way would be to have geodesics to be a map[id_1, id_2] = deque of pairs<id, time>
    for(std::size_t i=0; i<geodesics.size(); ++i)
    {
      const std::deque<std::pair<std::size_t, FT> >& geo = geodesics[i];
      std::size_t gs = geo.size();

/*
      std::cout << "considering geodesic between ";
      std::cout << points[geo[0].first].is_seed_holder << " and ";
      std::cout << points[geo[gs-1].first].is_seed_holder << std::endl;

      std::cout << geo[0].first << " " << geo[1].first << " ";
      std::cout << geo[gs-2].first << " " << geo[gs-1].first << std::endl;
*/
      if(points[geo.front().first].is_seed_holder == look_up_orig_1 &&
         points[geo.back().first].is_seed_holder == look_up_end_1)
      {
        ori_1 = geo[0].first;
        tangent_extr_1 = geo[1].first;
      }
      else if(points[geo.front().first].is_seed_holder == look_up_end_1 &&
              points[geo.back().first].is_seed_holder == look_up_orig_1)
      {
        ori_1 = geo[gs-1].first;
        tangent_extr_1 = geo[gs-2].first;
      }

      if(points[geo.front().first].is_seed_holder == look_up_orig_2 &&
         points[geo.back().first].is_seed_holder == look_up_end_2)
      {
        ori_2 = geo[0].first;
        tangent_extr_2 = geo[1].first;
      }
      else if(points[geo.front().first].is_seed_holder == look_up_end_2 &&
              points[geo.back().first].is_seed_holder == look_up_orig_2)
      {
        ori_2 = geo[gs-1].first;
        tangent_extr_2 = geo[gs-2].first;
      }

      if(ori_1 != -1 && ori_2 != -1)
        break; // found both geodesics
    }

    if(tangent_extr_1 == -1 || tangent_extr_2 == -1)
      std::cout << " didn't find the geodesics... " << std::endl;

    if(ori_1 != ori_2 && ori_1 != -1)
      std::cout << "screwed up the origins" << std::endl;

    std::cout << "ori: " << ori_1 << " " << ori_2 << std::endl;
    std::cout << "seed holders: " << points[ori_1].is_seed_holder << std::endl;
    std::cout << "t1: " << tangent_extr_1 << std::endl;
    std::cout << "t2: " << tangent_extr_2 << std::endl;

    if(tangent_extr_1 == tangent_extr_2)
    {
      std::cout << "degenerate case with equal tangent segments" << std::endl;
      return 0.;
    }

    // compute the angle and make sure it's within 0 pi
    Point_2 o_p = points[ori_1].point;
    Point_2 t1_p = points[tangent_extr_1].point;
    Point_2 t2_p = points[tangent_extr_2].point;

    const Metric& m = points[ori_1].metric;
    Point_2 t_o_p = m.transform(o_p);
    Point_2 t_t1p = m.transform(t1_p);
    Point_2 t_t2p = m.transform(t2_p);

    Vector_2 v1 = t_t1p - t_o_p;
    Vector_2 v2 = t_t2p - t_o_p;

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
                        const char* filename)
  {
    std::cout << "output: " << filename << std::endl;
    std::ofstream out(filename);
    std::size_t histo_n = histogram.size();
    std::cout << "histo_n: " << histo_n << std::endl;
    for(std::size_t i=0; i<histo_n; ++i)
    {
      FT val = min + (max-min)*((FT) i)/((FT) histo_n);
      out << i << "," << val << "," << histogram[i] << std::endl;
    }
  }

  void output_histogram(std::vector<FT>& values,
                        const char* filename)
  {
    FT min_value = *(std::min_element(values.begin(), values.end()));
    FT max_value = *(std::max_element(values.begin(), values.end()));

    std::cout << "Outputing: " << values.size() << " " << min_value << " " << max_value << std::endl;

    int histogram_size = 1000;
    std::vector<int> histogram(histogram_size, 0);
    FT limit_val = histogram_size - 1.;
    FT step_size = (max_value - min_value) / (FT) histogram_size;

    for(std::size_t i=0; i<values.size(); ++i)
      histogram[(std::min)(limit_val, std::floor((values[i]-min_value)/step_size))]++;

    output_histogram(histogram, min_value, max_value, filename);
  }

  void compute_Riemannian_angles(
         const std::vector<std::deque<std::pair<std::size_t/*id*/,
                                                FT/*time*/> > >& geodesics,
         const std::string str_base)
  {
    std::cout << "computing Riemannian angles" << std::endl;

    // Evaluate the angles of the Riemannian simplices

    // The angles are computed at the vertex in the metric of the vertex
    // and are defined as the angle between the tangent vectors of the
    // geodesics

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

        std::cout << "in triangle: " << tr[0] << " " << tr[1] << " " << tr[2] << std::endl;

        FT angle = compute_Riemannian_angle(orig, dest_1, dest_2, geodesics);
        values.push_back(angle);
      }
    }
    output_histogram(values, "histogram_riemannian_face_angles.cvs");
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

    std::vector<std::deque<std::pair<std::size_t/*id*/,FT/*time*/> > > geodesics;

    for(std::size_t seed_id=0, ss=seeds.size(); seed_id<ss; ++seed_id)
    {
      std::cout << "computing geodesics from : " << seed_id << std::endl;

      // below is kind of ugly... but it's clearer than switching everything to
      // 'FAR' earlier and then having the memory overwrite 'KNOWN' with 'FAR'
      // every loop iteration
      for(std::size_t i=0, ptss=points.size(); i<ptss; ++i)
        points[i].state = FAR;

      // since we're spreading from a seed over the other cells, we must keep
      // the previous colors in memory
      std::vector<Grid_point> grid_point_memory;

      spread_distances_from_one_seed(seed_id, neighbors, geodesics,
                                     grid_point_memory);
      rollback_points(grid_point_memory); // reset the points we have overwritten
    }

    compute_Riemannian_angles(geodesics, str_base);
    output_geodesics(geodesics, str_base);

    // simplify geodesics
    for(std::size_t i=0; i<geodesics.size(); ++i)
    {
//      geodesics[i] = simplify_geodesic(geodesics[i]);
    }

//    approximate_geodesics_with_Bezier(geodesics);
  }

  std::size_t probe_edge_midpoints_map(const std::size_t n_p,
                                       const std::size_t n_q,
                                       boost::unordered_map<std::pair<std::size_t, std::size_t>,
                                       std::size_t>& edge_mid_points,
                                       bool create_if_not_found = true)
  {
    typedef boost::unordered_map<std::pair<std::size_t, std::size_t>,
                                 std::size_t>                     Edge_mpts_map;

    std::size_t n_pq = -1;

    bool ordered = (n_p < n_q);
    int nmin = ordered ? n_p : n_q;
    int nmax = ordered ? n_q : n_p;

    std::pair<typename Edge_mpts_map::iterator, bool> is_insert_successful =
        edge_mid_points.insert(std::make_pair(std::make_pair(nmin, nmax), -1)); // -1 is a placeholder

    if(!is_insert_successful.second) // already exists in the map
      n_pq = is_insert_successful.first->second;
    else if(create_if_not_found) // the mid point hasn't been computed yet and we want to build it
    {
      Point_2 mid_point = CGAL::barycenter(points[n_p].point, 0.5,
                                           points[n_q].point, 0.5);
      n_pq = points.size();
      bool border_info = (points[n_p].is_on_domain_border && points[n_q].is_on_domain_border);

      Grid_point gm(this, mid_point, n_pq, border_info);
      points.push_back(gm);
      is_insert_successful.first->second = n_pq;
    }

    return n_pq;
  }

  std::size_t create_new_triangle(const std::size_t n_p, const std::size_t n_q,
                                  const std::size_t n_r)
  {
    Tri new_tri;
    new_tri[0] = n_p;
    new_tri[1] = n_q;
    new_tri[2] = n_r;
    triangles_buffer.push_back(new_tri);
    return triangles.size() + triangles_buffer.size() - 1;
  }

  void split_1(const std::size_t n_p, const std::size_t n_q, const std::size_t n_r,
               const std::size_t n_pq, const std::size_t i,
               boost::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t>& m)
  {
    CGAL_assertion(probe_edge_midpoints_map(n_p, n_q, m, false) == n_pq);

    // n_p, n_q & n_r belong to i but they are not necessarily the 0th, 1st, and 2nd coordinates
    Tri& triangle = triangles[i];
    Grid_point& gp = points[n_p];
    Grid_point& gq = points[n_q];
    Grid_point& gr = points[n_r];

    // p and r have the new i as incident triangle so only needs to be removed in q
    gq.remove_from_incident_triangles(i);

    triangle[0] = n_p;
    triangle[1] = n_pq;
    triangle[2] = n_r;

    std::size_t new_tri_id = create_new_triangle(n_pq, n_q, n_r);

    gr.incident_triangles.push_back(new_tri_id);
    gq.incident_triangles.push_back(new_tri_id);
    points[n_pq].incident_triangles.push_back(new_tri_id);
    points[n_pq].incident_triangles.push_back(i);

    // neighbors
    gp.remove_from_neighbors(n_q);
    gp.neighbors.insert(n_pq);

    gq.remove_from_neighbors(n_q);
    gq.neighbors.insert(n_pq);

    gr.neighbors.insert(n_pq);

    points[n_pq].neighbors.insert(n_r); // n_p & n_q already added when n_pq was created
  }

  void split_2(const std::size_t n_p, const std::size_t n_q, const std::size_t n_r,
               const std::size_t n_pq, const std::size_t n_qr, const std::size_t i,
               boost::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t>& m)
  {
    CGAL_assertion(probe_edge_midpoints_map(n_p, n_q, m, false) == n_pq);
    CGAL_assertion(probe_edge_midpoints_map(n_q, n_r, m, false) == n_qr);

    // n_p, n_q & n_r belong to i but they are not necessarily the 0th, 1st, and 2nd coordinates
    Tri& triangle = triangles[i];
    Grid_point& gp = points[n_p];
    Grid_point& gq = points[n_q];
    Grid_point& gr = points[n_r];

    // q has the new i as incident triangle so removed in p & r
    gp.remove_from_incident_triangles(i);
    gr.remove_from_incident_triangles(i);

    triangle[0] = n_pq;
    triangle[1] = n_q;
    triangle[2] = n_qr;

    std::size_t new_tri_1_id = create_new_triangle(n_p, n_pq, n_r);
    std::size_t new_tri_2_id = create_new_triangle(n_pq, n_qr, n_r);

    gp.incident_triangles.push_back(new_tri_1_id);
    gr.incident_triangles.push_back(new_tri_1_id);
    gr.incident_triangles.push_back(new_tri_2_id);
    points[n_pq].incident_triangles.push_back(new_tri_1_id);
    points[n_pq].incident_triangles.push_back(new_tri_2_id);
    points[n_pq].incident_triangles.push_back(i);
    points[n_qr].incident_triangles.push_back(new_tri_2_id);
    points[n_qr].incident_triangles.push_back(i);

    // neighbors
    gp.remove_from_neighbors(n_q);
    gp.neighbors.insert(n_pq);

    gq.remove_from_neighbors(n_p);
    gq.remove_from_neighbors(n_r);
    gq.neighbors.insert(n_pq);
    gq.neighbors.insert(n_qr);

    gr.remove_from_neighbors(n_q);
    gr.neighbors.insert(n_pq);

    points[n_pq].neighbors.insert(n_r);
    points[n_pq].neighbors.insert(n_qr);

    points[n_qr].neighbors.insert(n_pq);
  }

  void split_3(const std::size_t n_p, const std::size_t n_q, const std::size_t n_r,
               const std::size_t n_pq, const std::size_t n_qr, const std::size_t n_pr,
               const std::size_t i,
               boost::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t>& m)
  {
    CGAL_assertion(probe_edge_midpoints_map(n_p, n_q, m, false) == n_pq);
    CGAL_assertion(probe_edge_midpoints_map(n_q, n_r, m, false) == n_qr);
    CGAL_assertion(probe_edge_midpoints_map(n_p, n_r, m, false) == n_pr);

    // n_p, n_q & n_r belong to i but they are not necessarily the 0th, 1st, and 2nd coordinates
    Tri& triangle = triangles[i];
    Grid_point& gp = points[n_p];
    Grid_point& gq = points[n_q];
    Grid_point& gr = points[n_r];
    gp.remove_from_incident_triangles(i);
    gq.remove_from_incident_triangles(i);
    gr.remove_from_incident_triangles(i);

    triangle[0] = n_pq;
    triangle[1] = n_qr;
    triangle[2] = n_pr;

    std::size_t new_tri_1_id = create_new_triangle(n_p, n_pq, n_pr);
    std::size_t new_tri_2_id = create_new_triangle(n_pq, n_q, n_qr);
    std::size_t new_tri_3_id = create_new_triangle(n_qr, n_r, n_pr);

    gp.incident_triangles.push_back(new_tri_1_id);
    gq.incident_triangles.push_back(new_tri_2_id);
    gr.incident_triangles.push_back(new_tri_3_id);

    points[n_pq].incident_triangles.push_back(new_tri_1_id);
    points[n_pq].incident_triangles.push_back(new_tri_2_id);
    points[n_pq].incident_triangles.push_back(i);

    points[n_qr].incident_triangles.push_back(new_tri_3_id);
    points[n_qr].incident_triangles.push_back(new_tri_2_id);
    points[n_qr].incident_triangles.push_back(i);

    points[n_pr].incident_triangles.push_back(new_tri_1_id);
    points[n_pr].incident_triangles.push_back(new_tri_3_id);
    points[n_pr].incident_triangles.push_back(i);

    // neighbors
    gp.remove_from_neighbors(n_q);
    gp.remove_from_neighbors(n_r);
    gp.neighbors.insert(n_pq);
    gp.neighbors.insert(n_pr);

    gq.remove_from_neighbors(n_p);
    gq.remove_from_neighbors(n_r);
    gq.neighbors.insert(n_pq);
    gq.neighbors.insert(n_qr);

    gr.remove_from_neighbors(n_p);
    gr.remove_from_neighbors(n_q);
    gr.neighbors.insert(n_qr);
    gr.neighbors.insert(n_pr);

    points[n_pq].neighbors.insert(n_qr);
    points[n_pq].neighbors.insert(n_pr);

    points[n_qr].neighbors.insert(n_pq);
    points[n_qr].neighbors.insert(n_pr);

    points[n_pr].neighbors.insert(n_pq);
    points[n_pr].neighbors.insert(n_qr);
  }

  void refine_grid()
  {
    bool needs_to_be_refined = false;
    do
    {
      needs_to_be_refined = refine_grid_one_step();
    } while(needs_to_be_refined);
  }

  bool refine_grid_one_step()
  {
    std::cout << "refine grid one step" << std::endl;

    // tedious function so it's written in a very clear way... it's (as usual)
    // not the most efficient but the most straight forward

    // can't use ref or pointers to 'points' since they are regularly invalidated
    // by point insertions

    // todo
    // make it nice with subfunctions and efficient stuff when cleaning up the code

    boost::unordered_set<std::size_t> seeds_to_refine;
    boost::unordered_set<std::size_t> seeds_to_rebuild;
    std::vector<bool> refined_triangles(triangles.size(), false);
    boost::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t> edge_mid_points;

    // determine which seeds need to be refined
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& triangle = triangles[i];

      if(points[triangle[0]].closest_seed_id == points[triangle[1]].closest_seed_id &&
         points[triangle[1]].closest_seed_id == points[triangle[2]].closest_seed_id)
        continue;
      // here and below, 'triangle' is on a bisector (of dimension >= 1)

      for(int j=0; j<3; ++j)
      {
        const Grid_point& gp = points[triangle[j]];
        std::size_t ap_length = gp.ancestor_path_length;
        if(ap_length < min_ancestor_path_length)
        {
          seeds_to_refine.insert(gp.closest_seed_id);
          std::cout << "the canvas is thin (" << ap_length << ") at : "
                    << gp.index << " (" << gp.point << ")"
                    << " and cellid: " << gp.closest_seed_id << std::endl;
        }
      }
    }

    if(seeds_to_refine.empty())
      return false;

    std::cout << points.size() << " points before grid refinement" << '\n';
    std::cout << triangles.size() << " triangles before grid refinement" << std::endl;

    // base mesh refining, not seeds refining !
    for(std::size_t i=0, ts=triangles.size(); i < ts; ++i)
    {
      Tri& triangle = triangles[i];

      // the triangle must be refined if at least one of its colors is in seeds_to_refine

      // todo: a better algorithm would be to spread from a triangle for whom
      // a vertex needs a refinement (a geodesic Voronoi cell is convex)

      bool must_be_refined = false;
      for(int j=0; j<3; ++j)
      {
        const Grid_point& gp = points[triangle[j]];
        if(seeds_to_refine.find(gp.closest_seed_id) != seeds_to_refine.end())
        {
          must_be_refined = true;
          break;
        }
      }

      if(!must_be_refined)
        continue;

      refined_triangles[i] = true;

      // we use the triforce split to get similar triangles
      // the old index 'i' will correspond to the new triangle in the middle
      // and we build three new triangles

      // triforce split :
      //        a
      //       /|
      //      / |
      // f   /__| e
      //    /| /|
      // b /_|/_| c
      //     d

      // doing it explicitely (lengthy) because it's clearer than fancy loops
      const std::size_t n_p = triangle[0];
      const std::size_t n_q = triangle[1];
      const std::size_t n_r = triangle[2];
      CGAL_precondition(n_p < points.size() && n_q < points.size() && n_r < points.size());

      // clean off the neighbors, border neighbors and incident triangles info
      points[n_p].remove_from_incident_triangles(i);
      points[n_q].remove_from_incident_triangles(i);
      points[n_r].remove_from_incident_triangles(i);

      points[n_p].remove_from_neighbors(n_q);
      points[n_p].remove_from_neighbors(n_r);
      points[n_q].remove_from_neighbors(n_p);
      points[n_q].remove_from_neighbors(n_r);
      points[n_r].remove_from_neighbors(n_p);
      points[n_r].remove_from_neighbors(n_q);

      // keep the colors that we will need to paint again
      seeds_to_rebuild.insert(points[n_p].closest_seed_id);
      seeds_to_rebuild.insert(points[n_q].closest_seed_id);
      seeds_to_rebuild.insert(points[n_r].closest_seed_id);

      // find (or create) the edge points
      const std::size_t n_pq = probe_edge_midpoints_map(n_p, n_q, edge_mid_points);
      const std::size_t n_qr = probe_edge_midpoints_map(n_q, n_r, edge_mid_points);
      const std::size_t n_pr = probe_edge_midpoints_map(n_p, n_r, edge_mid_points);

      CGAL_postcondition(n_pq != static_cast<std::size_t>(-1) &&
                         n_qr != static_cast<std::size_t>(-1) &&
                         n_pr != static_cast<std::size_t>(-1));

      // create the three new triangles
      std::size_t new_tri_1_id = create_new_triangle(n_p, n_pq, n_pr);
      std::size_t new_tri_2_id = create_new_triangle(n_q, n_pq, n_qr);
      std::size_t new_tri_3_id = create_new_triangle(n_r, n_pr, n_qr);

      // change the triangle 'i' to be the middle triangle
      triangle[0] = n_pq;
      triangle[1] = n_qr;
      triangle[2] = n_pr;

      // build the new incident_triangle information
      points[n_p].incident_triangles.push_back(new_tri_1_id);
      points[n_q].incident_triangles.push_back(new_tri_2_id);
      points[n_r].incident_triangles.push_back(new_tri_3_id);

      points[n_pq].incident_triangles.push_back(i);
      points[n_pq].incident_triangles.push_back(new_tri_1_id);
      points[n_pq].incident_triangles.push_back(new_tri_2_id);

      points[n_qr].incident_triangles.push_back(i);
      points[n_qr].incident_triangles.push_back(new_tri_2_id);
      points[n_qr].incident_triangles.push_back(new_tri_3_id);

      points[n_pr].incident_triangles.push_back(i);
      points[n_pr].incident_triangles.push_back(new_tri_1_id);
      points[n_pr].incident_triangles.push_back(new_tri_3_id);

      // build the new neighboring information (border info is dealt with at point creation)
      points[n_p].neighbors.insert(n_pq);
      points[n_p].neighbors.insert(n_pr);
      points[n_q].neighbors.insert(n_pq);
      points[n_q].neighbors.insert(n_qr);
      points[n_r].neighbors.insert(n_qr);
      points[n_r].neighbors.insert(n_pr);

      points[n_pq].neighbors.insert(n_p);
      points[n_pq].neighbors.insert(n_pr);
      points[n_pq].neighbors.insert(n_qr);
      points[n_pq].neighbors.insert(n_q);

      points[n_qr].neighbors.insert(n_q);
      points[n_qr].neighbors.insert(n_pq);
      points[n_qr].neighbors.insert(n_pr);
      points[n_qr].neighbors.insert(n_r);

      points[n_pr].neighbors.insert(n_p);
      points[n_pr].neighbors.insert(n_pq);
      points[n_pr].neighbors.insert(n_qr);
      points[n_pr].neighbors.insert(n_r);
    }

    // to have a proper triangulation, we need to refine some triangles who were
    // not split but have gotten a point inserted on at least one of their edges
    std::cout << triangles_buffer.size() << " buffer size before handling border "
              << triangles.size() + triangles_buffer.size() << std::endl;

    CGAL_assertion(triangles.size() == refined_triangles.size());
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      if(refined_triangles[i]) // already refined (or new), nothing to do
        continue;

      Tri& triangle = triangles[i];

      // check if the edges were refined
      const std::size_t n_p = triangle[0];
      const std::size_t n_q = triangle[1];
      const std::size_t n_r = triangle[2];

      const std::size_t n_pq = probe_edge_midpoints_map(n_p, n_q, edge_mid_points,
                                                        false/* don't create*/);
      const std::size_t n_qr = probe_edge_midpoints_map(n_q, n_r, edge_mid_points,
                                                        false/* don't create*/);
      const std::size_t n_pr = probe_edge_midpoints_map(n_p, n_r, edge_mid_points,
                                                        false/* don't create*/);

      // Alright, now we want to subdivide the triangle properly
      // - 0 mid point --> nothing to do
      // - 1 mid point --> split in 2 triangles
      // - 2 mid points --> split in 3 triangles
      // - 3 mid points --> triforce split

      bool is_split_pq = (n_pq != static_cast<std::size_t>(-1));
      bool is_split_qr = (n_qr != static_cast<std::size_t>(-1));
      bool is_split_pr = (n_pr != static_cast<std::size_t>(-1));

      if(!is_split_pq && !is_split_qr && !is_split_pr)
        continue;

      seeds_to_rebuild.insert(points[n_p].closest_seed_id);
      seeds_to_rebuild.insert(points[n_q].closest_seed_id);
      seeds_to_rebuild.insert(points[n_r].closest_seed_id);

      if(is_split_pq) // pq
        if(is_split_qr) // pq & qr
          if(is_split_pr) // pq & qr & pr
            split_3(n_p, n_q, n_r, n_pq, n_qr, n_pr, i, edge_mid_points);
          else // pq & qr & !pr
            split_2(n_p, n_q, n_r, n_pq, n_qr, i, edge_mid_points);
        else // pq & !qr
          if(is_split_pr) // pq & !qr & pr
            split_2(n_r, n_p, n_q, n_pr, n_pq, i, edge_mid_points);
          else // pq & !qr & !pr
            split_1(n_p, n_q, n_r, n_pq, i, edge_mid_points);
      else // !pq
        if(is_split_qr) // !pq & qr
          if(is_split_pr) // !pq & qr & pr
            split_2(n_p, n_r, n_q, n_pr, n_qr, i, edge_mid_points);
          else // !pq & qr & !pr
            split_1(n_q, n_r, n_p, n_qr, i, edge_mid_points);
        else // !pq & !qr
          if(is_split_pr) // !pq & !qr & pr
            split_1(n_r, n_p, n_q, n_pr, i, edge_mid_points);
          // else // !pq & !qr & !pr --> nothing to do
    }

    // merge the buffer triangles in the main triangles
    triangles.insert(triangles.end(), triangles_buffer.begin(), triangles_buffer.end());
    triangles_buffer.clear();

    std::cout << points.size() << " points after grid refinement" << std::endl;
    std::cout << triangles.size() << " triangles after grid refinement" << std::endl;

CGAL_expensive_assertion_code(
    // some quick debug on incidency and adjacency
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& triangle = triangles[i];
      const Grid_point& gp = points[triangle[0]];
      const Grid_point& gq = points[triangle[1]];
      const Grid_point& gr = points[triangle[2]];

      CGAL_assertion(gp.neighbors.find(gq.index) != gp.neighbors.end());
      CGAL_assertion(gp.neighbors.find(gq.index) != gp.neighbors.end());
      CGAL_assertion(gq.neighbors.find(gp.index) != gq.neighbors.end());
      CGAL_assertion(gq.neighbors.find(gr.index) != gq.neighbors.end());
      CGAL_assertion(gr.neighbors.find(gq.index) != gr.neighbors.end());
      CGAL_assertion(gr.neighbors.find(gp.index) != gr.neighbors.end());

      CGAL_assertion(std::find(gp.incident_triangles.begin(),
                               gp.incident_triangles.end(), i) != gp.incident_triangles.end());
      CGAL_assertion(std::find(gq.incident_triangles.begin(),
                               gq.incident_triangles.end(), i) != gq.incident_triangles.end());
      CGAL_assertion(std::find(gr.incident_triangles.begin(),
                               gr.incident_triangles.end(), i) != gr.incident_triangles.end());
    }
)

    build_aabb_tree();

#ifdef USE_FULL_REBUILDS
    reset();
    locate_and_initialize_seeds();
#else
    refresh_grid_point_states();

//    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
//      points[i].state = KNOWN;

    // go an extra step (is it necessary) for points whose colors we are spreading:
    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      Grid_point& gp = points[i];
      if(seeds_to_rebuild.find(gp.closest_seed_id) != seeds_to_rebuild.end())
        gp.reset();
    }

//    CGAL_precondition(trial_points.empty());
//    std::cout << "Rebuild: " << std::endl;
//    boost::unordered_set<std::size_t>::iterator it = seeds_to_rebuild.begin(),
//                                                end = seeds_to_rebuild.end();
//    for(; it!=end; ++it)
//    {
//      const Point_2& p = seeds[*it];
//      locate_and_initialize(p, *it);
//    }
//    CGAL_postcondition(!trial_points.empty() && "empty prio queue ?");

    locate_and_initialize_seeds();
#endif
    clear_dual();
    spread_distances(false);
    return true;
  }

  void clear_dual()
  {
    edge_incident_tris.clear();

    edge_bisectors.clear();
    dual_edges.clear();
    dual_triangles.clear();

    size_queue.clear();
    distortion_queue.clear();
    intersection_queue.clear();
    quality_queue.clear();

    refinement_point = NULL;

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

//    std::cout << "queues: " << '\n';
//    std::cout << "size queue : " << size_queue.size() << '\n';
//    typename PQ::const_iterator pqcit = size_queue.begin();
//    for(; pqcit!=size_queue.end(); ++pqcit)
//      std::cout << *pqcit << '\n';

    if(!size_queue.empty())
       best_entry = *(size_queue.begin());
    else if(!intersection_queue.empty())
      best_entry = *(intersection_queue.begin());
    else if(!distortion_queue.empty())
      best_entry = *(distortion_queue.begin());
    else if(!quality_queue.empty())
      best_entry = *(quality_queue.begin());

    refinement_point = best_entry.get<1>();

    if(refinement_point)
    {
      std::cout << "naturally we picked : " << refinement_point->index << " [";
      std::cout << refinement_point->point << "]" << std::endl;
      std::cout << "second was: " << best_entry.get<2>() << std::endl;
    }

    if(!refinement_point)
    {
      std::cout << "Couldn't find a ref point, need more initial points" << std::endl;
      return false;
//      std::exit(EXIT_FAILURE);
    }
    else
      refine_seeds(refinement_point->point);
    return true;
  }

  void output_grid(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
//    std::ofstream out_distortion_bb((str_base + "_distortion.bb").c_str());

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';
    out << "Vertices" << '\n';
    out << points.size() << '\n';
    out_bb << "2 1 " << points.size() << " 2" << '\n';
//    out_distortion_bb << "2 1 " << points.size() << " 2" << '\n';

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << i+1 << '\n';

      out_bb << gp.distance_to_closest_seed << '\n';
//      out_bb << gp.closest_seed_id << '\n';
//      out_distortion_bb << gp.distortion_to_seed() << '\n';
//      out_bb << (gp.is_Voronoi_vertex?"1":"0") << '\n';
    }

    out << "Triangles" << '\n';
    out << triangles.size() << '\n';
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      std::set<std::size_t> materials;
      for(int j=0; j<3; ++j)
      {
        materials.insert(points[triangles[i][j]].closest_seed_id);
        out << triangles[i][j]+1 << " ";
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size()+1);
      out << mat << '\n';
    }

//    out_distortion_bb << "End" << '\n';
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
      if(it->get<1>()->distance_to_closest_seed < entry.get<1>()->distance_to_closest_seed)
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

  void test_simplex(const Grid_point* gp,
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
      FT gpd = gp->distance_to_closest_seed;
      if(gpd > max_size)
      {
        PQ_entry pqe = boost::make_tuple(dual_simplex, gp, gpd);
        insert_in_PQ(pqe, size_queue, 0);
        size_queue.insert(pqe);
        return;
      }
    }

    // distortion
    if(max_distortion > 1.)
    {
      FT gamma = gp->distortion_to_seed();
      if(gamma > max_distortion)
      {
        PQ_entry pqe = boost::make_tuple(dual_simplex, gp, gamma);
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
          PQ_entry pqe = boost::make_tuple(dual_simplex, gp, area);
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
          PQ_entry pqe = boost::make_tuple(dual_simplex, gp, 1./qual);
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
     test_simplex(&(points[*it]), dual_simplex);
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

  // Optimization --------------------------------------------------------------

  FT compute_tangent_angle(const Grid_point& gp, const Point_2& seed)
  {
    // to approximate the tangent of a geodesic, we take the vector [seed-gq]
    // where gq is a point on the geodesic path from seed to gp

    if(gp.point == seed)
      return 0;

    std::deque<const Grid_point*> geodesic_path; // starting at the grid point closest to seed
    geodesic_path.push_front(&gp);

    std::size_t anc_id = gp.ancestor;
    while(anc_id != static_cast<std::size_t>(-1))
    {
      const Grid_point& anc = points[anc_id];
      geodesic_path.push_front(&anc);
      anc_id = anc.ancestor;
    }

    // determine gq
    // Since we have cleared out the case gp = seed, we know that :
    // - geodesic_path.size() > 0
    // - even if we take gq_index = 0, we have a non degenerate vector
//    std::size_t gq_index = geodesic_path.size()-1; // TMP
    std::size_t gq_index = (std::min)(static_cast<std::size_t>(1e10),
                                      geodesic_path.size()-1);

    CGAL_assertion(gq_index < geodesic_path.size());
    const Grid_point* gq = geodesic_path[gq_index];

//    std::cout << "---" << std::endl;
//    std::cout << "gq_index: " << gq_index << " (total path length: " << geodesic_path.size() << ")" << std::endl;
//    std::cout << "dir vector is : " << gq->point.x() - seed.x() << " " << gq->point.y() - seed.y() << std::endl;
//    std::cout << "headed towards: " << seed << std::endl;

    // compute the angle with the vector (0,1) starting at seed
    FT theta = std::atan2(gq->point.y() - seed.y(), gq->point.x() - seed.x());

    return theta;
  }

  Point_2 map_point_to_tangent_space(const Grid_point& gp)
  {
    const FT r = gp.distance_to_closest_seed;
    const std::size_t seed_id = gp.closest_seed_id;
    const Point_2& seed = seeds[seed_id];

    // compute the tangent and the angle with the horizontal axis
    FT theta = compute_tangent_angle(gp, seed);

    // we don't actually care much for polar coordinates, so we keep them as
    // Cartesian coordinates instead
    FT x = r * std::cos(theta);
    FT y = r * std::sin(theta);

//    std::cout << " mapped: " << gp.point.x() << " " << gp.point.y() << " to "
//              << x << " " << y << " r/theta: " << r << " " << theta
//              << " (seed_id : " << seed_id << ")" << std::endl;

     return Point_2(x, y);
  }

  void map_grid_points_to_tangent_spaces(std::vector<Point_2>& mapped_points)
  {
    std::cout << "build map grids..." << std::endl;

    // transform grid points from the manifold to the 'tangent' space of their
    // closest seed with the log map.

    // The new coordinates are polar geodesic coordinates
    // - The radius is simply the geodesic distance
    // - The angle is given by the angle of the tangent of the geodesic at the seed

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      const Grid_point& gp = points[i];
      mapped_points[i] = map_point_to_tangent_space(gp);

#ifndef COMPUTE_PRECISE_VOR_VERTICES
    // If the grid point is also a Voronoi vertex, we keep it in a special memory
      if(gp.is_Voronoi_vertex)
        Voronoi_vertices[gp.closest_seed_id].push_back(gp);
#endif
    }
  }

  void map_precise_Voronoi_vertices_to_tangent_spaces(std::vector<Point_2>& mapped_points)
  {
    CGAL_assertion(Voronoi_vertices.size() == seeds.size());
    std::cout << "map precise" << std::endl;

    // count the number of precise Voronoi vertices
    std::size_t precise_Vor_vertices_count = 0;
    for(std::size_t seed_id=0, vvs=Voronoi_vertices.size(); seed_id<vvs; ++seed_id)
      precise_Vor_vertices_count += Voronoi_vertices[seed_id].size();

    // resize the mapped points to add the precise Vor vertices
    mapped_points.resize(points.size() + precise_Vor_vertices_count);

    // to number the precise Vor vertices (they are not base mesh points)
    std::size_t precise_Vor_vertex_index = points.size();

    for(std::size_t seed_id=0, vvs=Voronoi_vertices.size(); seed_id<vvs; ++seed_id)
    {
      Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];
      for(std::size_t i=0, vvs=Vor_vertices.size(); i<vvs; ++i)
      {
        // gp is NOT a point of the base mesh
        Grid_point& gp = Vor_vertices[i];
        gp.index = precise_Vor_vertex_index++;

        CGAL_assertion(gp.closest_seed_id == seed_id);

        const Point_2& mapped_p = map_point_to_tangent_space(gp);
        mapped_points[gp.index] = mapped_p;
      }
    }
  }

  Point_2 compute_centroid_as_closest_grid_point(const std::size_t seed_id,
                                                 const Point_2& mapped_centroid,
                                                 const std::vector<Point_2>& mapped_points)
  {
    // find the closest mapped point and return its corresponding unmapped grid point

    std::size_t closest_mapped_grid_point_id;
    FT min_sq_dist = FT_inf;
    K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();

    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      // only look for the closest points with the same seed id. (Is this correct ?)
      if(points[i].closest_seed_id != seed_id)
        continue;

      // to constrain to a border: check if the init point on the grid is a border
      // and ignore points here that are not on the border ? fixme

      // ugly hack for now (I know the border seeds all have an id lower than n)
//      if(seed_id < 8 && !points[i].is_on_domain_border)
//        continue;

      const Point_2& test_p = mapped_points[i];
      FT sq_d = sqd(mapped_centroid, test_p);

      if(sq_d < min_sq_dist)
      {
        closest_mapped_grid_point_id = i;
        min_sq_dist = sq_d;
      }
    }

    std::cout << "closest mapped point found: " << closest_mapped_grid_point_id << std::endl;
    return points[closest_mapped_grid_point_id].point;
  }

  Point_2 compute_centroid_with_barycentric_info(const std::size_t seed_id,
                                                 const Point_2& mapped_centroid,
                                                 const std::vector<Point_2>& mapped_points)
  {
    // compute the barycentric coordinates of the mapped_centroid in a mapped triangle
    // and then compute the new seed through the coordinates in the unmapped triangle

    FT unmapped_centroid_x = 0., unmapped_centroid_y = 0.;
    bool found = false;
    FT lambda_p = 0., lambda_q = 0., lambda_r = 0.;
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& tr = triangles[i];
      const Grid_point& gp = points[tr[0]];
      const Grid_point& gq = points[tr[1]];
      const Grid_point& gr = points[tr[2]];

      if(gp.closest_seed_id != seed_id ||
         gq.closest_seed_id != seed_id ||
         gr.closest_seed_id != seed_id)
        continue;

      const Point_2& p = mapped_points[tr[0]];
      const Point_2& q = mapped_points[tr[1]];
      const Point_2& r = mapped_points[tr[2]];
/*
      std::cout << "triangle (bar) : " << std::endl;
      std::cout << gp.index << " (" << p << ") seed : " << gp.closest_seed_id << std::endl;
      std::cout << gq.index << " (" << q << ") seed : " << gq.closest_seed_id << std::endl;
      std::cout << gr.index << " (" << r << ") seed : " << gr.closest_seed_id << std::endl;
*/
      compute_bary_weights(mapped_centroid, p, q, r, lambda_p, lambda_q, lambda_r);

      if(lambda_p >= 0. && lambda_q >= 0. && lambda_r >= 0.)
      {
        unmapped_centroid_x = lambda_p*gp.point.x() + lambda_q*gq.point.x() + lambda_r*gr.point.x();
        unmapped_centroid_y = lambda_p*gp.point.y() + lambda_q*gq.point.y() + lambda_r*gr.point.y();
        found = true;
        break;
      }
    }
    CGAL_assertion(found);

    return Point_2(unmapped_centroid_x, unmapped_centroid_y);
  }

  void add_to_centroids(const Grid_point& a, const Grid_point& b, const Grid_point& c,
                        std::vector<std::vector<std::pair<Point_2, FT> > >& centroids)
  {
    // this function inserts the centroid of abc in the centroid map
    // FOR THE SEED CORRESPONDING TO 'a' !!

    Point_2 centroid = CGAL::centroid(a.point, b.point, c.point);
    FT area = triangle_area_in_metric(a.point, b.point, c.point);
    CGAL_postcondition(area != 0.);

    std::size_t seed_id = a.closest_seed_id;
    CGAL_precondition(seed_id < seeds.size());
    centroids[seed_id].push_back(std::make_pair(centroid, area));
  }

  void compute_all_centroids(std::vector<std::vector<std::pair<Point_2, FT> > >& centroids)
  {
    std::size_t real_points_n = points.size();
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& triangle = triangles[i];

      Simplex colors;
      for(int i=0; i<3; ++i)
        colors.insert(points[triangle[i]].closest_seed_id);

      if(colors.size() == 1)
      {
        add_to_centroids(points[triangle[0]], points[triangle[1]], points[triangle[2]],
                         centroids);
      }
      else if(colors.size() == 2)
      {
        // compute the midpoint on the edges who have different colors
        int pos = 0;
        int third_edge_first_point_id = -1; // the edge for which both extremities have the same color
        boost::array<Grid_point, 2> mid_pts;
        for(std::size_t i=0; i<3; ++i)
        {
          CGAL_assertion(triangle[i] < points.size());
          const Grid_point& gp0 = points[triangle[i]];
          const Grid_point& gp1 = points[triangle[(i+1)%3]];

          if(gp0.closest_seed_id == gp1.closest_seed_id)
          {
            third_edge_first_point_id = i;
            continue;
          }

          CGAL_assertion(pos == 0 || pos == 1);
          mid_pts[pos++] = Vor_vertex_on_edge(gp0.index, gp1.index);
        }

        CGAL_postcondition(third_edge_first_point_id >= 0 &&
                           third_edge_first_point_id < 3);

        // 'triangle' is split in 3 smaller triangles.
        // a is the one with the different color
        // a-i1-i2 associated to the seed that 'a' is linked to
        // b-i2-i1 and c-b-i2 are associated to the seed that 'b' and 'c' are linked to

        //            a
        //           / |
        //          /  |
        //         /   |
        //       i0----i1
        //       /    /|
        //      /   /  |
        //     /  /    |
        //    / /      |
        // b //________| c

        const Grid_point& b = points[triangle[third_edge_first_point_id]];
        const Grid_point& c = points[triangle[(third_edge_first_point_id + 1)%3]];
        const Grid_point& a = points[triangle[(third_edge_first_point_id + 2)%3]];

        // Which of mid_pt_1 and mid_pt_2 are i1 and i2 depends on where the edge
        // with same colors is...
        // for example, if i = 1, we first looked at ab so mid_pt_1 is on [ab]
        const Grid_point& i0 = (third_edge_first_point_id == 1) ? mid_pts[0] : mid_pts[1];
        const Grid_point& i1 = (third_edge_first_point_id == 1) ? mid_pts[1] : mid_pts[0];

        add_to_centroids(a, i0, i1, centroids);
        add_to_centroids(b, i1, i0, centroids);
        add_to_centroids(c, i1, b, centroids);
      }
      else // colors.size() == 3
      {
        // compute the on the edges
        boost::array<Grid_point,3> mid_pts;
        for(int i=0; i<3; ++i)
          mid_pts[i] = Vor_vertex_on_edge(triangle[i], triangle[(i+1)%3]);

        // compute roughly the center point
        const Grid_point& v = compute_precise_Voronoi_vertex(triangle);

        // 'triangle' is split in 6 smaller triangles.
        // a-i1-v & a-i3-v --> associated to the seed a->closest_seed_id
        // b-i1-v & b-i2-v --> associated to the seed b->closest_seed_id
        // c-i2-v & c-i3-v --> associated to the seed c->closest_seed_id

        //            a
        //           / |
        //          /  |
        //         /   |
        //       i1    |
        //       /|    |
        //      /  |   |
        //     /   v---i3
        //    /    |   |
        // b /____i2___| c

        add_to_centroids(points[triangle[0]], v, mid_pts[0], centroids);
        add_to_centroids(points[triangle[0]], v, mid_pts[2], centroids);

        add_to_centroids(points[triangle[1]], v, mid_pts[0], centroids);
        add_to_centroids(points[triangle[1]], v, mid_pts[1], centroids);

        add_to_centroids(points[triangle[2]], v, mid_pts[1], centroids);
        add_to_centroids(points[triangle[2]], v, mid_pts[2], centroids);
      }
    }
    shave_off_virtual_points(real_points_n);
  }

  Point_2 compute_centroid_with_grid_triangles_precomputed(const std::size_t seed_id,
                         std::vector<std::vector<std::pair<Point_2, FT> > > centroids)
  {
    CGAL_precondition(seed_id < centroids.size());
    const std::vector<std::pair<Point_2, FT> >& seed_centroids = centroids[seed_id];

    FT total_area = 0.;
    FT centroid_x = 0., centroid_y = 0.;
    for(std::size_t i=0, scs=seed_centroids.size(); i<scs; ++i)
    {
      const std::pair<Point_2, FT>& centroid = seed_centroids[i];
      const Point_2& c = centroid.first;
      FT area = centroid.second;

      total_area += area;
      centroid_x += area * c.x();
      centroid_y += area * c.y();
    }

    CGAL_assertion(total_area != 0.);
    centroid_x /= total_area;
    centroid_y /= total_area;

    return Point_2(centroid_x, centroid_y);
  }

  Point_2 compute_centroid_with_grid_triangles(const std::size_t seed_id)
  {
    // compute the centroid through the sum c = sum_i (ci*area_i) / sum_i (area_i)
    // with the tiny triangles of the base mesh
    // it's not exact because we only consider the triangles with all vertices
    // having seed_id as closest seed (thus ignoring border triangles that
    // only have <= 2 vertices that have seed_id as closest_seed)

    FT total_area = 0.;
    FT centroid_x = 0., centroid_y = 0.;

    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& tr = triangles[i];
      const Grid_point& gp = points[tr[0]];
      const Grid_point& gq = points[tr[1]];
      const Grid_point& gr = points[tr[2]];

      if(gp.closest_seed_id != seed_id ||
         gq.closest_seed_id != seed_id ||
         gr.closest_seed_id != seed_id)
        continue;

//      FT area = triangle_area(gp.point, gq.point, gr.point);
      FT area = triangle_area_in_metric(gp.point, gq.point, gr.point);

      total_area += area;
      Point_2 local_centroid = CGAL::centroid(gp.point, gq.point, gr.point);

/*
      std::cout << "triangle (grid) : " << std::endl;
      std::cout << gp.index << " (" << gp.point << ")" << std::endl;
      std::cout << gq.index << " (" << gq.point << ")" << std::endl;
      std::cout << gr.index << " (" << gr.point << ")" << std::endl;
      std::cout << "area: " << area << std::endl;
      std::cout << "local centroid: " << local_centroid << std::endl;
*/
      centroid_x += area * local_centroid.x();
      centroid_y += area * local_centroid.y();
    }
//    std::cout << "total area: " << total_area << std::endl;

    CGAL_assertion(total_area != 0.);
    centroid_x /= total_area;
    centroid_y /= total_area;

    return Point_2(centroid_x, centroid_y);
  }

  Point_2 compute_centroid_with_mapped_grid_triangles(const std::size_t seed_id,
                                                      const std::vector<Point_2>& mapped_points)
  {
    // decompose the mapped Voronoi cells in mapped triangles and compute
    // the centroid on these mapped triangles as we would for Euclidean triangles

    FT total_area = 0.;
    FT centroid_x = 0., centroid_y = 0.;

    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& tr = triangles[i];
      const Grid_point& gp = points[tr[0]];
      const Grid_point& gq = points[tr[1]];
      const Grid_point& gr = points[tr[2]];

      if(gp.closest_seed_id != seed_id ||
         gq.closest_seed_id != seed_id ||
         gr.closest_seed_id != seed_id)
        continue;

      const Point_2& p = mapped_points[tr[0]];
      const Point_2& q = mapped_points[tr[1]];
      const Point_2& r = mapped_points[tr[2]];


      FT area = triangle_area(p, q, r);
      total_area += area;
      Point_2 local_centroid = CGAL::centroid(p, q, r);

/*
      std::cout << "triangle (mapped grid) : " << std::endl;
      std::cout << gp.index << " (" << p << ")" << std::endl;
      std::cout << gq.index << " (" << q << ")" << std::endl;
      std::cout << gr.index << " (" << r << ")" << std::endl;
      std::cout << "area: " << area << std::endl;
      std::cout << "local centroid: " << local_centroid << std::endl;
*/
      centroid_x += area * local_centroid.x();
      centroid_y += area * local_centroid.y();
    }
    std::cout << "total area: " << total_area << std::endl;

    CGAL_assertion(total_area != 0.);
    centroid_x /= total_area;
    centroid_y /= total_area;

    return Point_2(centroid_x, centroid_y);
  }

  Point_2 compute_centroid_with_voronoi_vertices(const std::size_t seed_id)
  {
    // in isotropic, the bisectors are segments between the voronoi vertices,
    // so we just decompose the full Voronoi cells in triangles of the form :
    // [seed, vor_vert_1, Vor_vert_2]

    const Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];
    if(Vor_vertices.size() < 3)
    {
      std::cout << "WARNING: NOT ENOUGH VORONOI VERTICES TO OPTIMIZE SEED: " << seed_id << std::endl;
      return Point_2(0.,0.);
    }

    const Point_2 old_seed = seeds[seed_id];
    FT total_area = 0;
    FT centroid_x = 0., centroid_y = 0.;

    Voronoi_vertices_container::const_iterator it = Vor_vertices.begin();
    Voronoi_vertices_container::const_iterator next_it = ++(Vor_vertices.begin());
    Voronoi_vertices_container::const_iterator end = Vor_vertices.end();
    for(; it!=end; ++it, ++next_it)
    {
      if(next_it == end)
        next_it = Vor_vertices.begin();

      const Point_2& pi = it->point;
      const Point_2& pj = next_it->point;

//      FT area = triangle_area(old_seed, pi, pj);
      FT area = triangle_area_in_metric(old_seed, pi, pj);
      total_area += area;

      Point_2 local_centroid = CGAL::centroid(old_seed, pi, pj);

/*
      std::cout << "triangle (voronoi vert) : " << std::endl;
      std::cout << seed_id << " (" << old_seed << ")" << std::endl;
      std::cout << *it << " (" << pi << ")" << std::endl;
      std::cout << *next_it << " (" << pj << ")" << std::endl;
      std::cout << "area: " << area << std::endl;
      std::cout << "local centroid: " << local_centroid << std::endl;
*/

      centroid_x += area * local_centroid.x();
      centroid_y += area * local_centroid.y();
    }

    CGAL_assertion(total_area != 0.);
    centroid_x /= total_area;
    centroid_y /= total_area;

    return Point_2(centroid_x, centroid_y);
  }

  Point_2 compute_mapped_centroid_with_mapped_Vor_vertices(const std::size_t seed_id,
                                                            const std::vector<Point_2>& mapped_points)
  {
    const Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];

    std::cout << "optimize seed " << seed_id << std::endl;
    std::cout << "currently at position: " << seeds[seed_id] << std::endl;
    std::cout << Vor_vertices.size() << " vor vertices" << std::endl;

    if(Vor_vertices.size() < 3)
    {
      std::cout << "WARNING: NOT ENOUGH VORONOI VERTICES TO OPTIMIZE SEED: " << seed_id << std::endl;
      return Point_2(0.,0.);
    }

    // compute the centroid of the polygon composed by the mapped Voronoi vertices
    // see wikipedia for the formula...
    FT area = 0;
    FT centroid_x = 0., centroid_y = 0.;

    Voronoi_vertices_container::const_iterator it = Vor_vertices.begin();
    Voronoi_vertices_container::const_iterator next_it = ++(Vor_vertices.begin());
    Voronoi_vertices_container::const_iterator end = Vor_vertices.end();
    for(; it!=end; ++it, ++next_it)
    {
      if(next_it == end)
        next_it = Vor_vertices.begin();

      const Point_2& pi = mapped_points[it->index];
      const Point_2& pj = mapped_points[next_it->index];

      std::cout << "current Vor vertex: " << it->index << " " << pi;
      std::cout << " which is an angle of: " << std::atan2(mapped_points[it->index].y(),
                                                           mapped_points[it->index].x()) << std::endl;

      area += pi.x()*pj.y() - pj.x()*pi.y();
      centroid_x += (pi.x() + pj.x()) * (pi.x()*pj.y() - pj.x()*pi.y());
      centroid_y += (pi.y() + pj.y()) * (pi.x()*pj.y() - pj.x()*pi.y());
    }
    area *= 0.5;
    CGAL_assertion(area != 0);

    FT denom = 1./(6.*area);

    centroid_x *= denom;
    centroid_y *= denom;
    return Point_2 (centroid_x, centroid_y);
  }

  FT optimize_seed(const std::size_t seed_id,
                   const std::vector<Point_2>& mapped_points,
                   const std::vector<std::vector<std::pair<Point_2, FT> > >& centroids,
                   const int counter)
  {
    const Point_2 old_seed = seeds[seed_id];
//    Point_2 mapped_centroid = compute_mapped_centroid_with_mapped_Vor_vertices(seed_id, mapped_points);
//    std::cout << "(mapped) centroid coordinates from mapped Vor vertices " << mapped_centroid  << std::endl;

    // -------------------------------------------------------------------------
    std::cout << "UNMAPPING POSSIBILITIES :" << std::endl;

//    Point_2 closest_grid_point_centroid = compute_centroid_as_closest_grid_point(seed_id, mapped_centroid, mapped_points);
//    std::cout << "unmapped centroid as closest unmapped grid point " << closest_grid_point_centroid << std::endl;

    // disabled till degenerate mapped triangles are fixed fixme
//    Point_2 bar_centroid = compute_centroid_with_barycentric_info(seed_id, mapped_centroid, mapped_points);
//    std::cout << "unmapped centroid from mapped polygon (barycentric) " << bar_centroid << std::endl;

    // --
      std::cout << "FOR REFERENCE, WITHOUT TRANSFORMATION GIVES :" << std::endl;

      Point_2 alt_centroid_1 = compute_centroid_with_grid_triangles(seed_id);
      std::cout << "centroid with grid triangles " << alt_centroid_1 << std::endl;

//      Point_2 alt_centroid_2 = compute_centroid_with_grid_triangles_precomputed(seed_id,
//                                                                                centroids);
//      std::cout << "centroid with grid triangles " << alt_centroid_2 << std::endl;

      // below is not a good idea if the metric field is not uniform
//      Point_2 alt_centroid_3 = compute_centroid_with_voronoi_vertices(seed_id);
//      std::cout << "centroid with voronoi vertices " << alt_centroid_3 << std::endl;
    // --

    std::cout << "ALTERNATIVE MAPPED CENTROID COMPUTATION : " << std::endl;

//    Point_2 mapped_grid_centroid = compute_centroid_with_mapped_grid_triangles(seed_id, mapped_points);
//    std::cout << "(mapped) centroid with mapped grid triangles " << mapped_grid_centroid << std::endl;

//    Point_2 closest_grid_point_centroid_alternate =
//        compute_centroid_as_closest_grid_point(seed_id, mapped_grid_centroid, mapped_points);
//    std::cout << "unmapped centroid as closest unmapped grid point (alternate)"
//              << closest_grid_point_centroid_alternate << std::endl;

    // disabled till degenerate mapped triangles is fixed fixme
//    Point_2 bar_centroid_2 = compute_centroid_with_barycentric_info(seed_id, mapped_grid_centroid, mapped_points);
//    std::cout << "unmapped centroid from mapped grid centroid (barycentric) " << bar_centroid_2 << std::endl;

    // -------------------------------------------------------------------------

    // output the mapping
//    output_mapping(seed_id, mapped_grid_centroid, mapped_points, counter);

    // need seed is now the closest mapped point, unmapped back to the manifold
    const Point_2& new_seed = alt_centroid_1;
    seeds[seed_id] = new_seed;
    seeds_m[seed_id] = mf->compute_metric(new_seed);

    // squared displacement between the old and the new seed
    K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();
    FT displacement = sqd(old_seed, new_seed);
    std::cout << "seed: " << seed_id << " displacement: " << displacement << std::endl;

    return displacement;
  }

  // we want to have the Voronoi vertices ordered in a cycle
  template<typename GP_iter>
  Point_2 barycenter(GP_iter it, GP_iter end,
                     const std::vector<Point_2>& mapped_points)
  {
    std::size_t n = end - it;
    FT x = 0., y =0.;
    for(; it!=end; ++it)
    {
      x += mapped_points[it->index].x();
      y += mapped_points[it->index].y();
    }
    CGAL_assertion(n != 0);
    x /= n;
    y /= n;

    return Point_2(x, y);
  }

  struct Voronoi_vertices_comparer
  {
    const std::vector<Point_2>& mapped_points;
    const Point_2& pol_bar;

    bool operator()(const Grid_point& left, const Grid_point& right)
    {
      const Point_2& left_p = mapped_points[left.index];
      const Point_2& right_p = mapped_points[right.index];

      FT left_theta = std::atan2(left_p.y() - pol_bar.y(), left_p.x() - pol_bar.x());
      FT right_theta = std::atan2(right_p.y() - pol_bar.y(), right_p.x() - pol_bar.x());

      return left_theta < right_theta;
    }

    Voronoi_vertices_comparer(const std::vector<Point_2>& mapped_points_,
                              const Point_2& polygon_barycenter_)
      :
        mapped_points(mapped_points_),
        pol_bar(polygon_barycenter_)
    { }
  };

  void sort_Voronoi_vertices(const std::vector<Point_2>& mapped_points)
  {
    CGAL_assertion(Voronoi_vertices.size() == seeds.size());

    for(std::size_t seed_id=0; seed_id<seeds.size(); ++seed_id)
    {
      // we want to sort the Voronoi vertices in a cycle around the seed
      Voronoi_vertices_container& seed_Vor_vertices = Voronoi_vertices[seed_id];
      Point_2 mapped_bar = barycenter(seed_Vor_vertices.begin(), seed_Vor_vertices.end(),
                                      mapped_points);

      std::sort(seed_Vor_vertices.begin(), seed_Vor_vertices.end(),
                Voronoi_vertices_comparer(mapped_points, mapped_bar));

      // just as debug :
      std::size_t old_size = seed_Vor_vertices.size();
      seed_Vor_vertices.erase(std::unique(seed_Vor_vertices.begin(),
                                           seed_Vor_vertices.end() ),
                              seed_Vor_vertices.end() );
      CGAL_assertion(seed_Vor_vertices.size() == old_size);
    }
  }

  FT compute_CVT_energy() const
  {
    FT e = 0.;
    FT third = 1./3.;

    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      const Tri& tr = triangles[i];
      const Grid_point& gp = points[tr[0]];
      const Grid_point& gq = points[tr[1]];
      const Grid_point& gr = points[tr[2]];

      FT dist = third * (gp.closest_seed_id + gq.closest_seed_id + gr.closest_seed_id);
      FT area = triangle_area_in_metric(gp.point, gq.point, gr.point);
      e += area * dist * dist;
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
      // the mapped Voronoi vertices (each point mapped according to its closest seed)
      std::vector<Point_2 > mapped_points(points.size());

      map_grid_points_to_tangent_spaces(mapped_points);
#ifdef COMPUTE_PRECISE_VOR_VERTICES
      map_precise_Voronoi_vertices_to_tangent_spaces(mapped_points);
#endif

      std::vector<std::vector<std::pair<Point_2, FT> > > centroids(seeds.size());
      compute_all_centroids(centroids);

      // must sort the Voronoi vertices so that they are ordered in a cycle
      // to be able to use the polygon centroid's formula
      sort_Voronoi_vertices(mapped_points);

      // Optimize each seed
      FT cumulated_displacement = 0;
      for(std::size_t i=0, ss=seeds.size(); i<ss; ++i)
        cumulated_displacement += optimize_seed(i, mapped_points,
                                                centroids, counter);

      FT e = compute_CVT_energy();

      std::cout << "at : " << counter << ", cumulated displacement : "
                << cumulated_displacement << std::endl;
      std::cout << "energy at : " << e << std::endl;

      reset();
      locate_and_initialize_seeds();
      spread_distances(true/*use_dual_shenanigans*/);

      std::ostringstream opti_out;
      opti_out << "optimized_" << str_base_mesh << "_tr_" << counter << std::ends;
      output_grid_data_and_dual(opti_out.str().c_str());

      is_optimized = (++counter > max_opti_n); // (cumulated_displacement < sq_bbox_diag_l * 1e-5);
      if(cumulated_displacement == 0.0)
        break;
    }
    while(!is_optimized);
  }

  // output stuff --------------------------------------------------------------
  void add_simplex_to_triangulation(const Grid_point* gp,
                                    const Simplex& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    std::size_t length = gp->ancestor_path_length;
    if(length < min_ancestor_path_length)
      std::cout << "the canvas is thin (" << length << ") at : " << gp->index << " ("
                << gp->point << ") dual simplex of size: " << dual_simplex.size()
                << " and cellid: " << gp->closest_seed_id << std::endl;

    // test against criteria should be moved somewhere else than during the dual computations todo
    if(dual_simplex.size() == 3 ||
       (dual_simplex.size() == 2 && gp->is_on_domain_border))
      test_simplex(gp, dual_simplex);

    add_simplex_to_triangulation(dual_simplex);
  }

  void add_simplex_to_triangulation(const Tri& grid_tri,
                                    const std::set<std::size_t>& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    int max_dist_p = -1;
    FT max = -FT_inf;
    for(int i=0; i<3; ++i)
    {
      const Grid_point& gp = points[grid_tri[i]];
      const FT d = gp.distance_to_closest_seed;

      // fixme we want to always take the pt on the border, it might not always
      // be the case
//      if(gp.is_on_domain_border)
//        std::cout << gp.index << " on border in add_simplex" << std::endl;

      if(d > max)
      {
        max = d;
        max_dist_p = grid_tri[i];
      }
    }

//    if(points[grid_tri[0]].is_on_domain_border ||
//       points[grid_tri[1]].is_on_domain_border ||
//       points[grid_tri[2]].is_on_domain_border)
//      std::cout << "picked " << points[max_dist_p].index << std::endl;

    // keep as dual point the farthest grid point -- at least until we have
    // precise Voronoi computations)
    return add_simplex_to_triangulation(&(points[max_dist_p]), dual_simplex);
  }

  void check_edelsbrunner()
  {
    std::cout << "it's Edel's time" << std::endl;

    // for each simplex, associate all the triangles (vector of size_t, the size_t
    // is the index of the triangle in 'triangles') of the base grid that correspond
    // to that simplex...
    // Then check that this tet set satisfies the closed ball property
    typedef std::map<Simplex, std::vector<std::size_t> > Dual_map;

    Dual_map dual_map;
    int failed_blob_counter = 0;

    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      Simplex dual_simplex;
      const Tri& triangle = triangles[i];

      for(int j=0; j<3; ++j)
        dual_simplex.insert(points[triangle[j]].closest_seed_id);

      std::pair<typename Dual_map::iterator, bool> is_insert_successful;
      std::vector<std::size_t> vec;
      vec.push_back(i);
      is_insert_successful = dual_map.insert(std::make_pair(dual_simplex, vec));

      if(!is_insert_successful.second)
        (is_insert_successful.first)->second.push_back(i);
    }

    // now, for each simplex, we want the triangles corresponding to that simplex
    // to be a single blob

    Dual_map::const_iterator dmit = dual_map.begin();
    Dual_map::const_iterator dmend = dual_map.end();
    for(; dmit!=dmend; ++dmit)
    {
/*
      const Simplex& simplex = dmit->first;
      std::cout << "looking at the simplex : ";
      typename Simplex::iterator sit = simplex.begin();
      typename Simplex::iterator send = simplex.end();
      for(; sit!=send; ++sit)
        std::cout << *sit << " ";
      std::cout << std::endl;
*/
      const std::vector<std::size_t>& grid_triangles = dmit->second;

      // collect all the edges that appear once
      std::set<Simplex> edges;

      for(std::size_t i=0, gts=grid_triangles.size(); i<gts; ++i)
      {
        const Tri& tr = triangles[grid_triangles[i]];
        for(int j=0; j<3; ++j)
        {
          Simplex edge;
          edge.insert(tr[j]); edge.insert(tr[(j+1)%3]);

          std::pair<std::set<Simplex>::iterator, bool> is_insert_successful;
          is_insert_successful = edges.insert(edge);

          if(!is_insert_successful.second) // already in 'edges', edge is internal, ignore it
            edges.erase(is_insert_successful.first);
        }
      }

//      std::cout << edges.size() << " edges" << std::endl;
      CGAL_assertion(!edges.empty());

      // check that all the border edges of the blob actually make a blob
      // note: 'bow-ties' configurations could appear... what to do ?
      std::vector<bool> visited_status(points.size(), true);
      std::map<std::size_t, std::deque<std::size_t> > neighbors;

      std::set<Simplex>::iterator it = edges.begin();
      std::set<Simplex>::iterator end = edges.end();
      for(; it!=end; ++it)
      {
        const Simplex& edge = *it;
        const std::size_t e1 = *(edge.begin());
        const std::size_t e2 = *++(edge.begin());

        std::pair<std::map<std::size_t, std::deque<std::size_t> >::iterator, bool>
            is_insert_successful;
        std::deque<std::size_t> vec;

        vec.push_back(e2);
        is_insert_successful = neighbors.insert(std::make_pair(e1, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(e2);

        vec.clear();
        vec.push_back(e1);
        is_insert_successful = neighbors.insert(std::make_pair(e2, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(e1);

        visited_status[e1] = false;
        visited_status[e2] = false;
      }

//      std::cout << "built neighbors " << neighbors.size() << std::endl;
      CGAL_assertion(!neighbors.empty());

      std::deque<std::size_t> points_to_visit;
      points_to_visit.push_back((neighbors.begin())->first);
      while(!points_to_visit.empty())
      {
//        std::cout << "pts to visit: " << points_to_visit.size() << std::endl;
        const std::size_t current_point = points_to_visit.front();
        points_to_visit.pop_front();

        visited_status[current_point] = true;

//        std::cout << current_point << " has ";
//        std::cout << neighbors[current_point].size() << " neighbors" << std::endl;

        std::deque<std::size_t>::iterator dit = neighbors[current_point].begin();
        std::deque<std::size_t>::iterator dend = neighbors[current_point].end();
        for(; dit!=dend; ++dit)
        {
          if(!visited_status[*dit])
            points_to_visit.push_back(*dit);
        }
      }

      for(std::size_t i=0, vss=visited_status.size(); i<vss; ++i)
      {
        if(!visited_status[i])
        {
//          std::cout << "Blob not achieved" << std::endl;
          failed_blob_counter++;
          break;
        }
      }
    }
    std::cout << "end edelsbrunner" << std::endl;
    std::cout << failed_blob_counter << " out of " << dual_map.size() << " failed" << std::endl;
  }

  void compute_dual(const bool are_Voronoi_vertices_needed = true)
  {
    // todo, we don't always need to compute the voronoi cell vertices
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

#ifdef REFINE_GRID
    refine_grid();
#endif

    for(std::size_t i=0, ts=triangles.size(); i < ts; ++i)
    {
      const Tri& triangle = triangles[i];

#ifndef COMPUTE_DUAL_FOR_ALL_DIMENSIONS
      // filter triangles with only one color
      if(points[triangle[0]].closest_seed_id == points[triangle[1]].closest_seed_id &&
         points[triangle[1]].closest_seed_id == points[triangle[2]].closest_seed_id)
        continue;

      if(points[triangle[0]].closest_seed_id != points[triangle[1]].closest_seed_id &&
          points[triangle[1]].closest_seed_id != points[triangle[2]].closest_seed_id &&
          points[triangle[0]].closest_seed_id != points[triangle[2]].closest_seed_id)
      { /* filter the case: 3 different colors */ }
      else
      {
      // filter triangles with _exactly_ two colors and none of the vertices
      // are on the border of the domain
        if(!points[triangle[0]].is_on_domain_border &&
           !points[triangle[1]].is_on_domain_border &&
           !points[triangle[2]].is_on_domain_border)
          continue;
      }
#endif
      Simplex dual_simplex;
      for(int j=0; j<3; ++j)
        dual_simplex.insert(points[triangle[j]].closest_seed_id);

      add_simplex_to_triangulation(triangle, dual_simplex);

      if(are_Voronoi_vertices_needed)
      {
#ifdef COMPUTE_PRECISE_VOR_VERTICES
        compute_precise_Voronoi_vertex(triangle, dual_simplex);
#else
        mark_voronoi_vertices(triangle, dual_simplex);
#endif
      }
    }

    if(are_Voronoi_vertices_needed)
    {
      // we have only computed the voronoi vertices corresponding to the dual of 3-simplex
      // to get the full cell we need to know the intersection of the voronoi cell
      // with the border of the domain
#ifdef COMPUTE_PRECISE_VOR_VERTICES
      compute_precise_Voronoi_vertices_on_border();
#else
      mark_voronoi_vertices_on_border();
#endif
    }
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

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';
    out << "Vertices" << '\n';
    out << seeds.size() << '\n';
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << i+1 << '\n';

    out << "Triangles" << '\n';
    out << dual_triangles.size() << '\n';

    for(boost::unordered_set<Tri>::iterator it = dual_triangles.begin();
                                            it != dual_triangles.end(); ++it)

    {
      const Tri& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
//      out << is_triangle_intersected(tr, dual_edges) << '\n';
      out << " 0" << '\n';
    }
    out << "End" << std::endl;
  }

  void output_grid_data_and_dual(const std::string str_base)
  {
    output_grid(str_base);
    output_straight_dual(str_base);
  }

  void output_mapping(const std::size_t seed_id,
                      const Point_2& mapped_centroid,
                      const std::vector<Point_2>& mapped_points,
                      const int counter) const
  {
    K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();

    // outputs the mapped points for a given seed_id
    std::ostringstream filename, filename_bb;
    filename << "mapping_" << seed_id << "_" << counter << ".mesh";
    filename_bb << "mapping_" << seed_id << "_" << counter << ".bb";
    std::ofstream out(filename.str().c_str());
    std::ofstream out_bb(filename_bb.str().c_str());

    std::size_t vertices_count = 0;
    std::ostringstream vertices_out, vertices_out_bb;

    // don't really wanna output all the vertices so keep a renumbering map
    std::map<std::size_t, std::size_t> local;

    // first output the 'real' grid points (mapped)
    for(std::size_t i=0, ps=points.size(); i<ps; ++i)
    {
      const Grid_point& gp = points[i];
      if(gp.closest_seed_id != seed_id)
        continue;

      local[gp.index] = vertices_count++;
      vertices_out << mapped_points[gp.index] << " " << gp.is_Voronoi_vertex << std::endl;
      vertices_out_bb << sqd(mapped_centroid, mapped_points[gp.index]) << std::endl;
    }

#ifdef COMPUTE_PRECISE_VOR_VERTICES
    // then output the 'virtual' grid points if there are any
    const Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];

    for(std::size_t i=0, vvs=Vor_vertices.size(); i<vvs; ++i)
    {
      const Grid_point& gp = Vor_vertices[i];
      CGAL_assertion(gp.closest_seed_id == seed_id);
      local[gp.index] = vertices_count++;

      vertices_out << mapped_points[gp.index] << " " << gp.is_Voronoi_vertex << std::endl;
      vertices_out_bb << sqd(mapped_centroid, mapped_points[gp.index]) << std::endl;
    }
#endif

    std::size_t triangles_count = 0;
    std::ostringstream triangles_out;

    //output the triangles for 'real' vertices
    for(std::size_t i=0, ts=triangles.size(); i<ts; ++i)
    {
      const Tri& tr = triangles[i];
      if(points[tr[0]].closest_seed_id != seed_id ||
         points[tr[1]].closest_seed_id != seed_id ||
         points[tr[2]].closest_seed_id != seed_id)
        continue;

      triangles_count++;
      triangles_out << local[tr[0]]+1 << " " << local[tr[1]]+1 << " " << local[tr[2]]+1 // +1 due to medit
                    << " 1" << '\n';
    }

    // output everything neatly to the file
    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 2" << '\n';
    out << "Vertices" << '\n';
    out << vertices_count << '\n';
    out << vertices_out.str().c_str() << '\n';

    out_bb << "2 1 " << points.size() << " 2" << '\n';
    out_bb << vertices_out_bb.str().c_str() << '\n';
    out_bb << "End" << '\n';

    out << "Triangles" << '\n';
    out << triangles_count << '\n';
    out << triangles_out.str().c_str() << '\n';

#ifdef COMPUTE_PRECISE_VOR_VERTICES
    // output edges between the Voronoi vertices
    std::size_t edges_count = Vor_vertices.size();

    out << "Edges" << '\n';
    out << edges_count << '\n';
    for(std::size_t i=0; i<edges_count; ++i)
    {
      std::size_t j = (i+1) % edges_count;
      out << local[Vor_vertices[i].index]+1 << " "
          << local[Vor_vertices[j].index]+1 << " 0" << '\n'; // +1 due to medit
    }
#endif

    out << "End" << std::endl;
  }

  Base_mesh()
    :
      points(),
      triangles(),
      triangles_buffer(),
      trial_points(),
      edge_bisectors(),
      dual_edges(),
      dual_triangles(),
      refinement_point(NULL),
      Voronoi_vertices()
  { }
};

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
  std::size_t anc_id = ancestor;
  while(anc_id != static_cast<std::size_t>(-1))
  {
    ++i;
    const Grid_point& anc = bm->points[anc_id];
    anc_id = anc.ancestor;

    CGAL_assertion(i < bm->points.size());
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
    std::cout << *cit << " ";
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

void Grid_point::remove_from_neighbors(const std::size_t n)
{
  Point_set::iterator it = neighbors.find(n);
  if(it != neighbors.end()) // might already have been removed from the other side
    neighbors.quick_erase(it);
//  else
//  {
//    std::cout << "WARNING: call to remove_from_neighbors didn't find the neighbor  (";
//    std::cout << n << " from " << index << ")" << std::endl;
//  }
}

void Grid_point::remove_from_incident_triangles(const std::size_t i)
{
  std::list<std::size_t>::iterator it = std::find(incident_triangles.begin(),
                                                  incident_triangles.end(), i);
  CGAL_precondition(it != incident_triangles.end());
  incident_triangles.erase(it);
}

void Grid_point::mark_descendants(std::size_t& count,
                                  boost::unordered_set<std::size_t>& potential_parents)
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
    cp.mark_descendants(count, potential_parents);
  }

#ifdef ADD_PARENTS_AFTER_DESCENDANT_RESET
  Point_set::iterator nit = neighbors.begin();
  Point_set::iterator nend = neighbors.end();

  for(; nit!=nend; ++nit)
  {
    const Grid_point& cp = bm->points[*nit];
    if(cp.state != ORPHAN)
      potential_parents.insert(*nit);
  }
#endif
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
  boost::unordered_set<std::size_t> potential_parents;
  mark_descendants(count, potential_parents);
  reset_descendants();

//  std::cout <<  count << " children in the lineage" << std::endl;

#ifdef ADD_PARENTS_AFTER_DESCENDANT_RESET
  // add the parents to the points we have to spread from...
//  std::cout << "add to trial points (before): " << bm->trial_points.size() << std::endl;
  typename boost::unordered_set<std::size_t>::iterator sit = potential_parents.begin();
  typename boost::unordered_set<std::size_t>::iterator lend = potential_parents.end();
  for(; sit!=lend; ++sit)
  {
    Grid_point& gp = bm->points[*sit];
    if(gp.state != TRIAL)
    {
//      std::cout << "push pot parent: " << *sit << std::endl;
      bm->trial_points.push_back(&gp);
      gp.state = TRIAL;
    }
  }

//  std::cout << "add to trial points (after): " << bm->trial_points.size() << std::endl;
  std::push_heap(bm->trial_points.begin(), bm->trial_points.end(),
                 Grid_point_comparer<Grid_point>());
#endif
}

bool Grid_point::compute_closest_seed(const std::size_t anc_id,
                                      const bool overwrite)
{
#if (verbose > 20)
  std::cout << "closest seed at " << index << " from " << anc_id;
  std::cout << " (current d/a: " << distance_to_closest_seed << " " << ancestor;
  std::cout << " seed: " << closest_seed_id << ")" << std::endl;
#endif

  CGAL_precondition(static_cast<std::size_t>(anc_id) < bm->points.size());
  const Grid_point& anc = bm->points[anc_id];
  CGAL_precondition(anc.ancestor_path_length == anc.anc_path_length());

      if(anc.state != KNOWN)
    std::cout << "WARNING: potential ancestor is not KNOWN" << std::endl;

  if(anc.find_in_ancestry(index))
    return false;

  FT d = FT_inf;

  // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
  // stored as 'this', ancestor1, ancestor2, etc.
  boost::array<std::size_t, k+1> ancestor_path;
  for(int i=0; i<k+1; ++i)
    ancestor_path[i] = -1;
  ancestor_path[0] = this->index;
  CGAL_assertion(ancestor_path[0] != static_cast<std::size_t>(-1));

  std::size_t n_curr_ancestor = anc_id;
  std::size_t best_i = -1;
  for(int i=1; i<=k; ++i)
  {
    // add the new segment to the ancestor path
    ancestor_path[i] = n_curr_ancestor;
    const Grid_point& curr_ancestor = bm->points[n_curr_ancestor];

    Vector2d ancestor_edge;
    ancestor_edge(0) = point.x() - curr_ancestor.point.x();
    ancestor_edge(1) = point.y() - curr_ancestor.point.y();
    FT ancestor_edge_length = ancestor_edge.norm();
    Vector2d normalized_anc_edge = ancestor_edge / ancestor_edge_length;

    // compute the distance for the current depth (i)
    FT dist_to_ancestor = 0.;
    for(int j=0; j<i; ++j) // we add a part for each edge in the path
    {
      // get the metric for the current edge

      CGAL_assertion(ancestor_path[j] < bm->points.size() &&
                     ancestor_path[j+1] < bm->points.size());
      const Grid_point& e0 = (j==0)?*this:bm->points[ancestor_path[j]];
      const Grid_point& e1 = bm->points[ancestor_path[j+1]];

      const Metric& m0 = e0.metric;
      const Metric& m1 = e1.metric;

      Vector2d curr_edge;
      curr_edge(0) = e0.point.x() - e1.point.x();
      curr_edge(1) = e0.point.y() - e1.point.y();

      // interpolate between both metric and transform the normalized edge
      // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
      const Eigen::Matrix2d& f = get_interpolated_transformation(m0, m1);
      Vector2d transformed_anc_edge = f*normalized_anc_edge;
      FT l = transformed_anc_edge.norm(); // length of the normalized anc edge in the metric

#if 0//def USE_MATHIJS_DISTANCE
      FT nl = ancestor_edge.transpose() * f * f * curr_edge;
      CGAL_precondition(l != 0.);
      dist_to_ancestor += nl / ((f * ancestor_edge).norm());
#else
      FT sp = curr_edge.dot(normalized_anc_edge);
      dist_to_ancestor += sp * l;
#endif
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
    std::cout << " to " << d << " from " << anc_id << " (" << anc.closest_seed_id
              << ")" << std::endl;
#endif
    if(overwrite)
    {
      // remove 'index' from the previous ancestor's children (if needed)
      if(ancestor != static_cast<std::size_t>(-1))
        bm->points[ancestor].remove_from_children(index);

      ancestor = anc_id;
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
        Grid_point& gp = bm->points[*(children.begin())];
        gp.deal_with_descendants();

        // checks for circular ancestry
        CGAL_postcondition(gp.anc_path_length() == gp.ancestor_path_length);
      }
#endif
    }
    return true;
  }
  return false;
}

PQ_state Grid_point::update_neighbors_distances(std::vector<Grid_point*>& trial_pq) const
{
  // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (verbose > 10)
  std::cout << "update neighbors of " << index << std::endl;
#endif
  CGAL_assertion(state == KNOWN);

  PQ_state pqs_ret = NOTHING_TO_DO;
  Point_set::const_iterator it = neighbors.begin(),
                            end = neighbors.end();
  for(; it!=end; ++it)
  {
    Grid_point& gp = bm->points[*it];

#if (verbose > 12)
    std::cout << "neighbor: " << *it << " has state: " << gp.state << std::endl;
#endif

    if(gp.state == KNOWN)
      continue;
    else if(gp.state == TRIAL)
    {
      // note that we don't insert in trial_pq since it's already in
      if(gp.compute_closest_seed(this->index))
        pqs_ret = REBUILD_TRIAL;
    }
    else
    {
      CGAL_assertion(gp.state == FAR || gp.state == ORPHAN);
      if(gp.compute_closest_seed(this->index))
      {
        gp.state = TRIAL;
        trial_pq.push_back(&gp);
        std::push_heap(trial_pq.begin(), trial_pq.end(), Grid_point_comparer<Grid_point>());
      }
    }
  }
  return pqs_ret;
}

FT Grid_point::distortion_to_seed() const
{
  FT gamma = 1.;
  //    std::cout << "init gamma: " << gamma << " " << index << std::endl;
  const Grid_point* curr = this;
  std::size_t anc_id = ancestor;

  while(anc_id != static_cast<std::size_t>(-1))
  {
    const Grid_point& anc = bm->points[anc_id];
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
    anc_id = anc.ancestor;
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
  is_seed_holder = -1;
  distance_to_closest_seed = FT_inf;
  closest_seed_id = -1;
  ancestor = -1;
  depth = 0;
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
  is_seed_holder = seed_id;
  ancestor = -1;
}
bool Grid_point::operator==(const Grid_point& gp) const
{
  // the only thing that matters is the point
  return (point == gp.point);
}

Grid_point::Grid_point()
  :
    bm(NULL),
    point(),
    index(static_cast<std::size_t>(-1)),
    metric(),
    neighbors(),
    incident_triangles(),
    is_on_domain_border(false),
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

Grid_point::Grid_point(Base_mesh* bm_,
                       const Point_2& point_,
                       const std::size_t index_,
                       const bool is_on_domain_border_)
  :
    bm(bm_),
    point(point_),
    index(index_),
    metric(mf->compute_metric(point)),
    neighbors(),
    incident_triangles(),
    is_on_domain_border(is_on_domain_border_),
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

Grid_point::Grid_point(const Grid_point& gp)
  :
    bm(gp.bm),
    point(gp.point),
    index(gp.index),
    metric(gp.metric),
    neighbors(gp.neighbors),
    incident_triangles(gp.incident_triangles),
    is_on_domain_border(gp.is_on_domain_border),
    state(gp.state),
    is_seed_holder(gp.is_seed_holder),
    distance_to_closest_seed(gp.distance_to_closest_seed),
    closest_seed_id(gp.closest_seed_id),
    ancestor(gp.ancestor),
    depth(gp.depth),
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
    return gp.point;
  }
};

void draw_metric_field(const MF mf, const Base_mesh& bm)
{
  typedef boost::transform_iterator<Point_extracter,
      std::vector<Grid_point>::const_iterator> Extracted_iterator;
  mf->draw(Extracted_iterator(bm.points.begin(), Point_extracter()),
           Extracted_iterator(bm.points.end(), Point_extracter()));
}

int main(int, char**)
{
  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);
  double duration;
  start = std::clock();
  std::srand(0);

//  mf = new Euclidean_metric_field<K>(1., 1.);
  mf = new Custom_metric_field<K>();
//  mf = new External_metric_field<K>("freefem.mesh", "freefem.sol");

  generate_grid();

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "generated canvas: " << duration << std::endl;

  Base_mesh bm;
  bm.initialize_base_mesh();

//  draw_metric_field(mf, bm);
//  exit(0);

  initialize_seeds();
  bm.locate_and_initialize_seeds();
  bm.spread_distances(true/*use_dual_shenanigans*/);

  if(n_refine > 0)
  {
    bm.output_grid_data_and_dual(str_base_mesh + "_pre_ref");

    for(int i=0; i<n_refine; ++i)
    {
      // can't compute the dual while spreading if we're refining since we have already a layer
      // of paint laying on the canvas...
      bool successful_insert = bm.refine_seeds_with_self_computed_ref_point();

      if(i%100 == 0)
      {
        std::ostringstream out;
        out << "ref_" << seeds.size();
        bm.output_straight_dual(out.str());
      }

      if(!successful_insert)
        break;
    }

    for(std::size_t i=0, ps=bm.points.size(); i<ps; ++i)
      bm.points[i].state = KNOWN;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "End refinement: " << duration << std::endl;
  }

  bm.output_grid_data_and_dual(str_base_mesh + "_tr");

  // optimize stuff
  if(max_opti_n > 0)
  {
    bm.optimize_seeds();
    bm.output_grid_data_and_dual("optimized_" + str_base_mesh + "_tr");
  }

  ignore_children = true;
  bm.compute_geodesics(str_base_mesh);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "duration: " << duration << std::endl;
}
