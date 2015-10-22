// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

#include <CGAL/Starset.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Mesh_2/Sizing_field_2.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <iostream>
#include <functional>
#include <limits>
#include <map>
#include <vector>
#include <set>
#include <utility>

#define COMPUTE_PRECISE_VOR_VERTICES

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;

using namespace CGAL::Anisotropic_mesh_2;

typedef typename K::Vector_2                                 Vector;
typedef Metric_base<K>                                       Metric;

// 'bit ugly to have the metric field running around here but we need to pass
// it around and around and around and it's faster this way !
typedef typename CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
//typedef typename CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;
MF mf;

// stuff to generate the grid using Mesh_2 since Aniso_mesh_2 is a turtle.
// Need to define our own refinement criterion based on metric shenanigans

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
    FT base = 10.;

    const Metric& m = mf->compute_metric(p);
    const FT width = 1./(m.get_max_eigenvalue());

    CGAL_assertion(m.get_max_eigenvalue() >= m.get_min_eigenvalue());

    FT discretization = 6.;
    FT metric_based_size = width / discretization;

    std::cout << "sizing field: " << base << " " << m.get_max_eigenvalue() << " " << metric_based_size << std::endl;

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

template<typename V>
struct Vertex_map_comparator
{
  bool operator()(const V& v1, const V& v2) const
  {
    if(v1.point().x() == v2.point().x())
      return v1.point().y() < v2.point().y();
    return v1.point().x() < v2.point().x();
  }
};

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

  std::map<typename CDT::Vertex, std::size_t,
           Vertex_map_comparator<typename CDT::Vertex> > vertex_map;

  std::ofstream out((str_base + ".mesh").c_str());
  out << std::setprecision(17);
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 2" << std::endl;
  out << "Vertices" << std::endl;
  out << n << std::endl;

  std::size_t v_counter = 1;
  for(typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                                             vit != cdt.finite_vertices_end(); ++vit)
  {
    vertex_map[*vit] = v_counter++;
    typename CDT::Vertex_handle vh = vit;
    out << *vit << " " << is_vertex_on_border<CDT>(cdt, vh) << std::endl;
  }

  out << "Triangles" << std::endl;
  out << m << std::endl;

  for(typename CDT::Face_iterator fit = cdt.finite_faces_begin();
                                  fit != cdt.finite_faces_end(); ++fit)
  {
    for(int i=0; i<=cdt.dimension(); ++i)
      out << vertex_map[*(fit->vertex(i))] << " ";
    out << "0" << std::endl;
  }
  out << "End" << std::endl;
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

  Vertex_handle va = cdt.insert(Point(-0.6,-0.6));
  Vertex_handle vb = cdt.insert(Point(0.3,-0.6));
  Vertex_handle vc = cdt.insert(Point(0.3,0.3));
  Vertex_handle vd = cdt.insert(Point(-0.6,0.3));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  FT shape = 0.125;
  FT distortion = 1.05; // > 1 to use
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(shape, distortion));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  output_cdt_to_mesh(cdt, "input_base_mesh");
}

// -----------------------------------------------------------------------------
// Geodesic stuff now !
// -----------------------------------------------------------------------------

typedef typename K::Point_2                                  Point_2;
typedef typename K::Point_3                                  Point_3;
typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;
typedef typename Eigen::Matrix<double, 2, 1>                 Vector2d;

typedef typename K::Segment_2                                Segment;
typedef typename K::Triangle_2                               Triangle_2;
typedef typename K::Segment_3                                Segment_3;
typedef typename K::Triangle_3                               Triangle_3;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_2                             EPoint;
typedef typename KExact::Segment_2                           ESegment;
typedef typename KExact::Triangle_2                          ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

#define REFINE_GRID
#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 6
const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 1;

// the metric field and the seeds
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

// seeds
const std::string str_seeds = "super_dense_base_mesh_tr_dual";

// base mesh
const std::string str_base_mesh = "rough_base_mesh";
CGAL::Bbox_2 base_mesh_bbox;

// refinement
int n_refine = 0;
std::size_t min_ancestor_path_length = 6;

// optimization
int max_opti_n = 0;
int max_depth = -1; // how precise are the Voronoi vertices (bigger = better)

//debug & info
int known_count=0, trial_count=0, far_count=0;
std::clock_t start;

// -----------------------------------------------------------------------------
// Utility functions below
// -----------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& os, const Simplex& s)
{
  typename Simplex::iterator it = s.begin();
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
      for(std::size_t i=0; i<lower_level_combis.size(); ++i)
        lower_level_combis[i][k_max-k] = *(sit);

      combis.insert(combis.end(), lower_level_combis.begin(),
                                  lower_level_combis.end());
    }
  }
  return combis;
}

bool is_triangle_intersected(const Tri& tri, const std::set<Edge>& edges)
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
    for(std::set<Edge>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
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
  typename K::Compute_area_2 o;
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

  typename K::Compute_area_2 o;
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
  FAR
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL
};

struct Base_mesh;

struct Grid_point
{
  typedef boost::unordered_set<std::size_t> Neighbors;

  // immuable stuff
  Base_mesh* bm;
  Point_2 point;
  std::size_t index;
  Metric metric;
  Neighbors neighbors;
  Neighbors border_neighbors;
  std::list<std::size_t> incident_triangles;
  bool is_on_domain_border;

  // stuff that depends on the seeds
  FMM_state state;
  FT distance_to_closest_seed;
  std::size_t closest_seed_id;
  std::size_t ancestor;
  bool is_Voronoi_vertex;

  void change_state(FMM_state new_state);
  std::size_t ancestor_path_length() const;
  void remove_from_neighbors(const std::size_t n_p);
  void remove_from_border_neighbors(const std::size_t n_p);
  void remove_from_incident_triangles(const std::size_t i);
  bool compute_closest_seed(const std::size_t n_anc, const bool overwrite = true);
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

    typedef typename K::Point_3            Point;
    typedef typename K::Triangle_3         Datum;

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
  typedef typename Primitive::Id                          Primitive_Id;
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
  mutable std::set<Edge> dual_edges;
  mutable std::set<Tri> dual_triangles;

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
    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cout << "WARNING: a new seed is overwriting the closest seed id";
      std::cout << " of a grid point! previous index is: " << gp->closest_seed_id;
      std::cout << " (seeds: " << seeds.size() << ")" << std::endl;
    }

    if(gp->state == TRIAL)
    {
      // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");
    }

    gp->initialize_from_point(dist, seed_id);
    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
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

    typename std::list<Primitive_Id>::const_iterator it = intersections.begin(),
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
          initialize_grid_point(&gp, d, seed_id);
        }
        break; // no need to check the other primitives
      }
    }
    CGAL_assertion(found);
  }

  void build_aabb_tree()
  {
    std::cout << "build aabb tree" << std::endl;
    tree.clear();

    // create an aabb tree of triangles to fasten it
    for(std::size_t i=0; i<triangles.size(); ++i)
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
    std::ifstream in((str_base_mesh + ".mesh").c_str());
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
          if(points[id1].is_on_domain_border && points[id2].is_on_domain_border)
            points[id1].border_neighbors.insert(id2);
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
    for(std::size_t i=0; i<seeds.size(); ++i)
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
    typename std::vector<Grid_point>::iterator it = points.begin();
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

    typename Grid_point::Neighbors::iterator it = gp->border_neighbors.begin();
    typename Grid_point::Neighbors::iterator end = gp->border_neighbors.end();
    for(; it!=end; ++it)
    {
      Grid_point* gq = &(points[*it]);
      CGAL_assertion(gq->is_on_domain_border);
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
    for(std::size_t i=0; i<points.size(); ++i)
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

    if(call_count > max_depth) // don't hardcode stuff like that fixme
      return gp;

    if(gp.closest_seed_id == points[n_q].closest_seed_id)
      return Vor_vertex_on_edge(gp.index, n_r, ++call_count);
    else
    {
      CGAL_assertion(gp.closest_seed_id == points[n_r].closest_seed_id);
      return Vor_vertex_on_edge(gp.index, n_q, ++call_count);
    }
  }

  void visit_neighbors_and_compute(const std::size_t n_p,
                                   std::vector<bool>& visited_points)
  {
    const Grid_point& gp = points[n_p];

    visited_points[n_p] = true;
    std::size_t seed_p = gp.closest_seed_id;

    if(n_p < 4) // corner hack fixme
      Voronoi_vertices[seed_p].push_back(gp);

    typename Grid_point::Neighbors::iterator it = gp.border_neighbors.begin();
    typename Grid_point::Neighbors::iterator end = gp.border_neighbors.end();
    for(; it!=end; ++it)
    {
      const Grid_point& gq = points[*it];
      CGAL_assertion(gq.is_on_domain_border);
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
        visit_neighbors_and_compute(gq.index, visited_points);
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
    for(std::size_t i=0; i<points.size(); ++i)
    {
      if(points[i].is_on_domain_border)
      {
        visit_neighbors_and_compute(i, visited_points);
        break;
      }
    }

    // clean off the virtual points created by Vor_vertex_on_edge
    shave_off_virtual_points(real_points_n);
  }

  void dual_shenanigans(Grid_point* gp,
                        const bool are_Voronoi_vertices_needed)
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

      if(are_Voronoi_vertices_needed)
      {
#ifdef COMPUTE_PRECISE_VOR_VERTICES
        compute_precise_Voronoi_vertex(tri, dual_simplex);
#else
        mark_voronoi_vertices(tri, dual_simplex);
#endif
      }
    }
  }

  void spread_distances_from_one_seed(const int seed_id)
  {
    // spread distances from the seed_id, until it has reached its neighboring
    // seeds

  }

  void spread_distances(const bool use_dual_shenanigans,
                        const bool are_Voronoi_vertices_needed = true /*fixme*/)
  {
    // fixme above : it's pointless to recompute the Voronoi vertices every time
    // but need a good criterion to decide when we need it (typically the last
    // time we spread_distance before optimizing)...

#if (verbose > 5)
    std::cout << "main loop" << std::endl;
#endif
    Grid_point* gp;
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
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << ", " << gp->point.y() << ") ";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from "
                << gp->closest_seed_id << " ancestor : " << gp->ancestor << std::endl;
#endif

      gp->state = KNOWN;
      known_count++; // tmp --> use change_state()

      if(use_dual_shenanigans)
        dual_shenanigans(gp, are_Voronoi_vertices_needed);

      PQ_state pqs = gp->update_neighbors_distances(trial_points);
      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
      is_t_empty = trial_points.empty();
    }

    debug();

    std::cout << "End of spread_distances. time: ";
    std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;

    if(are_Voronoi_vertices_needed && use_dual_shenanigans)
    {
      // we have only computed the voronoi vertices corresponding to the dual of 3-simplex
      // to get the full cell we need to know the intersection of the voronoi cell
      // with the border of the domain
#ifdef COMPUTE_PRECISE_VOR_VERTICES
      compute_precise_Voronoi_vertices_on_border();
#else
      mark_voronoi_vertices_on_border();
#endif
      std::cout << "After Voronoi border: ";
      std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
    }
  }

  void debug()
  {
    CGAL_assertion(trial_points.empty() );
    bool failed = false;

    for(std::size_t i=0; i<points.size(); ++i)
      points[i].state = KNOWN;

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];

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
        for(std::size_t j=0; j<points.size(); ++j)
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

  std::size_t probe_edge_midpoints_map(const std::size_t n_p,
                                       const std::size_t n_q,
                                       std::map<std::pair<std::size_t, std::size_t>,
                                       std::size_t>& edge_mid_points,
                                       bool create_if_not_found = true)
  {
    typedef std::map<std::pair<std::size_t, std::size_t>, std::size_t> Edge_mpts_map;

    std::size_t n_pq = -1;

    bool ordered = (n_p < n_q);
    int nmin = ordered ? n_p : n_q;
    int nmax = ordered ? n_q : n_p;

    std::pair<typename Edge_mpts_map::iterator, bool> is_insert_succesful =
        edge_mid_points.insert(std::make_pair(std::make_pair(nmin, nmax), -1)); // -1 is a placeholder

    if(!is_insert_succesful.second) // already exists in the map
      n_pq = is_insert_succesful.first->second;
    else if(create_if_not_found) // the mid point hasn't been computed yet and we want to build it
    {
      Point_2 mid_point = CGAL::barycenter(points[n_p].point, 0.5,
                                           points[n_q].point, 0.5);
      bool border_info = (points[n_p].is_on_domain_border && points[n_q].is_on_domain_border);
      n_pq = points.size();

      Grid_point gm(this, mid_point, n_pq, border_info);
      points.push_back(gm);
      is_insert_succesful.first->second = n_pq;

      if(border_info)
      {
        points[n_p].border_neighbors.insert(n_pq);
        points[n_q].border_neighbors.insert(n_pq);
        points.back().border_neighbors.insert(n_p);
        points.back().border_neighbors.insert(n_q);
      }
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
               std::map<std::pair<std::size_t, std::size_t>, std::size_t>& m)
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
               std::map<std::pair<std::size_t, std::size_t>, std::size_t>& m)
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
               std::map<std::pair<std::size_t, std::size_t>, std::size_t>& m)
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

    // tedious function so it's written in a very clear way... it's, as usual,
    // not the most efficient but the most straight forward

    // can't use ref or pointers to 'points' since they are regularly invalidated
    // by point insertions

    // todo
    // make it nice with subfunctions and efficient stuff when cleaning up the code

    std::set<std::size_t> seeds_to_refine;
    std::set<std::size_t> seeds_to_rebuild;
    std::vector<bool> refined_triangles(triangles.size(), false);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_mid_points;

    // determine which seeds need to be refined
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      const Tri& triangle = triangles[i];
      std::set<std::size_t> colors;

      for(int j=0; j<3; ++j)
        colors.insert(points[triangle[j]].closest_seed_id);

      if(colors.size() > 1) // bisector of at least an edge or higher dimensional simplex
      {
        for(int j=0; j<3; ++j)
        {
          const Grid_point& gp = points[triangle[j]];
          std::size_t ap_length = gp.ancestor_path_length();
          if(ap_length < min_ancestor_path_length)
          {
            seeds_to_refine.insert(gp.closest_seed_id);
            std::cout << "the canvas is thin (" << ap_length << ") at : "
                      << gp.index << " (" << gp.point
                      << ") dual simplex of size: " << colors.size()
                      << " and cellid: " << gp.closest_seed_id << std::endl;
          }
        }
      }
    }

    if(seeds_to_refine.empty())
      return false;

    std::cout << points.size() << " points before grid refinement" << std::endl;
    std::cout << triangles.size() << " triangles before grid refinement" << std::endl;

    // base mesh refining, not seeds refining !
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      Tri& triangle = triangles[i];

      // the triangle must be refined if at least one of its colors is in seeds_to_refine
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
      // the index 'i' will correspond to the new triangle in the middle
      // and we must build three new triangles

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

      // clean off the neighbors, border neighbors and incident triangles info
      points[n_p].remove_from_incident_triangles(i);
      points[n_q].remove_from_incident_triangles(i);
      points[n_r].remove_from_incident_triangles(i);

      points[n_p].remove_from_neighbors(points[n_q].index);
      points[n_p].remove_from_neighbors(points[n_r].index);
      points[n_q].remove_from_neighbors(points[n_p].index);
      points[n_q].remove_from_neighbors(points[n_r].index);
      points[n_r].remove_from_neighbors(points[n_p].index);
      points[n_r].remove_from_neighbors(points[n_q].index);

      if(points[n_p].is_on_domain_border && points[n_q].is_on_domain_border)
      {
        points[n_p].remove_from_border_neighbors(points[n_q].index);
        points[n_q].remove_from_border_neighbors(points[n_p].index);
      }

      if(points[n_q].is_on_domain_border && points[n_r].is_on_domain_border)
      {
        points[n_q].remove_from_border_neighbors(points[n_r].index);
        points[n_r].remove_from_border_neighbors(points[n_q].index);
      }

      if(points[n_p].is_on_domain_border && points[n_r].is_on_domain_border)
      {
        points[n_p].remove_from_border_neighbors(points[n_r].index);
        points[n_r].remove_from_border_neighbors(points[n_p].index);
      }

      // points[n_r]ab the colors that we will need to paint again
      seeds_to_rebuild.insert(points[n_p].closest_seed_id);
      seeds_to_rebuild.insert(points[n_q].closest_seed_id);
      seeds_to_rebuild.insert(points[n_r].closest_seed_id);

      // points[n_r]ab (or create) the edge points
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
    for(std::size_t i=0; i<triangles.size(); ++i)
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

    // some quick debug on incidency and adjacency
    for(std::size_t i=0; i<triangles.size(); ++i)
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

    build_aabb_tree();

#ifdef USE_FULL_REBUILDS
    reset();
    locate_and_initialize_seeds();
#else
    // if you allow other cells to be modified without reseting them fully,
    // you can create orphans, which must NOT happen...
    refresh_grid_point_states();
    for(std::size_t i=0; i<points.size(); ++i)
      points[i].state = KNOWN;

    // go an extra step (is it necessary) for points whose colors we are spreading:
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point& gp = points[i];
      if(seeds_to_rebuild.find(gp.closest_seed_id) != seeds_to_rebuild.end())
      {
//        std::cout << "reseting: " << i << "(id: " << gp.closest_seed_id << ")" << std::endl;
        gp.reset();
      }
    }

    std::cout << "must rebuild: ";
    std::set<std::size_t>::iterator it = seeds_to_rebuild.begin(),
                                    end = seeds_to_rebuild.end();
    for(; it!=end; ++it)
    {
      const Point_2& p = seeds[*it];
      std::cout << *it << " ";
      locate_and_initialize(p, *it);
    }
    std::cout << std::endl;
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

    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point& gp = points[i];
      gp.state = FAR;
      gp.distance_to_closest_seed = FT_inf;
      gp.closest_seed_id = -1;
      gp.ancestor = -1;
      gp.is_Voronoi_vertex = false;
    }

    known_count = 0;
    trial_count = 0;
    far_count = points.size();

    clear_dual();
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

    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << points.size() << std::endl;
    out_bb << "2 1 " << points.size() << " 2" << std::endl;
//    out_distortion_bb << "2 1 " << points.size() << " 2" << std::endl;

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << i+1 << std::endl;

      out_bb << gp.distance_to_closest_seed << std::endl;
//      out_bb << gp.closest_seed_id << std::endl;
//      out_distortion_bb << gp.distortion_to_seed() << std::endl;
//      out_bb << (gp.is_Voronoi_vertex?"1":"0") << std::endl;
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      std::set<std::size_t> materials;
      for(int j=0; j<3; ++j)
      {
        materials.insert(points[triangles[i][j]].closest_seed_id);
        out << triangles[i][j]+1 << " ";
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }

//    out_distortion_bb << "End" << std::endl;
    out_bb << "End" << std::endl;
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

  struct PQ_finder_comparer
  {
    PQ_entry entry;

    bool operator()(const PQ_entry& right)
    {
      const Simplex& s1 = entry.get<0>();
      const Simplex& s2 = right.get<0>();

      if(s1.size() != s2.size())
        return false;

      typename Simplex::const_iterator it1 = s1.begin();
      typename Simplex::const_iterator it2 = s2.begin();

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
    typename PQ::iterator it = std::find_if(queue.begin(), queue.end(),
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

    FT max_size = 0.; // <= 0 means unused
    FT max_distortion = 1.; // <= 1 means unused
    bool intersection_ref = true;
    FT min_qual = 0.5; // <= 0 means unsused

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
    typename Tri::const_iterator it = grid_tri.begin();
    typename Tri::const_iterator end = grid_tri.end();
    for(; it!=end; ++it)
     test_simplex(&(points[*it]), dual_simplex);
  }

  void add_simplex_to_triangulation(const Simplex& dual_simplex)
  {
    CGAL_assertion(dual_simplex.size() <= 3);

    typename Simplex::const_iterator it = dual_simplex.begin();
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

    std::size_t n_anc = gp.ancestor;
    while(n_anc != static_cast<std::size_t>(-1))
    {
      const Grid_point& anc = points[n_anc];
      geodesic_path.push_front(&anc);
      n_anc = anc.ancestor;
    }

    // determine gq
    // Since we have cleared out the case gp = seed, we know that :
    // - geodesic_path.size() > 0
    // - even if we take gq_index = 0, we have a non degenerate vector
    std::size_t gq_index = geodesic_path.size()-1; // geodesic_path.size() / 5; // TMP

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

    for(std::size_t i=0; i<points.size(); ++i)
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
    for(std::size_t seed_id=0; seed_id<Voronoi_vertices.size(); ++seed_id)
      precise_Vor_vertices_count += Voronoi_vertices[seed_id].size();

    // resize the mapped points to add the precise Vor vertices
    mapped_points.resize(points.size() + precise_Vor_vertices_count);

    // to number the precise Vor vertices (they are not base mesh points)
    std::size_t precise_Vor_vertex_index = points.size();

    for(std::size_t seed_id=0; seed_id<Voronoi_vertices.size(); ++seed_id)
    {
      Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];
      for(std::size_t i=0; i<Vor_vertices.size(); ++i)
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
    typename K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();

    for(std::size_t i=0; i<points.size(); ++i)
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
    for(std::size_t i=0; i<triangles.size(); ++i)
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
    for(std::size_t i=0; i<triangles.size(); ++i)
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
        {
          const Grid_point& gp0 = points[triangle[i]];
          const Grid_point& gp1 = points[triangle[(i+1)%3]];
          mid_pts[i] = Vor_vertex_on_edge(gp0.index, gp1.index);
        }

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
    for(std::size_t i=0; i<seed_centroids.size(); ++i)
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

    for(std::size_t i=0; i<triangles.size(); ++i)
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

    for(std::size_t i=0; i<triangles.size(); ++i)
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

    typename Voronoi_vertices_container::const_iterator it = Vor_vertices.begin();
    typename Voronoi_vertices_container::const_iterator next_it = ++(Vor_vertices.begin());
    typename Voronoi_vertices_container::const_iterator end = Vor_vertices.end();
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

    typename Voronoi_vertices_container::const_iterator it = Vor_vertices.begin();
    typename Voronoi_vertices_container::const_iterator next_it = ++(Vor_vertices.begin());
    typename Voronoi_vertices_container::const_iterator end = Vor_vertices.end();
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
    Point_2 mapped_centroid = compute_mapped_centroid_with_mapped_Vor_vertices(seed_id, mapped_points);
    std::cout << "(mapped) centroid coordinates from mapped Vor vertices " << mapped_centroid  << std::endl;

    // -------------------------------------------------------------------------
    std::cout << "UNMAPPING POSSIBILITIES :" << std::endl;

    Point_2 closest_grid_point_centroid = compute_centroid_as_closest_grid_point(seed_id, mapped_centroid, mapped_points);
    std::cout << "unmapped centroid as closest unmapped grid point " << closest_grid_point_centroid << std::endl;

    // disabled till degenerate mapped triangles are fixed fixme
//    Point_2 bar_centroid = compute_centroid_with_barycentric_info(seed_id, mapped_centroid, mapped_points);
//    std::cout << "unmapped centroid from mapped polygon (barycentric) " << bar_centroid << std::endl;

    // --
      std::cout << "FOR REFERENCE, WITHOUT TRANSFORMATION GIVES :" << std::endl;

      Point_2 alt_centroid_1 = compute_centroid_with_grid_triangles(seed_id);
      std::cout << "centroid with grid triangles " << alt_centroid_1 << std::endl;

      Point_2 alt_centroid_2 = compute_centroid_with_grid_triangles_precomputed(seed_id,
                                                                                centroids);
      std::cout << "centroid with grid triangles " << alt_centroid_2 << std::endl;

      // below is not a good idea if the metric field is not uniform
      Point_2 alt_centroid_3 = compute_centroid_with_voronoi_vertices(seed_id);
      std::cout << "centroid with voronoi vertices " << alt_centroid_3 << std::endl;
    // --

    std::cout << "ALTERNATIVE MAPPED CENTROID COMPUTATION : " << std::endl;

    Point_2 mapped_grid_centroid = compute_centroid_with_mapped_grid_triangles(seed_id, mapped_points);
    std::cout << "(mapped) centroid with mapped grid triangles " << mapped_grid_centroid << std::endl;

    Point_2 closest_grid_point_centroid_alternate =
        compute_centroid_as_closest_grid_point(seed_id, mapped_grid_centroid, mapped_points);
    std::cout << "unmapped centroid as closest unmapped grid point (alternate)"
              << closest_grid_point_centroid_alternate << std::endl;

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
    typename K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();
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

    for(std::size_t i=0; i<points.size(); ++i)
      CGAL_assertion(points[i].state == KNOWN);

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
      for(std::size_t i=0; i<seeds.size(); ++i)
        cumulated_displacement += optimize_seed(i, mapped_points, centroids, counter);

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

    std::size_t length = gp->ancestor_path_length();
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

    for(std::size_t i=0; i<triangles.size(); ++i)
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

      for(std::size_t i=0; i<grid_triangles.size(); ++i)
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

      for(std::size_t i=0; i<visited_status.size(); ++i)
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

    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      const Tri& triangle = triangles[i];
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
    typename K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();

    FT alpha = 4.*std::sqrt(3.);
    FT A = std::abs(typename K::Compute_area_2()(p, q, r));
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
      out << is_triangle_intersected(tr, dual_edges) << std::endl;
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
    typename K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();

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
    for(std::size_t i=0; i<points.size(); ++i)
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

    for(std::size_t i=0; i<Vor_vertices.size(); ++i)
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
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      const Tri& tr = triangles[i];
      if(points[tr[0]].closest_seed_id != seed_id ||
         points[tr[1]].closest_seed_id != seed_id ||
         points[tr[2]].closest_seed_id != seed_id)
        continue;

      triangles_count++;
      triangles_out << local[tr[0]]+1 << " " << local[tr[1]]+1 << " " << local[tr[2]]+1 // +1 due to medit
                    << " 1" << std::endl;
    }

    // output everything neatly to the file
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << vertices_count << std::endl;
    out << vertices_out.str().c_str() << std::endl;

    out_bb << "2 1 " << points.size() << " 2" << std::endl;
    out_bb << vertices_out_bb.str().c_str() << std::endl;
    out_bb << "End" << std::endl;

    out << "Triangles" << std::endl;
    out << triangles_count << std::endl;
    out << triangles_out.str().c_str() << std::endl;

#ifdef COMPUTE_PRECISE_VOR_VERTICES
    // output edges between the Voronoi vertices
    std::size_t edges_count = Vor_vertices.size();

    out << "Edges" << std::endl;
    out << edges_count << std::endl;
    for(std::size_t i=0; i<edges_count; ++i)
    {
      std::size_t j = (i+1) % edges_count;
      out << local[Vor_vertices[i].index]+1 << " "
          << local[Vor_vertices[j].index]+1 << " 0" << std::endl; // +1 due to medit
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

std::size_t Grid_point::ancestor_path_length() const
{
  std::size_t i = 1;
  std::size_t n_anc = ancestor;
  while(n_anc != static_cast<std::size_t>(-1))
  {
    ++i;
    const Grid_point& anc = bm->points[n_anc];
    n_anc = anc.ancestor;
  }
  return i;
}

void Grid_point::remove_from_neighbors(const std::size_t n_p)
{
  typename Neighbors::iterator it = neighbors.find(n_p);
  if(it != neighbors.end()) // might already have been removed from the other side
    neighbors.quick_erase(it);
}

void Grid_point::remove_from_border_neighbors(const std::size_t n_p)
{
  typename Neighbors::iterator it = border_neighbors.find(n_p);
  CGAL_assertion(it != border_neighbors.end()); // no other side here though !
  border_neighbors.quick_erase(it);
}

void Grid_point::remove_from_incident_triangles(const std::size_t i)
{
  std::list<std::size_t>::iterator it = std::find(incident_triangles.begin(),
                                                  incident_triangles.end(), i);
  CGAL_assertion(it != incident_triangles.end());
  incident_triangles.erase(it);
}

bool Grid_point::compute_closest_seed(const std::size_t n_anc,
                                      const bool overwrite)
{
#if (verbose > 20)
  std::cout << "closest seed at " << index << " from " << n_anc;
  std::cout << " (current d/a: " << distance_to_closest_seed << " " << ancestor;
  std::cout << " seed: " << closest_seed_id << ")" << std::endl;
#endif

  CGAL_precondition(n_anc < bm->points.size());
  const Grid_point& anc = bm->points[n_anc];
  CGAL_assertion(anc.state == KNOWN);

  const int k = 8; // depth of the ancestor edge
  FT d = FT_inf;

  // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
  // stored as 'this', ancestor1, ancestor2, etc.
  boost::array<std::size_t, k+1> ancestor_path;
  for(int i=1; i<k+1; ++i)
    ancestor_path[i] = -1;
  ancestor_path[0] = this->index;
  CGAL_assertion(ancestor_path[0] != static_cast<std::size_t>(-1));

  std::size_t n_curr_ancestor = n_anc;
  for(int i=1; i<=k; ++i)
  {
    // add the new segment to the ancestor path
    ancestor_path[i] = n_curr_ancestor;
    const Grid_point& curr_ancestor = bm->points[n_curr_ancestor];

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
#if (verbose > 25)
      std::cout << "new best at " << new_d << " for i: " << i << std::endl;
#endif
      d = new_d;
    }

    // check if we can go any farther up in the ancestor tree
    if(curr_ancestor.ancestor == static_cast<std::size_t>(-1))
      break;

    n_curr_ancestor = curr_ancestor.ancestor;
  }

  if(d < distance_to_closest_seed)
  {
#if (verbose > 20)
    std::cout << "improving " << distance_to_closest_seed << " from " << ancestor << " (" << closest_seed_id << ")";
    std::cout << " to " << d << " from " << n_anc << " (" << anc.closest_seed_id << ")" << std::endl;
#endif
    if(overwrite)
    {
      ancestor = n_anc;
      distance_to_closest_seed = d;
      closest_seed_id = anc.closest_seed_id;
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
  typename Neighbors::const_iterator it = neighbors.begin(),
                                     end = neighbors.end();
  for(; it!=end; ++it)
  {
    Grid_point& gp = bm->points[*it];

    if(gp.state == KNOWN)
      continue;
    else if(gp.state == TRIAL)
    {
      // note that we don't insert in trial_pq since it's already in
      if(gp.compute_closest_seed(this->index))
        pqs_ret = REBUILD_TRIAL;
    }
    else // gp.state == FAR
    {
      CGAL_assertion(gp.state == FAR);
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
  std::size_t n_anc = ancestor;

  while(n_anc != static_cast<std::size_t>(-1))
  {
    const Grid_point& anc = bm->points[n_anc];
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

void Grid_point::reset()
{
  state = FAR;
  distance_to_closest_seed = FT_inf;
  closest_seed_id = -1;
  ancestor = -1;
  is_Voronoi_vertex = false;
}

void Grid_point::initialize_from_point(const FT d,
                                       const std::size_t seed_id)
{
  closest_seed_id = seed_id;
  distance_to_closest_seed = d;
  state = TRIAL;
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
    border_neighbors(),
    incident_triangles(),
    is_on_domain_border(false),
    state(FAR),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    is_Voronoi_vertex(false)
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
    border_neighbors(),
    incident_triangles(),
    is_on_domain_border(is_on_domain_border_),
    state(FAR),
    distance_to_closest_seed(FT_inf),
    closest_seed_id(-1),
    ancestor(-1),
    is_Voronoi_vertex(false)
{ }

Grid_point::Grid_point(const Grid_point& gp)
  :
    bm(gp.bm),
    point(gp.point),
    index(gp.index),
    metric(gp.metric),
    neighbors(gp.neighbors),
    border_neighbors(gp.border_neighbors),
    incident_triangles(gp.incident_triangles),
    is_on_domain_border(gp.is_on_domain_border),
    state(gp.state),
    distance_to_closest_seed(gp.distance_to_closest_seed),
    closest_seed_id(gp.closest_seed_id),
    ancestor(gp.ancestor),
    is_Voronoi_vertex(gp.is_Voronoi_vertex)
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

//  mf = new Euclidean_metric_field<K>(5., 1.);
  mf = new Custom_metric_field<K>();

//  generate_grid();
//  exit(0);

  Base_mesh bm;
  bm.initialize_base_mesh();
//  draw_metric_field(mf, bm);
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

      std::ostringstream out;
      out << "ref_" << seeds.size();
      bm.output_grid_data_and_dual(out.str());

      if(!successful_insert)
        break;
    }

    for(std::size_t i=0; i<bm.points.size(); ++i)
      bm.points[i].state = KNOWN;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "End refinement: " << duration << std::endl;
  }

  bm.output_grid_data_and_dual(str_base_mesh + "_tr");
//  bm.check_edelsbrunner();

  // optimize stuff
  if(max_opti_n > 0)
  {
    bm.optimize_seeds();
    bm.output_grid_data_and_dual("bis_optimized_" + str_base_mesh + "_tr");
  }

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "duration: " << duration << std::endl;
}
