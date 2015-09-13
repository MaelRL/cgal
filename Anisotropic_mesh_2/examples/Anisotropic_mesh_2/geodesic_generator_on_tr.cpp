// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

#include <iostream>
#include <functional>
#include <limits>
#include <map>
#include <vector>
#include <set>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;

using namespace CGAL::Anisotropic_mesh_2;

typedef typename K::Vector_2                                 Vector;

typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Metric                                Metric;

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
class Delaunay_mesh_distortion_criteria_2 :
    public virtual CGAL::Delaunay_mesh_size_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits                          Geom_traits;
  FT distortionbound;

public:
  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>           Base;

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

      q.get<0>() = 0;
      if( distortion_bound != 0 )
      {
        const Metric& ma = mf->compute_metric(pa);
        const Metric& mb = mf->compute_metric(pb);
        const Metric& mc = mf->compute_metric(pc);

        FT gamma_ab = ma.compute_distortion(mb);
        FT gamma_ac = ma.compute_distortion(mc);
        FT gamma_bc = mb.compute_distortion(mc);

        FT max_gamma = (std::max)((std::max)(gamma_ab, gamma_ac), gamma_bc);
        q.get<0>() = max_gamma / distortion_bound;

//        std::cout << max_gamma / distortion_bound << std::endl;

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
      if( this->squared_size_bound != 0 )
      {
        q.get<1>() = max_sq_length / this->squared_size_bound;
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
           const FT size_bound,
           const FT dist_bound,
           const Geom_traits& traits)
      :
        Base::Is_bad(aspect_bound, size_bound, traits),
        distortion_bound(dist_bound)
    { }
  };

  Is_bad is_bad_object() const
  {
    return Is_bad(this->bound(), this->size_bound(), this->distortion_bound(),
                  this->traits /* from the bad class */);
  }

  Delaunay_mesh_distortion_criteria_2(const FT aspect_bound = 0.125,
                                      const FT size_bound = 0.,
                                      const FT distortion_bound = 0.,
                                      const Geom_traits& traits = Geom_traits())
    :
      Base(aspect_bound, size_bound, traits),
      distortionbound(distortion_bound)
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
  FT size = 0.002;
  FT distortion = 1.;
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(shape, size, distortion));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  output_cdt_to_mesh(cdt, "rough_base_mesh");
}

// -----------------------------------------------------------------------------
// Geodesic stuff now !
// -----------------------------------------------------------------------------

typedef typename Star::Point_2                               Point_2;
typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;
typedef typename Eigen::Matrix<double, 2, 1>                 Vector2d;

typedef typename K::Segment_2                                Segment;
typedef typename K::Triangle_2                               Triangle;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_2                             EPoint;
typedef typename KExact::Segment_2                           ESegment;
typedef typename KExact::Triangle_2                          ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

#define FILTER_SEEDS_OUTSIDE_GRID
#define verbose 6
const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 1;

// the metric field and the seeds
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

// base mesh
const std::string str_base_mesh = "rough_base_mesh";
CGAL::Bbox_2 base_mesh_bbox;

//refinement
int n_refine = 0;

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
  const Triangle triangle(p0, p1, p2);

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

void compute_bary_weights(const Point_2&p , const Point_2& a, const Point_2& b, const Point_2& c,
                          FT& lambda_a, FT& lambda_b, FT& lambda_c)
{
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
  std::ifstream in("adapted_base_mesh_tr_dual.mesh");
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

struct Grid_point
{
  typedef boost::unordered_set<Grid_point*> Neighbors;

  // immuable stuff
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
  const Grid_point* ancestor;
  bool is_Voronoi_vertex;

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
        const Eigen::Matrix2d& f = get_interpolated_transformation(m0, m1);
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
    boost::unordered_set<Grid_point*>::const_iterator it = neighbors.begin(),
                                                      end = neighbors.end();
    for(; it!=end; ++it)
    {
      Grid_point* gp = *it;
      if(!gp)
        continue;
      else if(gp->state == KNOWN)
        continue;
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

  FT distortion_to_seed() const
  {
    FT gamma = 1.;
//    std::cout << "init gamma: " << gamma << " " << index << std::endl;
    const Grid_point* curr = this;
    const Grid_point* anc = ancestor;

    while(anc)
    {
      const Metric& m1 = anc->metric;
      const Metric& m2 = curr->metric;
      FT loc_gamma = m1.compute_distortion(m2);
//      std::cout << "loc: " << loc_gamma << " " << anc->index << std::endl;
#if 1
      gamma *= loc_gamma;
#else
      gamma = (std::max)(loc_gamma, gamma);
#endif
//      std::cout << "gamma: " << gamma << std::endl;
      anc = anc->ancestor;
    }

    gamma = (std::min)(25., gamma);

//    std::cout << "final gamma:" << gamma << std::endl;
//    exit(0);
    return gamma;
  }

  void initialize_from_point(const FT d,
                             const std::size_t seed_id)
  {
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
    ancestor = NULL;
  }

  bool operator==(const Grid_point& gp)
  {
    // the only thing that matters is the point
    return (point == gp.point);
  }

  Grid_point(const Point_2& point_,
             const std::size_t index_,
             const bool is_on_domain_border_ = false)
    :
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
      ancestor(NULL),
      is_Voronoi_vertex(false)
  { }

  Grid_point(const Grid_point& gp)
    :
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
};

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

inline std::ostream& operator<<(std::ostream& os,
                                const boost::tuple<Simplex, const Grid_point*, FT>& pqe)
{
  std::cout << "PQ entry element:" << std::endl;
  std::cout << "Simplex: " << pqe.get<0>();
  std::cout << "gp: " << pqe.get<1>()->index << " [" << pqe.get<1>()->point << "] ";
  std::cout << pqe.get<1>()->distance_to_closest_seed << std::endl;
  std::cout << "FT: " << pqe.get<2>() << std::endl;

  return os;
}

struct Base_mesh
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
  std::vector<Tri> triangles;
  Grid_point_vector trial_points;

  // extra info: map[edge e] = triangles incident to e
  std::map<std::pair<std::size_t, std::size_t>, std::list<std::size_t> > edge_incident_tris;

  // Duality
  mutable Edge_bisec_map edge_bisectors;
  mutable std::set<Edge> dual_edges;
  mutable std::set<Tri> dual_triangles;

  // Refinement
  typedef boost::tuple<Simplex, const Grid_point*, FT>           PQ_entry;
//  typedef std::pair<const Grid_point*, FT>                   PQ_entry;
  typedef std::set<PQ_entry, PQ_entry_comparator<PQ_entry> > PQ;

  PQ size_queue;
  PQ distortion_queue;
  PQ intersection_queue;
  PQ quality_queue;

  mutable const Grid_point* refinement_point;
  mutable FT refinement_distance;

  void clear_dual()
  {
    edge_bisectors.clear();
    dual_edges.clear();
    dual_triangles.clear();

    size_queue.clear();
    distortion_queue.clear();
    intersection_queue.clear();
    quality_queue.clear();

    refinement_point = NULL;
    refinement_distance = 0.;
  }

  void initialize_grid_point(Grid_point* gp,
                             const FT dist,
                             const std::size_t seed_id)
  {
    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id;
      std::cerr << " (seeds: " << seeds.size() << ")" << std::endl;
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

    // something pretty todo... (aka triangulation_2's locate once the canvas derives from Tr_2)
    // for now: brutally loop the triangles and orientation tests!

    bool found = false;
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      const Tri& tr = triangles[i];
      const Point_2& p = points[tr[0]].point;
      const Point_2& q = points[tr[1]].point;
      const Point_2& r = points[tr[2]].point;

      CGAL::Sign o1 = CGAL::orientation(p, q, s);
      CGAL::Sign o2 = CGAL::orientation(q, r, s);
      CGAL::Sign o3 = CGAL::orientation(r, p, s);

//      std::cout << i << " " << o1 << " " << o2 << " " << o3 << std::endl;

      if((o1 >= 0 && o2 >=0 && o3 >= 0) ||
         (o1 <= 0 && o2 <=0 && o3 <= 0))
      {
        // just to be sure ! (useless stuff)
        Triangle t(p,q,r);
        Triangle t_p(q,r,s), t_q(r,p,s), t_r(p,q,s);
        FT asd = t_p.area() / t.area();
        FT bsd = t_q.area() / t.area();
        FT csd = t_r.area() / t.area();
        CGAL_assertion(asd >= 0. && bsd >= 0. && csd >= 0.);
#if (verbose > 10)
        std::cout << "locating seed " << seed_id
                  << " point: " << p.x() << " " << p.y() << std::endl;
        std::cout << "found triangle: " << i << std::endl;
        std::cout << tr[0] << " [" << p << "] " << std::endl;
        std::cout << tr[1] << " [" << q << "] " << std::endl;
        std::cout << tr[2] << " [" << r << "] " << std::endl;
#endif
        // end of useless stuff
        if(found)
          continue;

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

        found = true;
        //break;
      }
    }
    CGAL_assertion(found);
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
    std::cout << "base mesh made of " << nv << " vertices" << std::endl;

    for(std::size_t i=0; i<nv; ++i)
    {
      int border_info;
      in >> r_x >> r_y >> border_info;

      if(border_info != 0 && border_info != 1)
        CGAL_assertion(false && "you're not using a mesh with border info");

      Point_2 p(r_x, r_y);
      Grid_point gp(p, i, border_info);
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
          points[id1].neighbors.insert(&points[id2]);
          if(points[id1].is_on_domain_border && points[id2].is_on_domain_border)
            points[tr[i]].border_neighbors.insert(&points[id2]);
        }
        points[id1].incident_triangles.push_back(new_tri_id);
      }
    }

#if (verbose > 5)
    std::cout << "base mesh initialized" << std::endl;
#endif
  }

  void locate_and_initialize_seeds()
  {
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_2& p = seeds[i];
      locate_and_initialize(p, i);
    }
  }

  void print_states() const
  {
    std::cout << "known: " << known_count;
    std::cout << " trial: " << trial_count;
    std::cout << " far: " << far_count << std::endl;
  }

  void dual_shenanigans(Grid_point* gp)
  {
    // mega ugly hack to correctly compute the centroid of Voronoi cells who contain
    // a corner of the domain...
    if(gp->index < 4)
      gp->is_Voronoi_vertex = true;

    CGAL_assertion(gp->state == KNOWN);
    // Check the neighbors of gp for points that have the state KNOWN.
    // In these, consider those that have a closest_seed_id different than gp's
    // those give dual simplices.
    // If we have never seen that dual already, then we have also found which
    // point we will draw to in the geodesic_triangulation

    std::list<std::size_t>::const_iterator it = gp->incident_triangles.begin(),
                                           end = gp->incident_triangles.end();
    for(; it!=end; ++it)
    {
      const Tri& tri = triangles[*it];
      std::set<std::size_t> dual_simplex;

      for(std::size_t j=0; j<3; ++j)
      {
        const Grid_point& gq = points[tri[j]]; // it's possible that gq = gp
        if(gq.state != KNOWN)
          continue;

        std::size_t q_id = gq.closest_seed_id;
        CGAL_assertion(q_id != static_cast<std::size_t>(-1));
        dual_simplex.insert(q_id);

/* shenanigans to build pseudo geodesic triangulations...
        std::size_t p_id = gp->closest_seed_id;
        if(p_id != q_id) // we're on the bisector of [PQ]
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
*/
      }

      std::size_t ds_size = dual_simplex.size();
      if(ds_size == 2)
      {
        // doing it explicitely because it's clearer and there are many subcases atm

        Grid_point& gq = points[tri[0]];
        Grid_point& gr = points[tri[1]];
        Grid_point& gs = points[tri[2]];

        // we accept a point as a Voronoi vertex only if there are 2 vertices
        // on the domain border with different domain colors

        bool gq_on_border = gq.is_on_domain_border;
        bool gr_on_border = gr.is_on_domain_border;
        bool gs_on_border = gs.is_on_domain_border;

        // tedious tree ahead, wear a helmet...
        // it's the most explicit at the moment, is it worth obfuscating it to shorten it ?

        if(gq_on_border) // gq
        {
          if(gr_on_border) // gq && gr
          {
            if(gs_on_border) // gq && gr && gs
            {
              // really rare case: the triangle is a corner AND a bisector...
              // since we only have two colors, there are 2 points with the same color
              // we need to keep the one farthest away from the seed

              if(gq.closest_seed_id == gr.closest_seed_id) // gq & gr have the same color
              {
                if(gq.distance_to_closest_seed > gr.distance_to_closest_seed)
                  gq.is_Voronoi_vertex = true;
                else
                  gr.is_Voronoi_vertex = true;
                gs.is_Voronoi_vertex = true; // gs has the second color
              }
              else if(gq.closest_seed_id == gs.closest_seed_id) // gq & gs have the same color
              {
                if(gq.distance_to_closest_seed > gs.distance_to_closest_seed)
                  gq.is_Voronoi_vertex = true;
                else
                  gs.is_Voronoi_vertex = true;
                gr.is_Voronoi_vertex = true; // gr has the second color
              }
              else // gr & gs have the same color
              {
                if(gr.distance_to_closest_seed > gs.distance_to_closest_seed)
                  gr.is_Voronoi_vertex = true;
                else
                  gs.is_Voronoi_vertex = true;
                gq.is_Voronoi_vertex = true; // gq has the second color
              }
            }
            else // gq && gr && !gs
            {
              if(gq.closest_seed_id != gr.closest_seed_id)
              {
                gq.is_Voronoi_vertex = true;
                gr.is_Voronoi_vertex = true;
              }
            }
          }
          else // gq && !gr
          {
            if(gs_on_border) // gq && !gr && gs
            {
              if(gq.closest_seed_id != gs.closest_seed_id)
              {
                gq.is_Voronoi_vertex = true;
                gs.is_Voronoi_vertex = true;
              }
            }
            // can't do anything with gq && !gr && !gs
          }
        }
        else // !gq
        {
          if(gr_on_border) // !gq && gr
          {
            if(gs_on_border) // !gq && gr && gs
            {
              if(gr.closest_seed_id != gs.closest_seed_id)
              {
                gr.is_Voronoi_vertex = true;
                gs.is_Voronoi_vertex = true;
              }
            }
            // can't do anything with !gq && gr && !gs
          }
          // can't do anything with !gq && !gr (whatever gs is)
        }
      }
      else if(ds_size == 3)
      {
        for(std::size_t j=0; j<3; ++j)
        {
          Grid_point& gq = points[tri[j]]; // it's possible that gq = gp
          gq.is_Voronoi_vertex = true;
//          std::cout << gq.index << " (" << gq.point << ") vor vert for " << gq.closest_seed_id << std::endl;
        }
      }

      // the simplex of size ? is added to the dual simplices
      add_simplex_to_triangulation(gp, dual_simplex);
    }
  }

  void spread_distances()
  {
#if (verbose > 5)
    std::cout << "main loop" << std::endl;
#endif
    Grid_point* gp;

    bool is_t_empty = trial_points.empty();
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
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << ", " << gp->point.y() << ") ";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from " << gp->closest_seed_id << std::endl;
#endif

      gp->state = KNOWN;
      known_count++; // tmp --> use change_state()

      dual_shenanigans(gp);

      PQ_state pqs = gp->update_neighbors_distances(trial_points);

      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

      is_t_empty = trial_points.empty();
    }

    std::cerr << "End of spread_distances. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
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
      boost::unordered_set<Grid_point*>::const_iterator it = gp->neighbors.begin(),
                                                        end = gp->neighbors.end();
      for(; it!=end; ++it)
      {
        const Grid_point* gq = *it;
        if(gq)
          CGAL_assertion(!gp->compute_closest_seed(gq));
      }
    }

    std::cerr << "End of debug. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
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
      gp.ancestor = NULL;
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

  void refine_grid(const Point_2 new_seed)
  {
    vertices_nv = insert_new_seed(new_seed.x(), new_seed.y());
    refresh_grid_point_states(); // we can't spread from the new seed if all states are 'KNOWN'
    locate_and_initialize(new_seed, seeds.size()-1);

    spread_distances();
  }

  bool refine_grid_with_self_computed_ref_point()
  {
    std::cout << "state of the queues: "
              << " size: " << size_queue.size()
              << " distortion: " << distortion_queue.size()
              << " intersection: " << intersection_queue.size()
              << " quality: "  << quality_queue.size() << std::endl;

    PQ_entry best_entry;

    if(!size_queue.empty())
       best_entry = *(size_queue.begin());
    else if(!distortion_queue.empty())
      best_entry = *(distortion_queue.begin());
    else if(!intersection_queue.empty())
      best_entry = *(intersection_queue.begin());
    else if(!quality_queue.empty())
      best_entry = *(quality_queue.begin());

    refinement_point = best_entry.get<1>();

    if(refinement_point)
    {
      std::cout << "naturally we picked : " << refinement_point->index << "[";
      std::cout << refinement_point->point << "]" << std::endl;
      std::cout << "second was: " << best_entry.get<2>() << std::endl;
    }

    if(!refinement_point)
    {
      std::cerr << "Couldn't find a ref point, need more initial points" << std::endl;
      return false;
//      std::exit(EXIT_FAILURE);
    }
    else
      refine_grid(refinement_point->point);
    return true;
  }

  void output_grid(const std::string str_base) const
  {
#if (verbose > 5)
    std::cout << "output grid" << std::endl;
#endif

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
          Triangle tri(seeds[tr[0]], seeds[tr[1]], seeds[tr[2]]);
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

    const Grid_point* anc = gp.ancestor;
    while(anc)
    {
      geodesic_path.push_front(anc);
      anc = anc->ancestor;
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

  // we want to have the Voronoi vertices ordered in a cycle
  struct Voronoi_vertices_comparer
  {
    const std::vector<Point_2>& mapped_points;

    bool operator()(const std::size_t left, const std::size_t right)
    {
      const Point_2& left_p = mapped_points[left];
      const Point_2& right_p = mapped_points[right];

      FT left_theta = std::atan2(left_p.y(), left_p.x());
      FT right_theta = std::atan2(right_p.y(), right_p.x());

      return left_theta < right_theta;
    }

    Voronoi_vertices_comparer(const std::vector<Point_2>& mapped_points_)
      :
        mapped_points(mapped_points_)
    { }
  };

  typedef std::set<std::size_t, Voronoi_vertices_comparer>   Voronoi_vertices_container;
  typedef std::vector<Voronoi_vertices_container>            Voronoi_vertices_vector;

  void map_grid_points_to_tangent_spaces(std::vector<Point_2>& mapped_points,
                                         Voronoi_vertices_vector& Voronoi_vertices)
  {
    std::cout << "build map grids..." << std::endl;

    // transform grid points from the manifold to the tangent space of their
    // closest seed

    // conveniently, we can express it through polar geodesic coordinates
    // since we know the geodesic distance to the seed.
    // To determine the angle part of the polar coordinates, we have to compute
    // the tangent of the geodesic near the seed;

    // if the grid point is also a Voronoi vertex, we keep it in a special memory

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      const FT r = gp.distance_to_closest_seed;
      const std::size_t seed_id = gp.closest_seed_id;
      const Point_2& seed = seeds[seed_id];

      // compute the tangent and the angle with (Ox)
      FT theta = compute_tangent_angle(gp, seed);

      // we don't actually care much for polar coordinates, so we keep them as
      // Cartesian coordinates instead
      FT x = r * std::cos(theta);
      FT y = r * std::sin(theta);
      mapped_points[i] = Point_2(x, y);

//        std::cout << i << " mapped: " << gp.point.x() << " " << gp.point.y() << " to "
//                  << x << " " << y << " r/theta: " << r << " " << theta
//                  << " (seed_id at i : " << seed_id << ")" << std::endl;

      if(gp.is_Voronoi_vertex)
      {
        Voronoi_vertices[seed_id].insert(i);
      }
    }
  }

  Point_2 compute_centroid_as_closest_unmapped_grid_point(const std::size_t seed_id,
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

  Point_2 compute_centroid_with_grid_triangles(const std::size_t seed_id)
  {
    // careful now, this is for isotropic only !

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

      FT area = triangle_area(gp.point, gq.point, gr.point);

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

  Point_2 compute_centroid_with_voronoi_vertices(const std::size_t seed_id,
                                                 const Voronoi_vertices_vector& Voronoi_vertices)
  {
    // careful, this is for isotropic only !

    // in isotropic, the bisectors are segments between the voronoi vertices,
    // so we just decompose the full Voronoi cells in triangles of the form :
    // [seed, vor_vert_1, Vor_vert_2]

    const Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];
    if(Vor_vertices.size() < 3)
    {
      std::cerr << "WARNING: NOT ENOUGH VORONOI VERTICES TO OPTIMIZE SEED: " << seed_id << std::endl;
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

      const Point_2& pi = points[*it].point;
      const Point_2& pj = points[*next_it].point;

      FT area = triangle_area(old_seed, pi, pj);
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

  Point_2 compute_centroid_with_polygon_centroid(const std::size_t seed_id,
                                                 const Voronoi_vertices_vector& Voronoi_vertices)
  {
    // careful, this is for isotropic only !

    // use the formulae to compute the centroid of a polygon with the voronoi vertices

    const Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];

    if(Vor_vertices.size() < 3)
    {
      std::cerr << "WARNING: NOT ENOUGH VORONOI VERTICES TO OPTIMIZE SEED: " << seed_id << std::endl;
      return Point_2(0.,0.);
    }

    // compute the centroid of the polygon composed by the mapped Voronoi vertices
    // see wikipedia for the formulae...
    FT area = 0;
    FT centroid_x = 0., centroid_y = 0.;

    typename Voronoi_vertices_container::const_iterator it = Vor_vertices.begin();
    typename Voronoi_vertices_container::const_iterator next_it = ++(Vor_vertices.begin());
    typename Voronoi_vertices_container::const_iterator end = Vor_vertices.end();
    for(; it!=end; ++it, ++next_it)
    {
      if(next_it == end)
        next_it = Vor_vertices.begin();

      const Point_2& pi = points[*it].point;
      const Point_2& pj = points[*next_it].point;

      area += pi.x()*pj.y() - pj.x()*pi.y();
      centroid_x += (pi.x() + pj.x()) * (pi.x()*pj.y() - pj.x()*pi.y());
      centroid_y += (pi.y() + pj.y()) * (pi.x()*pj.y() - pj.x()*pi.y());
    }
    area *= 0.5;
    CGAL_assertion(area != 0);

    FT denom = 1./(6.*area);

    centroid_x *= denom;
    centroid_y *= denom;

    return Point_2(centroid_x, centroid_y);
  }

  FT optimize_seed(const std::size_t seed_id,
                   const std::vector<Point_2>& mapped_points,
                   const Voronoi_vertices_vector& Voronoi_vertices)
  {
    const Point_2 old_seed = seeds[seed_id];
    const Voronoi_vertices_container& Vor_vertices = Voronoi_vertices[seed_id];

    std::cout << "optimize seed " << seed_id << std::endl;
    std::cout << "currently at position: " << seeds[seed_id] << std::endl;
    std::cout << Vor_vertices.size() << " vor vertices" << std::endl;

    if(Vor_vertices.size() < 3)
    {
      std::cerr << "WARNING: NOT ENOUGH VORONOI VERTICES TO OPTIMIZE SEED: " << seed_id << std::endl;
      return 0.;
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

      const Point_2& pi = mapped_points[*it];
      const Point_2& pj = mapped_points[*next_it];

      std::cout << "current Vor vertex: " << *it << " " << pi;
      std::cout << " which is an angle of: " << std::atan2(mapped_points[*it].y(),
                                                           mapped_points[*it].x()) << std::endl;

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


    // -------------------------------------------------------------------------
    std::cout << "UNMAPPING POSSIBILITIES :" << std::endl;

    Point_2 closest_grid_point_centroid = compute_centroid_as_closest_unmapped_grid_point(seed_id, mapped_centroid, mapped_points);
    std::cout << "unmapped centroid as closest unmapped grid point " << closest_grid_point_centroid << std::endl;

    Point_2 bar_centroid = compute_centroid_with_barycentric_info(seed_id, mapped_centroid, mapped_points);
    std::cout << "unmapped centroid from mapped polygon (barycentric) " << bar_centroid << std::endl;

    // --
      std::cout << "FOR REFERENCE, ISOTROPIC GIVES :" << std::endl;

      // the following alt centroid computations are only relevant in isotropic
      Point_2 alt_centroid_1 = compute_centroid_with_grid_triangles(seed_id);
      std::cout << "centroid with grid triangles " << alt_centroid_1 << std::endl;

      Point_2 alt_centroid_2 = compute_centroid_with_voronoi_vertices(seed_id);
      std::cout << "centroid with voronoi vertices " << alt_centroid_2 << std::endl;

      Point_2 alt_centroid_3 = compute_centroid_with_polygon_centroid(seed_id, Voronoi_vertices);
      std::cout << "centroid with polygon centroid " << alt_centroid_3 << std::endl;
    // --

    std::cout << "ALTERNATIVE MAPPED CENTROID COMPUTATION : " << std::endl;

    Point_2 mapped_grid_centroid = compute_centroid_with_mapped_grid_triangles(seed_id, mapped_points);
    std::cout << "(mapped) centroid with mapped grid triangles " << mapped_grid_centroid << std::endl;

    Point_2 bar_centroid_2 = compute_centroid_with_barycentric_info(seed_id, mapped_grid_centroid, mapped_points);
    std::cout << "unmapped centroid from mapped grid centroid (barycentric) " << bar_centroid_2 << std::endl;

    // -------------------------------------------------------------------------

    // need seed is now the closest mapped point, unmapped back to the manifold
    const Point_2& new_seed = alt_centroid_2;
    seeds[seed_id] = new_seed;

    // squared displacement between the old and the new seed
    typename K::Compute_squared_distance_2 sqd = K().compute_squared_distance_2_object();
    FT displacement = sqd(old_seed, new_seed);
    std::cout << "seed: " << seed_id << " displacement: " << displacement << std::endl;

    // fixme
    // since we reset after computing, it's pointless to initialize now
    // need to keep it in memory for a little bit and then init after the reset...
//     initialize_grid_point(points[closest_mapped_grid_point_id], 0., seed_id);

    return displacement;
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
      std::vector<Point_2 > mapped_points(points.size());
      // the mapped Voronoi vertices for all the seeds
      Voronoi_vertices_vector Voronoi_vertices(seeds.size(),
                                     Voronoi_vertices_container(mapped_points));

      map_grid_points_to_tangent_spaces(mapped_points, Voronoi_vertices);

      FT cumulated_displacement = 0;
      for(std::size_t i=0; i<seeds.size(); ++i)
        cumulated_displacement += optimize_seed(i, mapped_points,
                                                Voronoi_vertices);

      std::cout << "at : " << counter << ", cumulated displacement : "
                << cumulated_displacement << std::endl;
      reset();
      locate_and_initialize_seeds();
      spread_distances();

      output_grid_data_and_dual("optimized_" + str_base_mesh + "_tr");

      is_optimized = (++counter > 100); // (cumulated_displacement < sq_bbox_diag_l * 1e-5);
    }
    while(!is_optimized);

  }

  // output stuff --------------------------------------------------------------

  void add_simplex_to_triangulation(const Tri& grid_tri,
                                    const std::set<std::size_t>& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    // test against criteria should be moved somewhere else than during the dual computations todo
    test_simplex(grid_tri, dual_simplex);

    add_simplex_to_triangulation(dual_simplex);
  }

  void add_simplex_to_triangulation(const Grid_point* gp,
                                    const Simplex& dual_simplex)
  {
    if(dual_simplex.size() <= 1)
      return;

    // test against criteria should be moved somewhere else than during the dual computations todo
    test_simplex(gp, dual_simplex);

    add_simplex_to_triangulation(dual_simplex);
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

  void compute_dual()
  {
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0; i<triangles.size(); ++i)
    {
      const Tri& triangle = triangles[i];
      Simplex dual_simplex;

      for(int j=0; j<3; ++j)
        dual_simplex.insert(points[triangle[j]].closest_seed_id);
      add_simplex_to_triangulation(triangle, dual_simplex);
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
    // not efficient (just keep the seeds' metric in memory instead) fixme
    const Point_2& p0 = seeds[tr[0]];
    const Point_2& p1 = seeds[tr[1]];
    const Point_2& p2 = seeds[tr[2]];
    const Metric& m0 = mf->compute_metric(p0);
    const Metric& m1 = mf->compute_metric(p1);
    const Metric& m2 = mf->compute_metric(p2);

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
#if (verbose > 5)
    std::cout << "output straight dual" << std::endl;
#endif

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
      out << is_triangle_intersected(tr, dual_edges) << std::endl;
    }
    out << "End" << std::endl;
  }

  void output_grid_data_and_dual(const std::string str_base)
  {
    output_grid(str_base);
    output_straight_dual(str_base);
  }

  Base_mesh()
    :
      points(),
      triangles(),
      trial_points(),
      edge_bisectors(),
      dual_edges(),
      dual_triangles(),
      refinement_point(NULL),
  { }
};

void initialize_seeds()
{
  vertices_nv = build_seeds();
  CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );
}

int main(int, char**)
{
  mf = new Euclidean_metric_field<K>(1., 1.);
//  mf = new Custom_metric_field<K>();

//  generate_grid();
//  exit(0);

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);
  double duration;
  start = std::clock();
  std::srand(0);

  Base_mesh bm;
  bm.initialize_base_mesh();
  initialize_seeds();
  bm.locate_and_initialize_seeds();

  bm.spread_distances();

  if(n_refine)
  {
    bm.output_grid_data_and_dual("pre_ref");

    for(int i=0; i<n_refine; ++i)
    {
      std::ostringstream out;
      out << "ref_" << i;
      bm.output_grid_data_and_dual(out.str());

      if(!bm.refine_grid_with_self_computed_ref_point())
        break;
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cerr << "End refinement: " << duration << std::endl;
  }

  bm.output_grid_data_and_dual(str_base_mesh + "_tr");
  bm.check_edelsbrunner();

  // optimize stuff
  bm.optimize_seeds();
  bm.output_grid_data_and_dual("optimized_" + str_base_mesh + "_tr");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
