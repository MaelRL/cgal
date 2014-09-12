#ifndef CGAL_ANISOTROPIC_MESH_2_CONSTRAIN_SURFACE_2_IMPLICIT_H
#define CGAL_ANISOTROPIC_MESH_2_CONSTRAIN_SURFACE_2_IMPLICIT_H

#define ZERO_DOT_PRODUCT 1e-5
#define ZERO_EIGENVALUE 1e-10

#include <CGAL/Constrain_surface_2.h>

#include <Eigen/Dense>

#include <CGAL/Implicit_mesh_domain_2.h>
#include <CGAL/make_mesh_2.h>

#include <CGAL/gl_draw/drawing_helper.h>

#include <CGAL/helpers/metric_helper.h>
#include <CGAL/helpers/c3t3_polyhedron_builder.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/bind.hpp>

namespace CGAL
{
namespace Anisotropic_mesh_2
{
template<typename PointerToMemberFunction, typename PointerToObject, typename FT, typename Point>
class Member_function_pointer_to_function_wrapper
{
public:
  typedef FT return_type;

  /// Constructor
  Member_function_pointer_to_function_wrapper(PointerToObject po, PointerToMemberFunction pf)
  : pf_(pf) , po_(po) {}

  // Default copy constructor and assignment operator are ok
  /// Destructor
  ~Member_function_pointer_to_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point& p) const
  {
    return ((*po_).*(pf_))(p);
  }
private:
  /// Function to wrap
  PointerToMemberFunction pf_;
  PointerToObject po_;
};  // end class Member_function_pointer_to_function_wrapper


template<typename K,
         typename Pt_container = std::vector<typename K::Point_2> >
class Constrain_surface_2_implicit :
  public Constrain_surface_2<K, Pt_container>
{
public:
  typedef Constrain_surface_2_implicit<K, Pt_container> Self;
  typedef Constrain_surface_2<K, Pt_container>          Base;
  typedef typename Base::Object_2                       Object_2;
  typedef typename Base::FT                             FT;
  typedef typename Base::Vector_2                       Vector_2;
  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::Plane_2                        Plane_2;
  typedef typename Base::Oriented_side                  Oriented_side;
  typedef typename Base::Colored_poly                   Colored_polyhedron;

  typedef Pt_container                                  Point_container;

  // for Mesh_2
  typedef FT (Self::*Function)(const Point_2& p) const;

  typedef Member_function_pointer_to_function_wrapper<Function, const Self*, FT, Point_2> Function_wrapper;
  typedef typename CGAL::Implicit_mesh_domain_2<Function_wrapper, K> Mesh_domain;

  // Triangulation
  typedef typename CGAL::Mesh_triangulation_2<Mesh_domain>::type Tr;
  typedef typename CGAL::Mesh_complex_2_in_triangulation_2<Tr> C3t3;

  // Criteria
  typedef typename CGAL::Mesh_criteria_2<Tr> Mesh_criteria;

protected:
  mutable C3t3 m_c3t3;
  FT error_bound;

public:
  C3t3 c3t3() const { return m_c3t3; }

  virtual FT get_bounding_radius() const = 0;
  virtual typename CGAL::Bbox_2 get_bbox() const = 0;
  virtual std::string name() const { return std::string("Implicit"); }

  virtual FT evaluate(const FT x, const FT y, const FT z) const = 0;
  FT implicit_function(const Point_2& p) const
  {
    return evaluate(p.x(), p.y());
  }

  virtual Oriented_side side_of_constraint(const Point_2 &p) const
  {
    FT w = evaluate(p.x(), p.y(), p.z());
    if (w < 0)
      return CGAL::ON_POSITIVE_SIDE;
    else if (w > 0)
      return CGAL::ON_NEGATIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

public:
  // Sizing field
  template<typename CSI>
  struct Aniso_sizing_field
  {
    const CSI* const csi;

    typedef typename CSI::FT                 FT;
    typedef typename CSI::Mesh_domain::Index Index;
    FT operator()(const Point_2& p, const int, const Index&) const
    {
      Vector_2 v0, v1, v2;
      double e0, e1, e2;
      csi->tensor_frame(p, v0, v1, v2, e0, e1, e2, 1e-5);
      e1 = std::max(e1, 1e-10);
      e2 = std::max(e2, 1e-10);
      FT ratio = 1./std::sqrt(std::max(e1/e2, e2/e1));
      FT ret = std::max(1e-3, ratio);

      ret *= 0.1*csi->get_bounding_radius();

      //std::cout << "ratio : " << ratio << " " << 1./ratio << std::endl;
      //std::cout << "ret : " << ret << std::endl;
      return ret;
    }

    Aniso_sizing_field(const CSI* const csi_):csi(csi_){}
  };

  virtual Point_container get_surface_points(unsigned int nb,
                                             double facet_distance_coeff /*= 0.05*/) const
  {
    std::vector<Point_2> all_points;
    typename C3t3::Triangulation::Finite_vertices_iterator v;
    for(v = m_c3t3.triangulation().finite_vertices_begin();
        v != m_c3t3.triangulation().finite_vertices_end();
        ++v)
    {
      if(m_c3t3.in_dimension(v) == 1 || m_c3t3.in_dimension(v) == 2)
        all_points.push_back(v->point());
    }

    if(std::size_t(nb) > all_points.size())
    {
#ifdef ANISO_VERBOSE
      std::cout << "C3T3 mesh did not contain enough points (" << all_points.size();
      std::cout << "). Generating a new c3t3...";
      std::cout << "Facet distance : " << facet_distance_coeff << std::endl;
#endif
      //not enough points on the c3t3 to fill the initial_points vector. Generating a new denser c3t3.
      Function fct = (Function)(&Self::implicit_function);
      FT r = this->get_bounding_radius();
      Function_wrapper fw(this, fct);
      Mesh_domain domain(fw, typename K::Sphere_2(CGAL::ORIGIN, r*r), error_bound);

      Aniso_sizing_field<Self> size(this);
      Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.,
                             CGAL::parameters::facet_size = r * 0.05, //size
                             CGAL::parameters::facet_distance = r * facet_distance_coeff,
                             CGAL::parameters::facet_topology = MANIFOLD);
                             // cell criteria are ignored

      // run Mesh_2
      m_c3t3 = CGAL::make_mesh_2<C3t3>(domain, criteria,
                                       CGAL::parameters::no_perturb(),
                                       CGAL::parameters::no_exude());
#if 1//def ANISO_OUTPUT_MESH_FOR_POLES
      std::ofstream out("mesh_2_temp.mesh");
      m_c3t3.output_to_medit(out);
#endif
      return get_surface_points(nb, facet_distance_coeff/3.);
    }
#ifdef ANISO_VERBOSE
    std::cout << "c3t3 has enough points : " << all_points.size() << " faces : " << m_c3t3.number_of_facets_in_complex() << std::endl;
#endif
    std::random_shuffle(all_points.begin(), all_points.end());
    return Point_container(all_points.begin(), (all_points.begin() + nb));
  }

  C3t3 run_mesh_2(const FT& approx, const bool verbose = false) const
  {
    Function fct = (Function)(&Self::implicit_function);
    FT r = this->get_bounding_radius();
    Function_wrapper fw(this, fct);
    Mesh_domain domain(fw, typename K::Sphere_2(CGAL::ORIGIN, r*r), error_bound);
    Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.,
                           CGAL::parameters::facet_size = r * 0.05,
                           CGAL::parameters::facet_distance = approx,
                           CGAL::parameters::facet_topology = MANIFOLD);
                           //no cell criteria
    C3t3 mesh = CGAL::make_mesh_2<C3t3>(domain, criteria,
                                   CGAL::parameters::no_perturb(),
                                   CGAL::parameters::no_exude());
    if(verbose)
      std::cout << "Mesh has "
                << mesh.triangulation().number_of_vertices () << " vertices.\n";
    return mesh;
  }

  virtual void compute_poles(std::set<Point_2>& poles) const
  {
    compute_triangulation_poles(m_c3t3, std::inserter(poles, poles.end()), get_bbox());
  }

  void gl_draw_intermediate_mesh_2(const Plane_2& plane) const
  {
    gl_draw_c3t3<C3t3, Plane_2>(m_c3t3, plane);
  }

  virtual Constrain_surface_2_implicit* clone() const = 0; // Uses the copy constructor

  Constrain_surface_2_implicit(const FT error = 1e-5)
    :
      Base(),
      error_bound(error)
  { }

  virtual ~Constrain_surface_2_implicit(){}
};

template<typename K>
class Constrain_surface_2_implicit_with_bounding_sphere :
  public Constrain_surface_2_implicit<K>
{

public:
  typedef Constrain_surface_2_implicit<K>                  Base;
  typedef typename Base::FT                                FT;
  typedef typename Base::Vector_2                          Vector_2;
  typedef typename Base::Point_2                           Point_2;
  typedef typename Base::Oriented_side                     Oriented_side;

  typedef typename std::vector<typename K::Point_2>        Point_container;

public:
  virtual FT get_bounding_radius() const = 0;
  virtual FT evaluate(const FT x, const FT y, const FT z) const = 0;

  virtual Oriented_side side_of_constraint(const Point_2 &p) const
  {
    FT br = get_bounding_radius();
    if (p.x() * p.x() + p.y() * p.y() + p.z() * p.z() > br * br)
      return CGAL::ON_NEGATIVE_SIDE;
    FT w = evaluate(p.x(), p.y(), p.z());
    if (w < 0)
      return CGAL::ON_POSITIVE_SIDE;
    else if (w > 0)
      return CGAL::ON_NEGATIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Constrain_surface_2_implicit_with_bounding_sphere() { }
  ~Constrain_surface_2_implicit_with_bounding_sphere() { }
};

}
}

#endif
