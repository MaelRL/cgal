#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_DOUBLE_ELLIPSOID_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_DOUBLE_ELLIPSOID_H

#include <CGAL/Constrain_surface_3_implicit.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Constrain_surface_3_double_ellipsoid : public Constrain_surface_3_implicit<K>
{
public:
  typedef Constrain_surface_3_implicit<K>  Base;
  typedef typename Base::FT                FT;
  typedef typename Base::Point_3           Point_3;
  typedef typename Base::Point_container   Point_container;

public:
  FT a, b, c; //first ellipsoid
  FT d, e, f; //second ellipsoid
  Point_3 s_center; //second ellipsoid center
  FT bounding_radius;

public:
  void set_a(const FT& aa) { a = aa; }
  void set_b(const FT& bb) { b = bb; }
  void set_c(const FT& cc) { c = cc; }
  void set_d(const FT& dd) { c = dd; }
  void set_e(const FT& ee) { c = ee; }
  void set_f(const FT& ff) { c = ff; }
  void set_s_center(const Point_3& p) { s_center = p; }

  FT get_a() const { return a; }
  FT get_b() const { return b; }
  FT get_c() const { return c; }
  FT get_d() const { return d; }
  FT get_e() const { return e; }
  FT get_f() const { return f; }
  Point_3 get_s_center() const { return s_center; }

  virtual std::string name() const { return std::string("Implicit double-ellipsoid"); }

  FT get_bounding_radius() const { return bounding_radius; }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    // a bit brutal
    FT aa = 6.; //std::max(a,d) + std::abs(s_center.x());
    FT bb = 3.; //std::max(b,e) + std::abs(s_center.y());
    FT cc = 4.; //std::max(c,f) + std::abs(s_center.z());
    aa *= 1.1; bb *= 1.1; cc *= 1.1;

    return CGAL::Bbox_3(-aa, -bb, -cc, aa, bb, cc);
  }

  FT evaluate(const FT x, const FT y, const FT z) const 
  {
    FT xx = x - s_center.x();
    FT yy = y - s_center.y();
    FT zz = z - s_center.z();
    FT first_elli = x*x/(a*a) + y*y/(b*b) + z*z/(c*c) - 1.0;
    FT second_elli = xx*xx/(d*d) + yy*yy/(e*e) + zz*zz/(f*f) - 1.0;
    return  1.0 - std::exp(-0.7*first_elli) - std::exp(-0.7*second_elli);
  }

//  virtual double global_max_curvature() const
//  {
//    std::cout << "using bad curv values" << std::endl;
//    return 1e30; //todo
//  }
//  virtual double global_min_curvature() const
//  {
//    std::cout << "using bad curv values" << std::endl;
//    return -1e30; //todo
//  }

  Point_container initial_points(const int nb = 8) const 
  {
    Point_container points;
    std::vector<Point_3> seeds;
    seeds.push_back(Point_3(0, 0, 0));
    return Base::initial_points(points, seeds, 0.2);
  }

  Constrain_surface_3_double_ellipsoid* clone() const //covariant return types
  {
    return new Constrain_surface_3_double_ellipsoid(*this);
  }

  Constrain_surface_3_double_ellipsoid(const FT a_ = 1.5, const FT b_ = 2.5, const FT c_ = 3.5,
                                       const FT d_ = 1.5, const FT e_ = 2.5, const FT f_ = 3.5,
                                       const Point_3 P = Point_3(3.5,0.0,0.0)) :
      a(a_), b(b_), c(c_), d(d_), e(e_), f(f_), s_center(P)
  {
    FT dist_to_s_center = std::abs(s_center.x()*s_center.x() + s_center.y()*s_center.y() + s_center.z()*s_center.z());
    FT max_elli1 = (std::max)((std::max)(a_, b_), c_);
    FT max_elli2 = (std::max)((std::max)(d_, e_), f_);
    bounding_radius = 6.0; //(dist_to_s_center + max_elli1 + max_elli2) * 1.1;
  }

  Constrain_surface_3_double_ellipsoid(const Constrain_surface_3_double_ellipsoid& de)
      : a(de.get_a()), b(de.get_b()), c(de.get_c()),
        d(de.get_d()), e(de.get_e()), f(de.get_f()),
        s_center((de.get_s_center()))
  {
    FT dist_to_s_center = std::abs(s_center.x()*s_center.x() + s_center.y()*s_center.y() + s_center.z()*s_center.z());
    FT max_elli1 = (std::max)((std::max)(a, b), c);
    FT max_elli2 = (std::max)((std::max)(d, e), f);
    bounding_radius = 6.0; //(dist_to_s_center + max_elli1 + max_elli2) * 1.1;
  }

  ~Constrain_surface_3_double_ellipsoid() { }
};

} //Anisotropic_mesh_3
} //CGAL

#endif
