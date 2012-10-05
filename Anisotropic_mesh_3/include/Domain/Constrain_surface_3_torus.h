

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_TORUS_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_TORUS_H

#include <CGAL/Constrain_surface_3_implicit.h>
#include <math.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Constrain_surface_3_torus : public Constrain_surface_3_implicit<K> 
{
public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT                          FT;
  typedef typename Base::Point_3                     Point_3;
  typedef typename Base::Oriented_side               Oriented_side;
  typedef typename Base::Point_container             Point_container;

public:
  FT R, r;
public:

  void set_r(const FT& rr) { r = rr; }
  void set_R(const FT& rr) { R = rr; }
  FT get_r() { return r; }
  FT get_R() { return R; }
  
  virtual std::string name() const
  {
    std::ostringstream o;
    o << "Implicit Torus (" << r << ", " << R << ")";
    return o.str();
  }

  virtual FT get_bounding_radius() const { return (R + r) * 1.01; }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    typedef typename CGAL::Bbox_3 Bbox; 
    return Bbox(-(R+r)*1.01, -(R+r)*1.01, -r*1.01, 
                 (R+r)*1.01,  (R+r)*1.01,  r*1.01);
  }

  virtual FT evaluate(const FT x, const FT y, const FT z) const
  {
    FT d = R - std::sqrt(x * x + y * y);
    return d*d + z*z - r*r;
  }

  Oriented_side side_of_constraint(const Point_3 &p) const {
    FT x = p.x(), y = p.y(), z = p.z();
    FT l2 = x * x + y * y;
    if (l2 < (R - r) * (R - r))
      return CGAL::ON_NEGATIVE_SIDE;
    FT l = R / sqrt(l2);
    FT val = (x * l - x) * (x * l - x) + (y * l - y) * (y * l - y) + z * z;
    if (val < r * r)
      return CGAL::ON_POSITIVE_SIDE;
    else if (val > r * r)
      return CGAL::ON_NEGATIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Point_container initial_points(const int nb = 40) const 
  {
    Point_container points;
    std::vector<Point_3> seeds;
    int nb_per_small_circle = 5; // 5 vertices are enough for one small circle
    int nb_per_large_circle = (int)std::ceil(nb / (nb_per_small_circle + 0.)); 

    CGAL::Random random;
    FT angle_i = (2. * CGAL_PI / nb_per_small_circle);
    FT angle_j = (2. * CGAL_PI / nb_per_large_circle);

    for(int j = 0; j < nb_per_large_circle; j++)
      for(int i = 0; i < nb_per_small_circle; i++)
      {
        double ri = random.get_double(-0.2, 0.2);
        double rj = random.get_double(-0.3, 0.3);
        Point_3 p = point_on_surface((i+ri) * angle_i,
          (j+rj) * angle_j);
        points.push_back(p);
        seeds.push_back(p);
      }

      //seeds.push_back(Point_3(R, 0, 0));
      //seeds.push_back(Point_3(0, R, 0));
      //seeds.push_back(Point_3(-R, 0, 0));
      //seeds.push_back(Point_3(0, -R, 0));
      return Base::initial_points(points,seeds, 0.05, nb);
  }

  virtual double global_max_curvature() const
  {
    return 1./r;
  }
  virtual double global_min_curvature() const
  {
    return 1./(R+r);
  }


  Point_3 point_on_surface(const FT& u/*small circle*/, 
                           const FT& v/*big circle*/) const
  {
    FT x = (R + r * std::cos(u)) * std::cos(v);
    FT y = (R + r * std::cos(u)) * std::sin(v);
    FT z = r * std::sin(u);
    return Point_3(x,y,z);
  }

  Constrain_surface_3_torus* clone() const // Covariant Return Types
  { 
    return new Constrain_surface_3_torus(*this); 
  }

  Constrain_surface_3_torus(FT R_ = 0.7, FT r_ = 0.3) 
    : R(R_), r(r_) {}

  Constrain_surface_3_torus(const Constrain_surface_3_torus& t)
    : R(t.R), r(t.r) {}
  
  ~Constrain_surface_3_torus() { };
};

#endif
