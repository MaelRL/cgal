#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_cyclide_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_cyclide_H

#include <CGAL/Constrain_surface_3_implicit.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Constrain_surface_3_cyclide : public Constrain_surface_3_implicit<K>
{
public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT                          FT;
  typedef typename Base::Point_3                     Point_3;
  typedef typename Base::Oriented_side               Oriented_side;
  typedef typename Base::Point_container             Point_container;

public:
  FT a, b, c, mu;

public:
  void set_a(const FT& aa) { a = aa; }
  void set_b(const FT& bb) { b = bb; }
  void set_mu(const FT& mmu) { mu = mmu; }
  FT get_a() const { return a; }
  FT get_b() const { return b; }
  FT get_mu() const { return mu; }
  
  virtual std::string name() const
  {
    std::ostringstream o;
    o << "Implicit Cyclide (" << a << ", " << b << ", " << mu << ")";
    return o.str();
  }

  virtual FT get_bounding_radius() const
  {
    return 15;
  }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    typedef typename CGAL::Bbox_3 Bbox;
    return Bbox(-10, -10, -10,
                10, 10, 10);
  }

  virtual FT evaluate(const FT x, const FT y, const FT z) const
  {
      return (x*x + y*y + z*z - mu*mu + b*b)*(x*x + y*y + z*z - mu*mu + b*b)
            - 4*(a*x - c*mu)*(a*x - c*mu) - 4*b*b*y*y;
  }

  Oriented_side side_of_constraint(const Point_3 &p) const
  {
    FT w = evaluate(p.x(), p.y(), p.z());
    if (w < 0)
      return CGAL::ON_POSITIVE_SIDE;
    else if (w > 0)
      return CGAL::ON_NEGATIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Point_3 point_on_surface(const FT& u/*small circle*/,
                           const FT& v/*big circle*/) const
  {
    FT denom = a - c*std::cos(v)*std::cos(u);

    FT num_x = mu*(c-a*std::cos(v)*std::cos(u))+b*b*std::cos(v);
    FT num_y = b*std::sin(v)*(a-mu*std::cos(u));
    FT num_z = b*std::sin(u)*(c*std::cos(v)-mu);

    FT x = num_x / denom;
    FT y = num_y / denom;
    FT z = num_z / denom;
    return Point_3(x,y,z);
  }

  Point_container initial_points(const int nb = 8) const
  {
    Point_container points;
    std::vector<Point_3> seeds;
    int nb_per_small_circle = 10; // 5 vertices are enough for one small circle
    int nb_per_large_circle = 100; //(int)std::ceil(nb / (nb_per_small_circle + 0.));

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
    return Base::initial_points(points,seeds, 0.05, nb);
  }

  Constrain_surface_3_cyclide* clone() const // Covariant Return Types
  { 
    return new Constrain_surface_3_cyclide(*this);
  }

  Constrain_surface_3_cyclide(FT a_ = 3., FT b_ = 1., FT mu_ = 2.)
  : a(a_), b(b_), mu(mu_), c(std::sqrt(a*a - b*b))
  { std::cout << "c is : " << c << std::endl; }

  Constrain_surface_3_cyclide(const Constrain_surface_3_cyclide& cy)
  : a(cy.a), b(cy.b), mu(cy.mu), c(cy.c)
  {}

  ~Constrain_surface_3_cyclide() { }
};

} //Anisotropic_mesh_3
} //CGAL

#endif
