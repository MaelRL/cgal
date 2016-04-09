#ifndef CGAL_ANISOTROPIC_MESH_2_CUSTOM_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_2_CUSTOM_METRIC_FIELD_H

#include <CGAL/Metric_field.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K>
class Custom_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>           Base;
  typedef typename Base::FT         FT;
  typedef typename Base::Metric     Metric;
  typedef typename Base::Point_2    Point_2;
  typedef typename Base::Vector_2   Vector_2;

public:
  Metric radial_shock(const Point_2& p) const
  {
    FT x = p.x();
    FT y = p.y();

    if(x == 0 && y == 0)
      return Metric();

    FT r = x*x+y*y;
    FT nx = x/std::sqrt(r);
    FT ny = y/std::sqrt(r);

    Vector_2 v1(-ny, nx);
    Vector_2 v2(nx, ny);

    FT s = std::sin(r/5.);
    FT c = std::cos(r/5.);

    FT e1 = (-5*s-2*r*(c-3*s*s))/std::pow(r*s*s,1./5.) ;
    FT e2 = s/std::sqrt(r*s*s);

    return this->build_metric(v1, v2, std::sqrt(e1), std::sqrt(e2));
  }

  Metric starred_shock(const Point_2& p) const
  {
    if(p == CGAL::ORIGIN)
      return Metric();

    FT h = 0.05;
    FT phi = 1.;

    FT x = p.x();
    FT y = p.y();

    FT g1 = 3*x*x - x*y/5. -y*y + 3*x*x/10.;
    FT g2 = -x*x/10. - 2*x*y;

    FT n = std::sqrt(g1*g1 + g2*g2);

    FT l1 = (std::max)(0.1, 1./(1. + phi*n));
    FT l2 = 1.;

    Vector_2 v1 = Vector_2(g1/n, g2/n);
    Vector_2 v2 = Vector_2(-g2/n, g1/n);

    // probably need a sqrt on l1 and l2 (but not h)
    return Metric(v1, v2, l1/h, l2/h, this->epsilon);
  }

  Metric hyperbolic_shock(const Point_2 &p,
                          const FT delta = 0.6) const
  {
    FT h = 0.1; // 0.05

    FT x = p.x();
    FT y = p.y();
    FT tanhder = tanh((2.0 * x - sin(5.0 * y)) / delta);
    tanhder = (1.0 - tanhder * tanhder) / delta;
    FT x1 = 2. * tanhder + 3.0 * x * x + y * y;
    FT y1 = -tanhder * cos(5.0 * y) * 5.0 + 2.0 * x * y;

    FT r = sqrt(x1 * x1 + y1 * y1);
    FT x2 = -y1 / r;
    FT y2 = x1 / r;
    FT l1 = std::sqrt(x1*x1 + y1*y1);
    FT l2 = std::sqrt(x2*x2 + y2*y2);

    CGAL_assertion(l1 != 0. && l2 != 0.);

    Vector_2 v1 = (1./l1) * Vector_2(x1, y1);
    Vector_2 v2 = (1./l2) * Vector_2(x2, y2);
    return Metric(v1, v2, l1 / h, l2 / h, this->epsilon);
  }

  Metric yang_liu_cube_shock_1D(const Point_2& p) const
  {
    double x = p.x();
    double y = p.y();

    double h1 = 1./(0.0025+0.2*(1-std::exp(-std::abs(x-0.6))));
    double h2 = 5.;
    h2 = 1./(0.0025+0.2*(1-std::exp(-std::abs(y-0.3))));

    Vector_2 v1(1.,0.);
    Vector_2 v2(0.,1.);

    return this->build_metric(v1, v2, h1, h2);
  }

  Metric yang_liu_cube_shock(const Point_2& p) const
  {
    double h = 0.3;

    double x = p.x();
    double y = p.y();
    double r = std::sqrt(x*x + y*y);

    if(r<0.1)
      return this->build_metric(Vector_2(1.,0.), Vector_2(0.,1.),
                                1./h, 1./h);

    x /= r;
    y /= r;

    double lambda = std::exp(-0.5*std::abs(r*r-0.5));
//    double h1 = 0.025 + (1-std::exp(-0.01 * (r*r-4)*(r*r-4) )); // to get something smooth
    double h1 = 0.1 + (1-lambda); // artificial smoothing
//    double h1 = 0.025 + (1-lambda); // original
    double h2 = 1.;

    //lambda = std::exp(-0.01*std::abs(r*r-49));
    //h1 /= 10*lambda; h2 /= 10*lambda;

    Vector_2 v1(x, y);
    Vector_2 v2(-y, x);

    return this->build_metric(v1, v2, 1./(h*h1), 1./(h*h2));
  }

  Metric cos_1D(const Point_2& p) const
  {
    double h = 0.3;
//    Vector_2 v1(1, 0.);
//    Vector_2 v2(0., 1.);
    Vector_2 v1(p.x(), p.y());
    Vector_2 v2(p.y(), -p.x());

    double h1 = (std::max)(1e-5, std::abs(std::cos(3*p.x())));
    double h2 = 1.;

    return this->build_metric(v1, v2, 1./(h*h1), 1./(h*h2));
  }

  virtual Metric compute_metric(const Point_2 &p) const
  {
    return starred_shock(p);
    return radial_shock(p);
    return hyperbolic_shock(p);
    return yang_liu_cube_shock(p);
    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              std::sqrt(6.*p.x()*p.x() + 1.),
                              std::sqrt(6.*p.y()*p.y() + 1.));

    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              std::sqrt(6.*p.y()*p.y() + 1.),
                              std::sqrt(6.*p.x()*p.x() + 1.));

    return cos_1D(p);
    return yang_liu_cube_shock_1D(p);

    FT denom = std::sqrt(9.); // = ratio at x=+/-1
    double h = 1.0;
    double x = p.x();
    double y = p.y();

    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              1.1+std::exp(std::exp(x)), 1.);

    return this->build_metric(Vector_2(1, 0), // for DU WANG, Transf2
                              Vector_2(0, 1),
                              std::sqrt(1./(x*x*x*x) + 4./(x*x*x*x*x*x)),
                              std::sqrt(9.*y*y*y*y+4.*y*y));
    return this->build_metric(Vector_2(1, 0), // for TC
                              Vector_2(0, 1),
                              1./(h*x), y/(h));
    return this->build_metric(Vector_2(1, 0), // for DU WANG, Transf1
                              Vector_2(0, 1),
                              std::sqrt(1./(x*x*x*x) + 1./(x*x*x*x*x*x)),
                              std::sqrt(9.*y*y*y*y+y*y));
    return this->build_metric(Vector_2(1, 0), // for DU WANG, Transf3s
                              Vector_2(0, 1),
                              std::sqrt(1. + 4./(x*x*x*x*x*x)),
                              std::sqrt(1.+4.*std::exp(4*y)));

    //linear interpolation
    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              1./((denom*std::abs(x)+(1-std::abs(x)))*h), 1./h);

    //logexp 1D
    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              std::pow(1./denom, std::abs(x))/h, 1./h);

    //logexp 1D extended into uniform
    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              std::max(std::pow(1./denom, std::abs(x)),1./denom), 1.);
  }

  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: custom" << std::endl;
  }

public:
  Custom_metric_field(FT epsilon_ = 1e-6) : Metric_field<K>(epsilon_) { }
};

} //Aniso
} //CGAL

#endif
