#ifndef CGAL_ANISOTROPIC_MESH_3_CUSTOM_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_CUSTOM_METRIC_FIELD_H

#include <CGAL/Metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Custom_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>           Base;
  typedef typename Base::FT         FT;
  typedef typename Base::Metric     Metric;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::Vector_3   Vector_3;

public:
  Metric yang_liu_sphere_sin(const Point_3& p) const
  {
    double x = p.x();

    Vector_3 v1(2*std::cos(6*x), 1., 0.);
    v1 = v1 / std::sqrt(v1*v1);
    Vector_3 v2(0.,0.,1.);
    Vector_3 v3 = CGAL::cross_product(v1, v2);

    double h1 = std::sqrt(1000);
    double h2 = std::sqrt(10);
    double h3 = std::sqrt(10);

    return this->build_metric(v1, v2, v3, h1, h2, h3);
  }

  Metric yang_liu_cube_shock_1D(const Point_3& p) const
  {
    double h = 0.1;

    double x = p.x();
    double y = p.y();
    double z = p.z();

    double lambda_1 = std::exp(-std::abs(x - 0.6));
    double h1 = 0.0025 + (1 - lambda_1);

    double lambda_2 = std::exp(-std::abs(y - 0.6));
    double h2 = 0.0025 + (1 - lambda_2);

    double lambda_3 = std::exp(-std::abs(z - 0.6));
    double h3 = 0.0025 + (1 - lambda_3);

    Vector_3 v1(1.,0.,0.);
    Vector_3 v2(0.,1.,0.);
    Vector_3 v3(0.,0.,1.);

    return this->build_metric(v1, v2, v3,
                              1./(h*h1), 1./(h*h2), 1./(h*h3));
  }

  Metric yang_liu_bumpy(const Point_3& p) const
  {
    double h = 0.05;

    double x = p.x();
    double y = p.y();
    double z = p.z();

    double r = std::sqrt(x*x+y*y);

    Vector_3 v1(x/r, y/r, 0.);
    Vector_3 v2(-y/r, x/r, 0.);
    Vector_3 v3(0., 0., 1.);

    double lambda_1 = std::exp(-0.05*(r*r - 2.56));
    double h1 = 0.5 + 1 - lambda_1;
    double h2 = h1 * (1 + 5*r);
    double h3 = h2;

    return this->build_metric(v1, v2, v3,
                              1./(h*h1), 1./(h*h2), 1./(h*h3));
  }

  Metric yang_liu_cube_shock(const Point_3& p) const
  {
    double h = 0.05; // 0.3 for the 1.5m canvas

    //Yang Liu 3D Metric
    double x = p.x() + 0.1; // to avoid the origin...
    double y = p.y() + 0.1;
    double z = p.z() + 0.1;

    double r = std::sqrt(x*x + y*y + z*z);

//    std::cout << "r: " << r << std::endl;

    if(r < 0.05)
      return this->build_metric(Vector_3(1.,0.,0.),
                                Vector_3(0.,1.,0.),
                                Vector_3(0.,0.,1.),
                                1./h, 1./h, 1./h);

    x /= r;
    y /= r;
    z /= r;

    double lambda = std::exp(-0.5 * std::abs(r*r - 0.15));
//    double h1 = 0.025 + (1-lambda);
    double h1 = 0.5 + (1 - lambda);
    double h2 = 1.;
    double h3 = 1.;

//    std::cout << "at: " << p << " r: " << r << " || h1: " << h1 << std::endl;

    CGAL_assertion(std::abs(y) > 1e-15);

    double v2y = -y*z / std::sqrt((x*x + y*y + z*z)*(x*x + y*y));
    double v2x = x*v2y/y;
    double v2z = std::sqrt(1-(x*x/(y*y)+1)*v2y*v2y);

    Vector_3 v1(x, y, z);
    Vector_3 v2(v2x, v2y, v2z);
    Vector_3 v3 = CGAL::cross_product(v1, v2);


    return this->build_metric(v1, v2, v3, 1./(h*h1), 1./(h*h2), 1./(h*h3));
  }

  virtual Metric compute_metric(const Point_3& p) const
  {
    return yang_liu_cube_shock(p);
    return yang_liu_cube_shock_1D(p);
    return yang_liu_sphere_sin(p);

    FT denom = std::sqrt(9.); // = ratio at x=+/-1
    double h = 1.0;

    //linear interpolation
    /*
    return this->build_metric(Vector_3(1, 0, 0),
                              Vector_3(0, 1, 0),
                              Vector_3(0, 0, 1),
                              1./((denom*std::abs(p.x())+(1-std::abs(p.x())))*h), 1./h, 1./h);
    */

    //logexp 1D
    return this->build_metric(Vector_3(1, 0, 0),
                              Vector_3(0, 1, 0),
                              Vector_3(0, 0, 1),
                              std::pow(1./denom, std::abs(p.x()))/h, 1./h, 1./h);

    //logexp 1D extended into uniform
    return this->build_metric(Vector_3(1, 0, 0),
                              Vector_3(0, 1, 0),
                              Vector_3(0, 0, 1),
                              std::max(std::pow(1./denom, std::abs(p.x())),1./denom), 1., 1.);

/*
    //Various interpolations
    FT l = 10;
    FT lambda = 1. - std::abs(p.x()) / l;
    FT at0 = 1.; //value at lambda = 0
    FT at1 = 1. / l; //value at lambda = 1
    FT a  = at0 + (at1-at0)*lambda*lambda*(3.-2.*lambda); //smoothstep

    //other interpolation: ev =  l^(-[2*]lambda) to have delta(x,y) == delta(y, z) if d(x,y) == d(y,z)
    a = std::pow(l, lambda);
    a = 1./(a);

    return this->build_metric(Vector_3(1, 0, 0),
                              Vector_3(0, 1, 0),
                              Vector_3(0, 0, 1),
                              a, 1., 1.);

    FT x = 39.9667;
    l = x/2.;
    FT half_l = l/2.;
    FT max_ev = (1./half_l)*(1./half_l); // sqrt(max_ev) = desired_max_ratio = 2/l

    Point_3 pOy(0., p.y(), 0.); //projection on Oy
    FT sq_d = CGAL::squared_distance(p, pOy);

    FT ez = 1., ey = 1.;
    Vector_3 vy(0., 1., 0.);
    Vector_3 vx(pOy, p);

    if(sq_d == 0.)
      vx = Vector_3(1.,0.,0.);
    else
      vx = std::sqrt(1./(vx*vx)) * vx;

    Vector_3 vz = CGAL::cross_product(vx, vy);

    FT ex;
    if(sq_d < half_l*half_l)
      ex = max_ev;
    else if(sq_d > (l+half_l)*(l+half_l))
      ex = 1.;
    else
    {
      FT lambda = (std::sqrt(sq_d)-l/2.)/l;

      //linear interpolation
      //ex = lambda + (1-lambda)*max_ev;

      //log exp interpolation
      ex = std::exp((1-lambda)*std::log(max_ev)); // if d==l/2., lambda = 0, ex = (2./l)²
                                                    // if d==3l/2, lambda = 1, ex = 1
      //smoothstep interpolation
      FT aa = max_ev; //value at lambda = 0
      FT bb = 1.; //value at lambda = 1
      ex = aa + (bb-aa)*lambda*lambda*(3.-2.*lambda);
    }

    ey = 0.005;

    Metric M = this->build_metric(vx, vy, vz, ex, ey, ez);
    return M;
*/
  }

  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: custom" << std::endl;
  }

public:
  Custom_metric_field(FT epsilon_ = 1e-6, FT en_factor = 1.)
  : Metric_field<K>(epsilon_, en_factor) { }
};

} //Aniso
} //CGAL

#endif
