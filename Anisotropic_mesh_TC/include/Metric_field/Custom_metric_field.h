#ifndef CGAL_ANISOTROPIC_MESH_TC_CUSTOM_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_TC_CUSTOM_METRIC_FIELD_H

#include <CGAL/Metric_field.h>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Kd>
class Custom_metric_field : public Metric_field<Kd>
{
public:
  typedef Metric_field<Kd>                            Base;
  typedef typename Kd::Dimension                      Dim;

  typedef typename Base::FT                           FT;
  typedef typename Base::Metric                       Metric;
  typedef typename Base::Point_d                      Point_d;
  typedef typename Base::Vector_d                     Vector_d;
  typedef typename Base::Metric::E_Vector             E_Vector;

public:
  Metric shock_1D(const Point_d& p) const
  {
    int d = Dim::value;
    typename Kd::Construct_vector_d constr_vec = Kd().construct_vector_d_object();
    int dir_id = 0; // Ox axis

    std::vector<FT> evals(d);
    std::vector<Vector_d> evecs;
    for(int i=0; i<d; ++i)
    {
      if(i == dir_id)
        evals[i] = 1./(0.0025+0.2*(1-std::exp(-std::abs(p[i]-0.6))));
      else
        evals[i] = 5.;


      E_Vector vi = E_Vector::Zero();
      vi(i) = 1.;
      evecs.push_back(constr_vec(d, vi.data(), vi.data() + d));
    }

    evals[1] = 1./(0.0025+0.2*(1-std::exp(-std::abs(p[1]-0.3)))); // tmp
    return this->build_metric(evecs, evals);
  }

  Metric yang_liu_cube_shock(const Point_d& p) const // fixme this is only 2D
  {
    int d = Dim::value;
    typename Kd::Construct_vector_d constr_vec = Kd().construct_vector_d_object();
    double h = 0.3;
    std::vector<FT> evals(d);
    std::vector<Vector_d> evecs;

    double x = p[0];
    double y = p[1];
    double r = std::sqrt(x*x + y*y);

    x /= r;
    y /= r;

    double lambda = std::exp(-0.5*std::abs(r*r-1));
    double h1 = 0.1 + (1-lambda);
    double h2 = 1.;

    E_Vector v1 = E_Vector::Zero(), v2 = E_Vector::Zero();
    v1(0) = x; v1(1) = y;
    v2(0) = -y; v2(1) = x;
    evecs.push_back(constr_vec(d, v1.data(), v1.data() + d));
    evecs.push_back(constr_vec(d, v2.data(), v2.data() + d));
    evals[0] = 1./(h*h1);
    evals[1] = 1./(h*h2);

    return this->build_metric(evecs, evals);
  }

  Metric interpolation(const Point_d& p) const
  {
    int d = Dim::value;
    typename Kd::Construct_vector_d constr_vec = Kd().construct_vector_d_object();

    std::vector<FT> evals(d);
    std::vector<Vector_d> evecs;
    FT denom = std::sqrt(9.); // = ratio at x=+/-1
    double h = 1.0;

    for(int i=0; i<d; ++i)
    {
      evals[i] = 1./h;

      E_Vector vi = E_Vector::Zero();
      vi(i) = 1.;
      evecs.push_back(constr_vec(d, vi.data(), vi.data() + d));
    }

    FT x = p[0];

    evals[0] = 1./((denom*std::abs(x)+(1-std::abs(x)))*h); //linear interpolation
    //evals[0] = std::pow(1./denom, std::abs(x))/h; //logexp 1D
    //evals[0] = std::max(std::pow(1./denom, std::abs(x)),1./denom); //logexp 1D extended into uniform

    return this->build_metric(evecs, evals);
  }

  virtual Metric compute_metric(const Point_d &p) const
  {
    return yang_liu_cube_shock(p);
    return shock_1D(p);
    return interpolation(p);
  }

  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: custom" << std::endl;
  }

public:
  Custom_metric_field(FT epsilon_ = 1e-6, FT en_factor = 0.999)
  : Metric_field<Kd>(epsilon_, en_factor) { }
};

} // Anisotropic_mesh_TC
} // CGAL

#endif
