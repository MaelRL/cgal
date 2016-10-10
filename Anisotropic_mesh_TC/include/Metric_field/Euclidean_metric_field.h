#ifndef CGAL_ANISOTROPIC_MESH_TC_EUCLIDEAN_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_TC_EUCLIDEAN_METRIC_FIELD

#include <CGAL/Metric_field.h>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Kd>
class Euclidean_metric_field : public Metric_field<Kd>
{
public:
  typedef Metric_field<Kd>                            Base;
  typedef typename Kd::Dimension                      Dim;

  typedef typename Base::FT                           FT;
  typedef typename Base::Metric                       Metric;
  typedef typename Base::Point_d                      Point_d;
  typedef typename Base::Vector_d                     Vector_d;

  typedef typename Base::Metric::E_Vector             E_Vector;

  mutable std::vector<FT> m_evals;

public:
  virtual Metric compute_metric(const Point_d& p) const
  {
    const int d = Dim::value;
    std::vector<Vector_d> evecs;
    typename Kd::Construct_vector_d constr_vec = Kd().construct_vector_d_object();

#if 1 // temporary shock...
    FT x = p[0];
    FT y = p[1];
    FT delta = 0.6;
    FT tanhder = tanh((2.0 * x - sin(5.0 * y)) / delta);
    tanhder = (1.0 - tanhder * tanhder) / delta;
    FT x1 = 2. * tanhder + 3.0 * x * x + y * y;
    FT y1 = -tanhder * cos(5.0 * y) * 5.0 + 2.0 * x * y;

    FT r = sqrt(x1 * x1 + y1 * y1);
    FT x2 = -y1 / r;
    FT y2 = x1 / r;
    FT l1 = std::sqrt(x1*x1 + y1*y1);
    FT l2 = std::sqrt(x2*x2 + y2*y2);

    E_Vector v1 = (1./l1) * E_Vector(x1, y1);
    E_Vector v2 = (1./l2) * E_Vector(x2, y2);

    evecs.push_back(constr_vec(d, v1.data(), v1.data()+d));
    evecs.push_back(constr_vec(d, v2.data(), v2.data()+d));
    m_evals[0] = l1; m_evals[1] = l2;
#else
    for(int i=0; i<d; ++i)
    {
      E_Vector vi = E_Vector::Zero();
      vi(i) = 1.;
      evecs.push_back(constr_vec(d, vi.data(), vi.data() + d));
    }
#endif
//    FT x = p[0];
//    m_evals[0] = 1.1+ std::exp(std::exp(x)); // TMP

    return this->build_metric(evecs, m_evals);
  }

  // this function is used to report the setting of the metric
  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: euclidean" << std::endl;
  }

public:
  Euclidean_metric_field(const std::vector<FT>& evals_ = std::vector<FT>(Dim::value, 1.),
                         FT epsilon_ = 1e-6)
  :
    Metric_field<Kd>(epsilon_), m_evals(evals_) { }
};

} // Anisotropic_mesh_TC
} // CGAL

#endif
