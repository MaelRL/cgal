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

  std::vector<FT> m_evals;

public:
  virtual Metric compute_metric(const Point_d&) const
  {
    const int d = Dim::value;
    std::vector<Vector_d> evecs;
    typename Kd::Construct_vector_d constr_vec = Kd().construct_vector_d_object();
    for(int i=0; i<d; ++i)
    {
      E_Vector vi = E_Vector::Zero();
      vi(i) = 1.;
      evecs.push_back(constr_vec(d, vi.data(), vi.data() + d));
    }

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
