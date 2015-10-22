#ifndef CGAL_ANISOTROPIC_MESH_3_EUCLIDEAN_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_EUCLIDEAN_METRIC_FIELD

#include <CGAL/Metric_field.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Euclidean_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>           Base;
  typedef typename Base::FT         FT;
  typedef typename Base::Metric     Metric;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::Vector_3   Vector_3;

private:
  double m_a, m_b, m_c;
public:
  double& a() { return m_a; }
  double& b() { return m_b; }
  double& c() { return m_c; }

public:
  virtual Metric compute_metric(const Point_3&) const
  {
    return this->build_metric(Vector_3(1, 0, 0),
                              Vector_3(0, 1, 0),
                              Vector_3(0, 0, 1),
                              m_a, m_b, m_c);
  }

  // this function is used to report the setting of the metric
  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: euclidean" << std::endl;
  }

public:
  Euclidean_metric_field(const double a = 1.,
                         const double b = 1.,
                         const double c = 1.,
                         FT epsilon_ = 1e-6,
                         const double en_factor = 1.)
  : Metric_field<K>(epsilon_, en_factor), m_a(a), m_b(b), m_c(c) { }
};

} // Anisotropic_mesh_3
} // CGAL

#endif
