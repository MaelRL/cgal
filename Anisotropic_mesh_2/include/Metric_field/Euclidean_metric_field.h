#ifndef CGAL_ANISOTROPIC_MESH_2_EUCLIDEAN_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_2_EUCLIDEAN_METRIC_FIELD

#include <CGAL/Metric_field.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K>
class Euclidean_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>           Base;
  typedef typename Base::FT         FT;
  typedef typename Base::Metric     Metric;
  typedef typename Base::Point_2    Point_2;
  typedef typename Base::Vector_2   Vector_2;

private:
  double m_a, m_b;
public:
  double& a() { return m_a; }
  double& b() { return m_b; }

public:
  virtual Metric compute_metric(const Point_2 &p) const
  {
    return this->build_metric(Vector_2(1, 0),
                              Vector_2(0, 1),
                              m_a, m_b);
  }

  // this function is used to report the setting of the metric
  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: euclidean" << std::endl;
  }

public:
  Euclidean_metric_field(const double& a = 1.,
                         const double& b = 1.,
                         FT epsilon_ = 1e-6)
  :
    Metric_field<K>(epsilon_), m_a(a), m_b(b) { }
};

} // Anisotropic_mesh_2
} // CGAL

#endif
