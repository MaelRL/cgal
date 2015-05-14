#ifndef CGAL_ANISOTROPIC_MESH_2_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_2_METRIC_FIELD_H

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Metric.h>
#include <CGAL/helpers/metric_helper.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K, typename KExact = K>// = CGAL::Exact_predicates_exact_constructions_kernel>
class Metric_field
{
public:
  typedef Metric_base<K, KExact>  Metric;
  typedef typename K::FT          FT;
  typedef typename K::Point_2     Point_2;
  typedef typename K::Vector_2    Vector_2;

public:
  FT epsilon;

  mutable double m_max_sq_eigenvalue;
  mutable double m_min_sq_eigenvalue;
  mutable bool m_cache_max_sq_eigenvalue;
  mutable bool m_cache_min_sq_eigenvalue;

  virtual double global_max_sq_eigenvalue() const
  {
    std::cout << "max sq ev not defined for this mf" << std::endl;
    return 0.;
  }

  virtual double global_min_sq_eigenvalue() const
  {
    std::cout << "min sq ev not defined for this mf" << std::endl;
    return 0.;
  }

  virtual void min_max_sq_eigenvalues(const Point_2&, double&, double&) const
  {
    std::cout << "min/max sq evs not defined for this mf" << std::endl;
  }

  virtual Metric compute_metric(const Point_2 &p) const = 0;

  Metric uniform_metric(const Point_2& /*p*/) const
  {
    return Metric(Vector_2(1, 0, 0), Vector_2(0, 1, 0), 1., 1., epsilon);
  }

  Metric build_metric(const Vector_2& v0, const Vector_2& v1,
                      const double& e0, const double& e1) const
  {
    return Metric(v0, v1, e0, e1, epsilon);
  }

  // this function is used to report the setting of the metric
  virtual void report(typename std::ofstream &fx) const = 0;

  Metric_field(FT epsilon_ = 1e-6)
    :
      epsilon(epsilon_),
      m_max_sq_eigenvalue(0.),
      m_min_sq_eigenvalue(1e30),
      m_cache_max_sq_eigenvalue(false),
      m_cache_min_sq_eigenvalue(false)
  { }

  virtual ~Metric_field( ) { }
};

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_METRIC_FIELD_H
