#ifndef CGAL_ANISOTROPIC_MESH_TC_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_TC_METRIC_FIELD_H

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Metric.h>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Kd>
class Metric_field
{
public:
  typedef Metric_base<Kd>                      Metric;
  typedef typename Kd::FT                      FT;
  typedef typename Kd::Point_d                 Point_d;
  typedef typename Kd::Vector_d                Vector_d;

public:
  FT epsilon;
  double en_factor;

  virtual Metric compute_metric(const Point_d& p) const = 0;

  Metric uniform_metric(const Point_d& /*p*/) const
  {
    return Metric();
  }

  Metric build_metric(const std::vector<Vector_d>& evecs,
                      const std::vector<FT>& evals) const
  {
    // en_factor was used to avoid degeneracies, but it's more complicated than
    // automatically slightly changing one of the eigenvalues. What if (for example)
    // the degen is in the first d-1 eigenvalues being the same ?
    // todo

//    int d = Kd::Dimension::value;
//    if(d>1)
//      evals[d-1] *= en_factor;

    return Metric(evecs, evals, epsilon);
  }

  // this function is used to report the setting of the metric
  virtual void report(typename std::ofstream &fx) const = 0;

  Metric_field(FT epsilon_ = 1e-6,
               FT en_factor_ = 0.999) // to avoid the degen case of two equal evs
    :
      epsilon(epsilon_),
      en_factor(en_factor_)
  { }

  virtual ~Metric_field( ) { }
};

} // Anisotropic_mesh_TC
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_TC_METRIC_FIELD_H
