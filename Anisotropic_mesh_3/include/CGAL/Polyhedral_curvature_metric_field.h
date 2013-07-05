#ifndef CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD

#include <iostream>
#include <fstream>
#include <utility>


#include <CGAL/Metric_field.h>
#include <CGAL/Constrain_surface_3_polyhedral.h>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

    template<typename K>
    class Polyhedral_curvature_metric_field : public Metric_field<K> {
    public:
      typedef Metric_field<K>                   Base;
      typedef typename K::FT                    FT;
      typedef typename Base::Metric             Metric;
      typedef typename Base::Point_3            Point_3;
      typedef typename Base::Vector_3           Vector_3;
      typedef Constrain_surface_3_polyhedral<K> Constrain_surface;

    public:
      mutable FT maxl;
      const Constrain_surface* const m_pConstrain;

      virtual void report(typename std::ofstream &fx) const 
      {
        fx << "type: polyhedral surface" << std::endl;
        fx << "epsilon: " << this->epsilon << std::endl;
        fx << "maxl: " << maxl << std::endl;
      }

      virtual Metric compute_metric(const Point_3 &p) const 
      {
        Vector_3 v0, v1, v2;
        double e0, e1, e2;
        ((Constrain_surface *)(m_pConstrain))->tensor_frame(p, v0, v1, v2, e0, e1, e2);
      
        return this->build_metric(v0, v1, v2, std::sqrt(e0), std::sqrt(e1), std::sqrt(e2));
      }

      Polyhedral_curvature_metric_field(const Constrain_surface& surface_, 
                                        const FT epsilon_ = 1.0,
                                        const double& en_factor_ = 1.0)
        : Metric_field<K>(epsilon_, en_factor_), 
          m_pConstrain(&surface_)
        { }
    };

  }
}

#endif
