#ifndef CGAL_ANISOTROPIC_MESH_3_IMPLICIT_CURVATURE_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_IMPLICIT_CURVATURE_METRIC_FIELD

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Metric_field.h>
#include <CGAL/Constrain_surface_3_implicit.h>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

    template<typename K>
    class Implicit_curvature_metric_field : public Metric_field<K> {
    public:
      typedef Metric_field<K> Base;
      typedef typename K::FT FT;
      typedef typename Base::Metric Metric;
      typedef typename Base::Point_3 Point_3;
      typedef typename Base::Vector_3 Vector_3;

      typedef Constrain_surface_3_implicit<K> Constrain_surface;
      mutable FT maxl;
    public:
      const Constrain_surface &surface;

      virtual void report(typename std::ofstream &fx) const 
      {
        fx << "type: implicit surface" << std::endl;
        fx << "epsilon: " << this->epsilon << std::endl;
        fx << "maxl: " << maxl << std::endl;
      }

      virtual Metric compute_metric(const Point_3 &p) const 
      {
        Vector_3 en, e1, e2;
        double vp1, vp2;
        surface.tensor_frame(p, en, e1, e2, vp1, vp2);

        //double vpn = surface.global_max_curvature();//global max of c_max
        double vpn = (std::max)(vp1, vp2);

        return this->build_metric(en, e1, e2, vpn, vp1, vp2);
      }

      Implicit_curvature_metric_field(const Constrain_surface &surface_, 
                                      const FT epsilon_ = 1.0,
                                      const double& en_factor_ = 1.) 
        : Metric_field<K>(epsilon_, en_factor_), surface(surface_) { }
    };

  }
}

#endif
