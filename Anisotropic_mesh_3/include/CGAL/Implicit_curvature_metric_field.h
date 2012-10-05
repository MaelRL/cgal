#ifndef CGAL_ANISOTROPIC_MESH_3_IMPLICIT_CURVATURE_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_IMPLICIT_CURVATURE_METRIC_FIELD

#include <iostream>
#include <fstream>
#include <utility>

#include "Metric_field.h"
#include "Constrain_surface_3_implicit.h"

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
      FT lambda;
      mutable FT maxl;
    public:
      const Constrain_surface &surface;

      virtual void report(typename std::ofstream &fx) const 
      {
        fx << "type: implicit surface" << std::endl;
        fx << "lambda: " << lambda << std::endl;
        fx << "maxl: " << maxl << std::endl;
      }

      Vector_3 adjust_length(const Vector_3 &v) const 
      {
        FT l = sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
        FT newl = l * lambda + 1.0;
        l = newl / l;
        if (newl > maxl)
          maxl = newl;
        return Vector_3(v.x() * l, v.y() * l, v.z() * l);
      }

      static Vector_3 get_unit(const Vector_3 &v)
      {
        FT vx = v.x(), vy = v.y(), vz = v.z();
        FT vl = sqrt(vx * vx + vy * vy + vz * vz);
        return Vector_3(vx / vl, vy / vl, vz / vl);
      }

      static Vector_3 get_orthogonal(const Vector_3 &v)
      {
        FT vx = v.x(), vy = v.y(), vz = v.z();
        FT vl = sqrt(vx * vx + vy * vy + vz * vz);
        vx /= vl; vy /= vl; vz /= vl;
        FT dx = 1.0, dy = 0, dz = 0;
        if (fabs(vx) > 0.9)
        {
          dx = 0.0;
          dy = 1.0;
          if (fabs(vy) > 0.9)
          {
            dy = 0.0;
            dz = 1.0;
          }
        }
        Vector_3 d(dx, dy, dz);
        return CGAL::cross_product(v, d);
      }

      bool is_zero(const double &d, const double& epsilon = 1e-10) const
      {
        return (std::abs(d) < 1e-10);
      }

      virtual Metric compute_metric(const Point_3 &p) const 
      {
        Vector_3 en, e1, e2;
        double vp1, vp2;
        surface.tensor_frame(p, en, e1, e2, vp1, vp2, 1e-5);
#ifdef ANISO_USE_EIGEN
//        return Metric(n, e1, e2, 5*lambda, (vp1*lambda + 1.), (vp2*lambda + 1.));

        double vpn = surface.global_max_curvature();//global max of c_max
        double tmp_epsilon = lambda; // to get it from command-line 'main'

        //double epsilon = 0.01 * global_max_curvature

        return Metric(en, e1, e2, vpn, vp1, vp2, tmp_epsilon);
#else
        Vector_3 e0 = n;
        FT e1len = e1.x() * e1.x() + e1.y() * e1.y() + e1.z() * e1.z();
        FT e2len = e2.x() * e2.x() + e2.y() * e2.y() + e2.z() * e2.z();
        if (e1len < 1e-8) {
          if (e2len < 1e-8) {
            // two zero eigen values
            e1 = get_orthogonal(e0);
            e2 = CGAL::cross_product(e0, e1);
          } else {
            e1 = get_unit(CGAL::cross_product(e2, e0));
            e2 = adjust_length(e2);
          }
        } else {
          if (e2len < 1e-6) {
            e1 = adjust_length(e1);
            e2 = get_unit(CGAL::cross_product(e0, e1));
          } else {
            e1 = adjust_length(e1);
            e2 = adjust_length(e2);
          }
        }
        FT rate = 5.0 * lambda;
        e0 = Vector_3(e0.x() * rate, e0.y() * rate, e0.z() * rate);

        double epsilon = 1e-6;
        return Metric(e0, e1, e2, rate, vp1, vp2, epsilon);
#endif
      }

      Implicit_curvature_metric_field(const Constrain_surface &surface_, 
                                      const FT lambda_ = 1.0) 
        : Metric_field<K>(), surface(surface_), lambda(lambda_) { }
    };

  }
}

#endif
