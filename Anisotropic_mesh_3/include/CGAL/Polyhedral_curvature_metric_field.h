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
      typedef Metric_field<K> Base;
      typedef typename K::FT	FT;
      typedef typename Base::Metric	Metric;
      typedef typename Base::Point_3	Point_3;
      typedef typename Base::Vector_3 Vector_3;
      typedef Constrain_surface_3_polyhedral<K>	Constrain_surface;
      FT lambda;
      mutable FT maxl;
    public:
      const Constrain_surface* const m_pConstrain;

      virtual void report(typename std::ofstream &fx) const 
      {
        fx << "type: polyhedral surface" << std::endl;
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
        if (CGAL::abs(vx) > 0.9)
        {
          dx = 0.0;
          dy = 1.0;
          if (CGAL::abs(vy) > 0.9)
          {
            dy = 0.0;
            dz = 1.0;
          }
        }
        Vector_3 d(dx, dy, dz);
        return CGAL::cross_product(v, d);
      }

      virtual Metric compute_metric(const Point_3 &p) const 
      {
        Vector_3 e0, e1, e2;
        double vp0, vp1, vp2;
        ((Constrain_surface *)(m_pConstrain))->tensor_frame(p, e0, e1, e2, vp0, vp1, vp2);
#ifdef ANISO_USE_EIGEN

        //double en = ((Constrain_surface *)(m_pConstrain))->global_max_curvature();
        
        double tmp_epsilon = lambda; // to get it from command-line 'main'

        return Metric(e0, e1, e2, vp0, vp1, vp2, tmp_epsilon);
#else
        FT e0len = sqrt(e0.x() * e0.x() + e0.y() * e0.y() + e0.z() * e0.z());
        FT e1len = sqrt(e1.x() * e1.x() + e1.y() * e1.y() + e1.z() * e1.z());
        FT e2len = sqrt(e2.x() * e2.x() + e2.y() * e2.y() + e2.z() * e2.z());
        e0len = e0len * lambda;
        e1len = e1len * lambda;
        e2len = e2len * lambda;
        FT ub = 500.0;
        if (e0len > ub) e0len = ub;
        if (e1len > ub) e1len = ub;
        if (e2len > ub) e2len = ub;
        if (e0len < 1e-8) {
          if (e1len < 1e-8) {
            if (e2len < 1e-8) {
              e0 = Vector_3(1, 0, 0);
              e1 = Vector_3(0, 1, 0);
              e2 = Vector_3(0, 0, 1);
            } else {
              e0 = get_orthogonal(e2);
              e1 = get_unit(CGAL::cross_product(e2, e0));
              e2 = get_unit(e2) * (e2len + 1.0);
            }
          } else {
            if (e2len < 1e-8) {
              e0 = get_orthogonal(e1);
              e2 = get_unit(CGAL::cross_product(e0, e1));
              e1 = get_unit(e1) * (e1len + 1.0);
            } else {
              e0 = get_unit(CGAL::cross_product(e1, e2));
              e1 = get_unit(e1) * (e1len + 1.0);
              e2 = get_unit(e2) * (e2len + 1.0);
            }
          }
        } else {
          if (e1len < 1e-8) {
            if (e2len < 1e-8) {
              e1 = get_orthogonal(e0);
              e2 = get_unit(CGAL::cross_product(e0, e1));
              e0 = get_unit(e0) * (e0len + 1.0);
            } else {
              e1 = get_unit(CGAL::cross_product(e2, e0));
              e0 = get_unit(e0) * (e0len + 1.0);
              e2 = get_unit(e2) * (e2len + 1.0);
            }
          } else {
            if (e2len < 1e-8) {
              e2 = get_unit(CGAL::cross_product(e0, e1));
              e0 = get_unit(e0) * (e0len + 1.0);
              e1 = get_unit(e1) * (e1len + 1.0);
            } else {
              e0 = get_unit(e0) * (e0len + 1.0);
              e1 = get_unit(e1) * (e1len + 1.0);
              e2 = get_unit(e2) * (e2len + 1.0);
            }
          }
        }
        double epsilon = 0.;
        return Metric(e0, e1, e2, 1., 1., 1., epsilon);
#endif
      }

      Polyhedral_curvature_metric_field(const Constrain_surface& surface_, 
                                        const FT lambda_ = 1.0) 
        : Metric_field<K>(),  
          m_pConstrain(&surface_), 
          lambda(lambda_) { }
    };

  }
}

#endif
