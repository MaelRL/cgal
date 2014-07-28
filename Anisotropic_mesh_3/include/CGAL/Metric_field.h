// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Kan-Le Shi

#ifndef CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD_H

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Metric.h>
#include <CGAL/helpers/metric_helper.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename KExact = K>// = CGAL::Exact_predicates_exact_constructions_kernel>
class Metric_field
{
public:
  typedef Metric_base<K, KExact>  Metric;
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_3;
  typedef typename K::Vector_3    Vector_3;

public:
  FT epsilon;
  double en_factor;

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

  virtual void min_max_sq_eigenvalues(const Point_3&, double&, double&) const
  {
    std::cout << "min/max sq evs not defined for this mf" << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const = 0;

  Metric uniform_metric(const Point_3& p) const
  {
    return Metric(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 1., 1., 1., epsilon);
  }

  Metric build_metric(const Vector_3& v0, const Vector_3& v1, const Vector_3& v2,
                      const double& e0, const double& e1, const double& e2) const
  {
    return Metric(v0, v1, v2, en_factor*e0, e1, e2, epsilon);
  }

  //metric_si is the metric at point si. The function scales it to the point p
  Metric scale_metric_to_point(const Metric metric_si,
                               const Point_3 & si,
                               const Point_3 & p,
                               const FT beta = 1.1 /*later to be moved as member of metric_base or somewhere else.*/) const
  {
    Eigen::Matrix3d scaled_transformation = scale_matrix_to_point<K>(metric_si.get_transformation(), si, p, beta);

    Vector_3 scaled_v0, scaled_v1, scaled_v2;
    double e_0, e_1, e_2;
    get_eigen_vecs_and_vals<K>(scaled_transformation, scaled_v0, scaled_v1, scaled_v2, e_0, e_1, e_2);

    //return Metric(intersected_eigen_transformation);
    Metric scaled_metric = build_metric(scaled_v0, scaled_v1, scaled_v2, e_0, e_1, e_2);

    return scaled_metric;
  }

  Metric intersection(Metric metric_p, Metric metric_q) const
  {
    Eigen::Matrix3d eigen_transformation_p = metric_p.get_transformation();
    Eigen::Matrix3d eigen_transformation_q = metric_q.get_transformation();
    Eigen::Matrix3d M_p = eigen_transformation_p.transpose()*eigen_transformation_p;
    Eigen::Matrix3d M_q = eigen_transformation_q.transpose()*eigen_transformation_q;

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
    std::cout << "now intersecting..." << std::endl;
    std::cout << "p" << std::endl << eigen_transformation_p << std::endl;
    std::cout << "p evalues : " << metric_p.get_third_eigenvalue() << " " << metric_p.get_max_eigenvalue() << " " << metric_p.get_min_eigenvalue() << std::endl;
    std::cout << "q" << std::endl << eigen_transformation_q << std::endl;
    std::cout << "q evalues : " << metric_q.get_third_eigenvalue() << " " << metric_q.get_max_eigenvalue() << " " << metric_q.get_min_eigenvalue() << std::endl;
    std::cout << "M_p" << std::endl << M_p << std::endl;
    std::cout << "M_q" << std::endl << M_q << std::endl;
#endif
    Eigen::Matrix3d intersected_metric = matrix_intersection<K>(M_p, M_q);

    Vector_3 intersected_v0, intersected_v1, intersected_v2;
    double e0, e1, e2;
    get_eigen_vecs_and_vals<K>(intersected_metric,
                               intersected_v0, intersected_v1, intersected_v2,
                               e0, e1, e2);

    //need F, not M
    e0 = std::sqrt(e0);
    e1 = std::sqrt(e1);
    e2 = std::sqrt(e2);

    Metric ret_met = build_metric(intersected_v0, intersected_v1, intersected_v2, e0, e1, e2);
    std::cout << "final intersected metric : " << std::endl << ret_met.get_transformation() << std::endl;

    return ret_met;
  }

  // this function is used to report the setting of the metric
  virtual void report(typename std::ofstream &fx) const = 0;

  Metric_field(FT epsilon_ = 1.0,
               FT en_factor_ = 0.999) // to avoid the degen case of two equal evs
    :
      epsilon(epsilon_),
      en_factor(en_factor_),
      m_max_sq_eigenvalue(0.),
      m_min_sq_eigenvalue(1e30),
      m_cache_max_sq_eigenvalue(false),
      m_cache_min_sq_eigenvalue(false)
  { }

  virtual ~Metric_field( ) { }
};

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD_H
