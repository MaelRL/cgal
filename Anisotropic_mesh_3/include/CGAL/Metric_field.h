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

#ifndef CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Metric.h>

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

      virtual Metric compute_metric(const Point_3 &p) const = 0;

      Metric uniform_metric(const Point_3& p) const
      {
        return build_metric(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 
                            1., 1., 1.);
      }

      Metric build_metric(const Vector_3& vn, const Vector_3& v1, const Vector_3& v2,
                          const double& en, const double& e1, const double& e2) const
      {
        return Metric(vn, v1, v2, this->en_factor*en, e1, e2, epsilon);
      }

      //metric_si is the metric at point si. The function scales it to the point p
      Metric scale_metric_to_point(const Metric metric_si,
                                   const Point_3 & si,
                                   const Point_3 & p,
                                   const FT beta = 1.1 /*later to be moved as member of metric_base or somewhere else.*/) const
      {
        std::cout << "scaling from " << si << " to " << p << std::endl;
        Eigen::Vector3d si_p(p.x()-si.x(), p.y()-si.y(), p.z()-si.z());
        Eigen::RowVector3d tsi_p = si_p.transpose();
        double sq_dist = tsi_p * metric_si.get_transformation() * si_p;
        double scale_value = 1 + std::sqrt(sq_dist)*std::log10(beta);
        scale_value = 1./(scale_value * scale_value);

        std::cout << "sq dist " << sq_dist << std::endl;
        std::cout << "scale value " << scale_value << std::endl;
        std::cout << "matrix at si : " << std::endl << metric_si.get_transformation() << std::endl;

        Eigen::Matrix3d scaled_transformation = scale_value * metric_si.get_transformation();

        Eigen::EigenSolver<Eigen::Matrix3d> scaled_es(scaled_transformation, true);
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& scaled_vecs = scaled_es.eigenvectors();
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& scaled_vals = scaled_es.eigenvalues();

        Vector_3 scaled_v0 = get_eigenvector(scaled_vecs.col(0));
        Vector_3 scaled_v1 = get_eigenvector(scaled_vecs.col(1));
        Vector_3 scaled_v2 = get_eigenvector(scaled_vecs.col(2));
        double e_0 = std::abs(std::real(scaled_vals[0]));
        double e_1 = std::abs(std::real(scaled_vals[1]));
        double e_2 = std::abs(std::real(scaled_vals[2]));

        //return Metric(intersected_eigen_transformation);
        Metric scaled_metric = build_metric(scaled_v0, scaled_v1, scaled_v2, e_0, e_1, e_2);

        return scaled_metric;
      }

      Vector_3 get_eigenvector(const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType::ConstColXpr& v) const
      {
        return Vector_3(std::real(v[0]), std::real(v[1]), std::real(v[2]));
      }

      Metric intersection(Metric metric_p, Metric metric_q) const
      {
        Eigen::Matrix3d eigen_transformation_p = metric_p.get_transformation();
        Eigen::Matrix3d eigen_transformation_q = metric_q.get_transformation();
        Eigen::Matrix3d M_p = eigen_transformation_p.transpose()*eigen_transformation_p;
        Eigen::Matrix3d M_q = eigen_transformation_q.transpose()*eigen_transformation_q;

        Eigen::Matrix3d inverse_M_p;
        //brute force inverse computation : probably more efficient to use the fact the eigenvectors
        //are the same for the inverse and the eigenvalues are the inverses. (todo)
        bool invertible;
        M_p.computeInverseWithCheck(inverse_M_p, invertible);
        if(!invertible)
          std::cout << "M_p not invertible..." << std::endl;

        std::cout << "now intersecting..." << std::endl;
        std::cout << "p" << std::endl << eigen_transformation_p << std::endl;
        std::cout << "p evalues : " << metric_p.get_third_eigenvalue() << " " << metric_p.get_max_eigenvalue() << " " << metric_p.get_min_eigenvalue() << std::endl;
        std::cout << "q" << std::endl << eigen_transformation_q << std::endl;
        std::cout << "q evalues : " << metric_q.get_third_eigenvalue() << " " << metric_q.get_max_eigenvalue() << " " << metric_q.get_min_eigenvalue() << std::endl;
        std::cout << "M_p" << std::endl << M_p << std::endl;
        std::cout << "M_p^-1" << std::endl << inverse_M_p << std::endl;
        std::cout << "M_q" << std::endl << M_q << std::endl;

        Eigen::Matrix3d N = inverse_M_p * M_q;

        std::cout << "matrix N : " << std::endl << N << std::endl;

        Eigen::EigenSolver<Eigen::Matrix3d> es(N, true);
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();

        std::cout << "vecs : " << std::endl << vecs << std::endl;
        std::cout << "vals : " << std::endl << vals << std::endl;

        Eigen::Vector3d v0 = vecs.col(0).real();
        Eigen::Vector3d v1 = vecs.col(1).real();
        Eigen::Vector3d v2 = vecs.col(2).real();

        double lambda_0 = v0.transpose() * M_p * v0;
        double lambda_1 = v1.transpose() * M_p * v1;
        double lambda_2 = v2.transpose() * M_p * v2;

        double mu_0 = v0.transpose() * M_q * v0;
        double mu_1 = v1.transpose() * M_q * v1;
        double mu_2 = v2.transpose() * M_q * v2;

        double e0 = (std::max)(lambda_0, mu_0);
        double e1 = (std::max)(lambda_1, mu_1);
        double e2 = (std::max)(lambda_2, mu_2);

        std::cout << "lambdas, mu" << std::endl;
        std::cout << lambda_0 << " " << lambda_1 << " " << lambda_2 << std::endl;
        std::cout << mu_0 << " " << mu_1 << " " << mu_2 << std::endl;

        Eigen::Matrix3d intersected_eigen_diag = Eigen::Matrix3d::Zero();
        intersected_eigen_diag(0,0) = e0;
        intersected_eigen_diag(1,1) = e1;
        intersected_eigen_diag(2,2) = e2;

        std::cout << "intersected eigen diag : " << std::endl << intersected_eigen_diag << std::endl;

        Eigen::Matrix3d real_vecs = vecs.real();
        Eigen::Matrix3d inverse_real_vecs;
        real_vecs.computeInverseWithCheck(inverse_real_vecs, invertible);
        if(!invertible)
          std::cout << "errrrr...." << std::endl;

        std::cout << "real vecs : " << std::endl << real_vecs << std::endl;
        std::cout << "inverse real vecs : " << std::endl << inverse_real_vecs << std::endl;

        Eigen::Matrix3d intersected_eigen_transformation = inverse_real_vecs.transpose() * intersected_eigen_diag * inverse_real_vecs;
        std::cout << "intersected : " << std::endl << intersected_eigen_transformation << std::endl;

        Eigen::EigenSolver<Eigen::Matrix3d> intersected_es(intersected_eigen_transformation, true);
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& intersected_vecs = intersected_es.eigenvectors();
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& intersected_vals = intersected_es.eigenvalues();

        std::cout << "intersected vecs : " << std::endl << intersected_vecs << std::endl;
        std::cout << "intersected ev : " << std::endl << intersected_vals << std::endl;

        Vector_3 intersected_v0 = get_eigenvector(intersected_vecs.col(0));
        Vector_3 intersected_v1 = get_eigenvector(intersected_vecs.col(1));
        Vector_3 intersected_v2 = get_eigenvector(intersected_vecs.col(2));
        double e_0 = std::sqrt(std::abs(std::real(intersected_vals[0]))); // sqrt since those are evalues of M and we want F
        double e_1 = std::sqrt(std::abs(std::real(intersected_vals[1])));
        double e_2 = std::sqrt(std::abs(std::real(intersected_vals[2])));

        //return Metric(intersected_eigen_transformation);
        Metric temp_p = build_metric(intersected_v0, intersected_v1, intersected_v2, e_0, e_1, e_2);
        std::cout << "temp p : " << std::endl << temp_p.get_transformation() << std::endl;

        return temp_p;
      }

      // this function is used to report the setting of the metric
      virtual void report(typename std::ofstream &fx) const = 0;

      Metric_field(FT epsilon_ = 1.0,
                   const double& en_factor_ = 1.)
       : epsilon(epsilon_),
         en_factor(en_factor_)
      { }

      virtual ~Metric_field ( ) { }
    };
  }
}

#endif
