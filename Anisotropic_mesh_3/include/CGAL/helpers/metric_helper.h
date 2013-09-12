#ifndef METRIC_HELPER_H
#define METRIC_HELPER_H

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Kernel>
typename Kernel::Vector_3 get_eigenvector(const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType::ConstColXpr& v)
{
  return typename Kernel::Vector_3(std::real(v[0]), std::real(v[1]), std::real(v[2]));
}

template<typename Kernel>
void get_eigen_vecs_and_vals(const Eigen::Matrix3d& matrix,
                             typename Kernel::Vector_3& v0,
                             typename Kernel::Vector_3& v1,
                             typename Kernel::Vector_3& v2,
                             double& e_0, double& e_1, double& e_2)
{
  Eigen::EigenSolver<Eigen::Matrix3d> es(matrix, true);
  const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();
  const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();

  v0 = get_eigenvector<Kernel>(vecs.col(0));
  v1 = get_eigenvector<Kernel>(vecs.col(1));
  v2 = get_eigenvector<Kernel>(vecs.col(2));

  e_0 = std::abs(std::real(vals[0]));
  e_1 = std::abs(std::real(vals[1]));
  e_2 = std::abs(std::real(vals[2]));
}

template<typename Kernel>
Eigen::Matrix3d scale_matrix_to_point(const Eigen::Matrix3d& matrix_si,
                                      const typename Kernel::Point_3 & si,
                                      const typename Kernel::Point_3 & p,
                                      const typename Kernel::FT beta = 1.1)
{
  Eigen::Vector3d si_p(p.x()-si.x(), p.y()-si.y(), p.z()-si.z());
  Eigen::RowVector3d tsi_p = si_p.transpose();
  double sq_dist = tsi_p * matrix_si * si_p;
  double scale_value = 1 + std::sqrt(sq_dist)*std::log10(beta);
  scale_value = 1./(scale_value * scale_value);

  Eigen::Matrix3d scaled_matrix = scale_value * matrix_si;

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "scaling from " << si << " to " << p << std::endl;
  std::cout << "sq dist " << sq_dist << std::endl;
  std::cout << "scale value " << scale_value << std::endl;
  std::cout << "matrix at si : " << std::endl << matrix_si << std::endl;
#endif

  return scaled_matrix;
}

template<typename Kernel>
Eigen::Matrix3d matrix_intersection(const Eigen::Matrix3d& M_p, const Eigen::Matrix3d& M_q)
{
  Eigen::Matrix3d inverse_M_p;
  //brute force inverse computation : probably more efficient to use the fact the eigenvectors
  //are the same for the inverse and the eigenvalues are the inverses. (todo)
  bool invertible;
  M_p.computeInverseWithCheck(inverse_M_p, invertible);
  if(!invertible)
    std::cout << "M_p not invertible..." << std::endl;

  Eigen::Matrix3d N = inverse_M_p * M_q;

  Eigen::EigenSolver<Eigen::Matrix3d> es(N, true);
  const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();

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

  Eigen::Matrix3d intersected_eigen_diag = Eigen::Matrix3d::Zero();
  intersected_eigen_diag(0,0) = e0;
  intersected_eigen_diag(1,1) = e1;
  intersected_eigen_diag(2,2) = e2;

  Eigen::Matrix3d real_vecs = vecs.real();
  Eigen::Matrix3d inverse_real_vecs;
  real_vecs.computeInverseWithCheck(inverse_real_vecs, invertible);
  if(!invertible)
    std::cout << "real_vecs not invertible...." << std::endl;

  Eigen::Matrix3d intersected_eigen_transformation = inverse_real_vecs.transpose() * intersected_eigen_diag * inverse_real_vecs;

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "matrix N : " << std::endl << N << std::endl;
  std::cout << "vecs : " << std::endl << vecs << std::endl;
  std::cout << "vals : " << std::endl << vals << std::endl;
  std::cout << "lambdas, mu" << std::endl;
  std::cout << lambda_0 << " " << lambda_1 << " " << lambda_2 << std::endl;
  std::cout << mu_0 << " " << mu_1 << " " << mu_2 << std::endl;
  std::cout << "intersected eigen diag : " << std::endl << intersected_eigen_diag << std::endl;
  std::cout << "real vecs : " << std::endl << real_vecs << std::endl;
  std::cout << "inverse real vecs : " << std::endl << inverse_real_vecs << std::endl;
  std::cout << "intersected : " << std::endl << intersected_eigen_transformation << std::endl;
#endif

  return intersected_eigen_transformation;
}

} //namespace Aniso
} //namespace CGAL

#endif
