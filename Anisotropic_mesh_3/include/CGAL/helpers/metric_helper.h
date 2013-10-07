#ifndef CGAL_ANISOTROPIC_MESH_3_METRIC_HELPER_H
#define CGAL_ANISOTROPIC_MESH_3_METRIC_HELPER_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


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
Eigen::Matrix3d build_UDUt (const typename Kernel::Vector_3& v_0,
                            const typename Kernel::Vector_3& v_1,
                            const typename Kernel::Vector_3& v_2,
                            const typename Kernel::FT& e_0,
                            const typename Kernel::FT& e_1,
                            const typename Kernel::FT& e_2)
{
  Eigen::Matrix3d eigen_m;
  eigen_m(0,0) = v_0.x();  eigen_m(0,1) = v_1.x();  eigen_m(0,2) = v_2.x();
  eigen_m(1,0) = v_0.y();  eigen_m(1,1) = v_1.y();  eigen_m(1,2) = v_2.y();
  eigen_m(2,0) = v_0.z();  eigen_m(2,1) = v_1.z();  eigen_m(2,2) = v_2.z();

  Eigen::Matrix3d eigen_diag = Eigen::Matrix3d::Zero();
  eigen_diag(0,0) = e_0;
  eigen_diag(1,1) = e_1;
  eigen_diag(2,2) = e_2;

  Eigen::Matrix3d eigen_mtransp = eigen_m.transpose();
  Eigen::Matrix3d ret_mat = eigen_m * eigen_diag * eigen_mtransp;

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "build UDUt evs" << std::endl;
  std::cout << e_0 << " " << e_1 << " " << e_2 << std::endl;
  std::cout << v_0 << std::endl << v_1 << std::endl << v_2 << std::endl;
  std::cout << "out UDUt" << std::endl;
  std::cout << ret_mat << std::endl;
#endif
  return ret_mat;
}

template<typename Kernel>
void get_eigen_vecs_and_vals(const Eigen::Matrix3d& matrix,
                             typename Kernel::Vector_3& v0,
                             typename Kernel::Vector_3& v1,
                             typename Kernel::Vector_3& v2,
                             double& e_0, double& e_1, double& e_2,
                             const bool abs_value = true)
{
  Eigen::EigenSolver<Eigen::Matrix3d> es(matrix, true);
  const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();
  const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();

  v0 = get_eigenvector<Kernel>(vecs.col(0));
  v1 = get_eigenvector<Kernel>(vecs.col(1));
  v2 = get_eigenvector<Kernel>(vecs.col(2));

  e_0 = std::real(vals[0]);
  e_1 = std::real(vals[1]);
  e_2 = std::real(vals[2]);

  if(abs_value)
  {
    e_0 = std::abs(e_0);
    e_1 = std::abs(e_1);
    e_2 = std::abs(e_2);
  }

#ifdef ANISO_DEBUG_NEG_EIGENVALUES
  if(e_0 < 0 || e_1 < 0 || e_2 < 0)
  {
    std::cout << "WARNING, NEGATIVE EIGENVALUES RETURNED BY EIGEN_VECS" << std::endl;
    std::cout << e_0 << " "<< e_1 << " "<< e_2 << std::endl;
  }
#endif
}

template<typename Kernel>
Eigen::Matrix3d matrix_log(const Eigen::Matrix3d& m, const typename Kernel::FT& alpha = 1)
{
  typename Kernel::Vector_3 v_0, v_1, v_2;
  double e_0, e_1, e_2;

  get_eigen_vecs_and_vals<Kernel>(m, v_0, v_1, v_2, e_0, e_1, e_2, false /*no abs*/);

  e_0 = alpha * std::log(e_0);
  e_1 = alpha * std::log(e_1);
  e_2 = alpha * std::log(e_2);

  Eigen::Matrix3d ret_mat = build_UDUt<Kernel>(v_0, v_1, v_2, e_0, e_1, e_2);

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "in log : " << alpha << std::endl;
  std::cout << m << std::endl;
  std::cout << "out log" << std::endl;
  std::cout << ret_mat << std::endl;
#endif
  return ret_mat;
}

template<typename Kernel>
Eigen::Matrix3d matrix_exp(const Eigen::Matrix3d& m)
{
  typename Kernel::Vector_3 v_0, v_1, v_2;
  double e_0, e_1, e_2;

  get_eigen_vecs_and_vals<Kernel>(m, v_0, v_1, v_2, e_0, e_1, e_2, false /*no abs*/);

  e_0 = std::exp(e_0);
  e_1 = std::exp(e_1);
  e_2 = std::exp(e_2);

  Eigen::Matrix3d ret_mat = build_UDUt<Kernel>(v_0, v_1, v_2, e_0, e_1, e_2);

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "in exp" << std::endl;
  std::cout << m << std::endl;
  std::cout << "out exp" << std::endl;
  std::cout << ret_mat << std::endl;
#endif
  return ret_mat;
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
  std::cout << "scaled matrix : " << std::endl << scaled_matrix << std::endl;
#endif

  return scaled_matrix;
}

template<typename Kernel>
Eigen::Matrix3d matrix_intersection(const Eigen::Matrix3d& M_p,
                                    const Eigen::Matrix3d& M_q,
                                    const bool verbose = false)
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

  Eigen::Matrix3d intersected_M = inverse_real_vecs.transpose() * intersected_eigen_diag * inverse_real_vecs;

  Eigen::Matrix3d diag_lambda = Eigen::Matrix3d::Zero();
  diag_lambda(0,0) = lambda_0;
  diag_lambda(1,1) = lambda_1;
  diag_lambda(2,2) = lambda_2;
  Eigen::Matrix3d checkM_p = inverse_real_vecs.transpose() * diag_lambda * inverse_real_vecs;

  Eigen::Matrix3d diag_mu= Eigen::Matrix3d::Zero();
  diag_mu(0,0) = mu_0;
  diag_mu(1,1) = mu_1;
  diag_mu(2,2) = mu_2;
  Eigen::Matrix3d checkM_q = inverse_real_vecs.transpose() * diag_mu * inverse_real_vecs;

  Eigen::Matrix3d tPMpP = real_vecs.transpose() * M_p * real_vecs;
  Eigen::Matrix3d tPMqP = real_vecs.transpose() * M_p * real_vecs;

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> ges(M_q,M_p);
  Eigen::Matrix3d gvecs = es.eigenvectors().real();

  Eigen::Matrix3d tPMpP2 = gvecs.transpose() * M_p * gvecs;
  Eigen::Matrix3d tPMqP2 = gvecs.transpose() * M_p * gvecs;

  if(verbose)
  {
#if 1//def ANISO_DEBUG_MATRIX_OPERATIONS
    std::cout.precision(15);
    std::cout << "in : " << std::endl << M_p << std::endl << inverse_M_p << std::endl << M_q << std::endl;
    std::cout << "check tPMP" << std::endl << tPMpP << std::endl << tPMqP << std::endl;
    std::cout << "check tPMP with generalized eigenvalues" << std::endl << tPMpP2 << std::endl << tPMqP2 << std::endl;
    std::cout << "check lambdamuP : " << std::endl << checkM_p << std::endl << checkM_q << std::endl;
    std::cout << "check inverse : " << std::endl << inverse_M_p*M_p << std::endl;
    std::cout << "matrix N : " << std::endl << N << std::endl;
    std::cout << "The eigenvalues of the pencil (A,B) are:" << std::endl << ges.eigenvalues() << std::endl;
    std::cout << "The matrix of eigenvectors, V, is:" << std::endl << ges.eigenvectors() << std::endl;
    std::cout << "lambdas, mu" << std::endl;
    std::cout << lambda_0 << " " << lambda_1 << " " << lambda_2 << std::endl;
    std::cout << mu_0 << " " << mu_1 << " " << mu_2 << std::endl;
    std::cout << "intersected eigen diag : " << std::endl << intersected_eigen_diag << std::endl;
    std::cout << "real vecs : " << std::endl << real_vecs << std::endl;
    std::cout << "inverse real vecs : " << std::endl << inverse_real_vecs << std::endl;
    std::cout << "check inverse : "<< std::endl << real_vecs*inverse_real_vecs << std::endl;
    std::cout << "intersected : " << std::endl << intersected_M << std::endl;
#endif
  }

  return intersected_M;
}

template<typename K>
Eigen::Matrix3d interpolate_colors(const std::vector<std::pair<Eigen::Matrix3d, typename K::FT> >& w_metrics)
{
#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "interpolating " << w_metrics.size() << " matrices" << std::endl;
#endif
  Eigen::Matrix3d log_sum = Eigen::Matrix3d::Zero();
  for(std::size_t i=0; i<w_metrics.size(); ++i)
    log_sum += matrix_log<K>(w_metrics[i].first, w_metrics[i].second);

  return matrix_exp<K>(log_sum); //Prod(Miˆalphai)
}

template<typename K>
void get_metric_normal_index(const typename K::Vector_3& normal, const std::vector<typename K::Vector_3>& v, int& index)
{
  typedef typename K::FT FT;
  FT sp0 = std::abs(v[0] * normal);
  FT sp1 = std::abs(v[1] * normal);
  FT sp2 = std::abs(v[2] * normal);
  FT max_sp = (std::max)(sp0, (std::max)(sp1, sp2));

  index = 0;
  if(sp1 == max_sp)
    index = 1;
  else if(sp2 == max_sp)
    index = 2;

  //todo maybe, swap the normal if normal and the metric vector closest to normal have negative scalar prod
}

} //namespace Aniso
} //namespace CGAL

#endif
