#ifndef CGAL_ANISOTROPIC_MESH_2_METRIC_HELPER_H
#define CGAL_ANISOTROPIC_MESH_2_METRIC_HELPER_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename Kernel>
typename Kernel::Vector_2 get_eigenvector(const Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType::ConstColXpr& v)
{
  return typename Kernel::Vector_2(std::real(v[0]), std::real(v[1]));
}

template<typename Kernel>
Eigen::Matrix2d build_UDUt(const typename Kernel::Vector_2& v_0,
                           const typename Kernel::Vector_2& v_1,
                           const typename Kernel::FT& e_0,
                           const typename Kernel::FT& e_1)
{
  Eigen::Matrix2d eigen_m;
  eigen_m(0,0) = v_0.x();  eigen_m(0,1) = v_1.x();
  eigen_m(1,0) = v_0.y();  eigen_m(1,1) = v_1.y();

  Eigen::Matrix2d eigen_diag = Eigen::Matrix2d::Zero();
  eigen_diag(0,0) = e_0;
  eigen_diag(1,1) = e_1;

  Eigen::Matrix2d eigen_mtransp = eigen_m.transpose();
  Eigen::Matrix2d ret_mat = eigen_m * eigen_diag * eigen_mtransp;

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "build UDUt evs" << std::endl;
  std::cout << e_0 << " " << e_1 << " " << std::endl;
  std::cout << v_0 << std::endl << v_1 << std::endl;
  std::cout << "out UDUt" << std::endl;
  std::cout << ret_mat << std::endl;
#endif
  return ret_mat;
}

template<typename Kernel>
void get_eigen_vecs_and_vals(const Eigen::Matrix2d& matrix,
                             typename Kernel::Vector_2& v0,
                             typename Kernel::Vector_2& v1,
                             double& e_0, double& e_1,
                             const bool abs_value = true)
{
  Eigen::EigenSolver<Eigen::Matrix2d> es(matrix, true);
  const Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType& vecs = es.eigenvectors();
  const Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType& vals = es.eigenvalues();

  v0 = get_eigenvector<Kernel>(vecs.col(0));
  v1 = get_eigenvector<Kernel>(vecs.col(1));

  e_0 = std::real(vals[0]);
  e_1 = std::real(vals[1]);

  if(abs_value)
  {
    e_0 = std::abs(e_0);
    e_1 = std::abs(e_1);
  }

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  if(e_0 < 0 || e_1 < 0)
  {
    std::cout << "WARNING, NEGATIVE EIGENVALUES RETURNED BY EIGEN_VECS" << std::endl;
    std::cout << e_0 << " "<< e_1 << " " << std::endl;
  }
#endif
}

template<typename Kernel>
typename Kernel::FT compute_distortion(const Eigen::Matrix2d& matrix1,
                                       const Eigen::Matrix2d& matrix2)
{
  typename Kernel::Vector_2 v0, v1;
  typename Kernel::FT e0, e1;
  get_eigen_vecs_and_vals<Kernel>(matrix1, v0, v1, e0, e1);
  e0 = std::sqrt(e0); e1 = std::sqrt(e1);
  Eigen::Matrix2d Fp = build_UDUt<Kernel>(v0, v1, e0, e1);
  Eigen::Matrix2d Fpm1 = build_UDUt<Kernel>(v0, v1, 1./e0, 1./e1);

  get_eigen_vecs_and_vals<Kernel>(matrix2, v0, v1, e0, e1);
  e0 = std::sqrt(e0); e1 = std::sqrt(e1);
  Eigen::Matrix2d Fq = build_UDUt<Kernel>(v0, v1, e0, e1);
  Eigen::Matrix2d Fqm1 = build_UDUt<Kernel>(v0, v1, 1./e0, 1./e1);

  double eigen_dis1 = (Fp * Fqm1).operatorNorm();
  double eigen_dis2 = (Fq * Fpm1).operatorNorm();
  return max(eigen_dis1, eigen_dis2);
}

template<typename Kernel>
Eigen::Matrix2d matrix_log(const Eigen::Matrix2d& m, const typename Kernel::FT& alpha = 1)
{
  typename Kernel::Vector_2 v_0, v_1;
  double e_0, e_1;

  get_eigen_vecs_and_vals<Kernel>(m, v_0, v_1, e_0, e_1, false /*no abs*/);

  e_0 = alpha * std::log(e_0);
  e_1 = alpha * std::log(e_1);

  Eigen::Matrix2d ret_mat = build_UDUt<Kernel>(v_0, v_1, e_0, e_1);

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "in log : " << alpha << std::endl;
  std::cout << m << std::endl;
  std::cout << "out log" << std::endl;
  std::cout << ret_mat << std::endl;
#endif
  return ret_mat;
}

template<typename Kernel>
Eigen::Matrix2d matrix_exp(const Eigen::Matrix2d& m)
{
  typename Kernel::Vector_2 v_0, v_1;
  double e_0, e_1;

  get_eigen_vecs_and_vals<Kernel>(m, v_0, v_1, e_0, e_1, false /*no abs*/);

  e_0 = std::exp(e_0);
  e_1 = std::exp(e_1);

  Eigen::Matrix2d ret_mat = build_UDUt<Kernel>(v_0, v_1, e_0, e_1);

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "in exp" << std::endl;
  std::cout << m << std::endl;
  std::cout << "out exp" << std::endl;
  std::cout << ret_mat << std::endl;
#endif
  return ret_mat;
}


template<typename Kernel>
Eigen::Matrix2d scale_matrix_to_point(const Eigen::Matrix2d& matrix_si,
                                      const typename Kernel::Point_2 & si,
                                      const typename Kernel::Point_2 & p,
                                      const typename Kernel::FT beta = 1.1)
{
  Eigen::Vector2d si_p(p.x()-si.x(), p.y()-si.y());
  Eigen::RowVector2d tsi_p = si_p.transpose();
  double sq_dist = tsi_p * matrix_si * si_p;
  double scale_value = 1 + std::sqrt(sq_dist)*std::log10(beta);
  scale_value = 1./(scale_value * scale_value);

  Eigen::Matrix2d scaled_matrix = scale_value * matrix_si;

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
Eigen::Matrix2d matrix_intersection(const Eigen::Matrix2d& M_p,
                                    const Eigen::Matrix2d& M_q,
                                    const bool verbose = false)
{
  Eigen::Matrix2d inverse_M_p;
  //brute force inverse computation : probably more efficient to use the fact the eigenvectors
  //are the same for the inverse and the eigenvalues are the inverses. (todo)
  bool invertible;
  M_p.computeInverseWithCheck(inverse_M_p, invertible);
  if(!invertible)
    std::cout << "M_p not invertible..." << std::endl;
  Eigen::Matrix2d N = inverse_M_p * M_q;

  Eigen::EigenSolver<Eigen::Matrix2d> es(N, true);
  const Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType& vecs = es.eigenvectors();

  Eigen::Vector2d v0 = vecs.col(0).real();
  Eigen::Vector2d v1 = vecs.col(1).real();

  double lambda_0 = v0.transpose() * M_p * v0;
  double lambda_1 = v1.transpose() * M_p * v1;

  double mu_0 = v0.transpose() * M_q * v0;
  double mu_1 = v1.transpose() * M_q * v1;

  double e0 = (std::max)(lambda_0, mu_0);
  double e1 = (std::max)(lambda_1, mu_1);

  Eigen::Matrix2d intersected_eigen_diag = Eigen::Matrix2d::Zero();
  intersected_eigen_diag(0,0) = e0;
  intersected_eigen_diag(1,1) = e1;

  Eigen::Matrix2d real_vecs = vecs.real();
  Eigen::Matrix2d inverse_real_vecs;
  real_vecs.computeInverseWithCheck(inverse_real_vecs, invertible);
  if(!invertible)
    std::cout << "real_vecs not invertible...." << std::endl;

  Eigen::Matrix2d intersected_M = inverse_real_vecs.transpose() * intersected_eigen_diag * inverse_real_vecs;
  if(verbose)
  {
#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
    Eigen::Matrix2d diag_lambda = Eigen::Matrix2d::Zero();
    diag_lambda(0,0) = lambda_0;
    diag_lambda(1,1) = lambda_1;
    Eigen::Matrix2d checkM_p = inverse_real_vecs.transpose() * diag_lambda * inverse_real_vecs;

    Eigen::Matrix2d diag_mu = Eigen::Matrix2d::Zero();
    diag_mu(0,0) = mu_0;
    diag_mu(1,1) = mu_1;
    Eigen::Matrix2d checkM_q = inverse_real_vecs.transpose() * diag_mu * inverse_real_vecs;

    Eigen::Matrix2d tPMpP = real_vecs.transpose() * M_p * real_vecs;
    Eigen::Matrix2d tPMqP = real_vecs.transpose() * M_p * real_vecs;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> ges(M_q,M_p);
    Eigen::Matrix2d gvecs = es.eigenvectors().real();

    Eigen::Matrix2d tPMpP2 = gvecs.transpose() * M_p * gvecs;
    Eigen::Matrix2d tPMqP2 = gvecs.transpose() * M_p * gvecs;

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
    std::cout << lambda_0 << " " << lambda_1 << " " << std::endl;
    std::cout << mu_0 << " " << mu_1 << " " << std::endl;
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
Eigen::Matrix2d logexp_interpolate(const std::vector<std::pair<Eigen::Matrix2d, typename K::FT> >& w_metrics)
{
#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "interpolating (logexp)" << w_metrics.size() << " matrices" << std::endl;
#endif
  Eigen::Matrix2d log_sum = Eigen::Matrix2d::Zero();
  for(std::size_t i=0; i<w_metrics.size(); ++i)
    log_sum += matrix_log<K>(w_metrics[i].first, w_metrics[i].second);

  return matrix_exp<K>(log_sum); // exp(sum(alpha_i*log(Mi)) (== Prod(Mi^alpha_i) if the Mi commute)
}

template<typename K>
Eigen::Matrix2d linear_interpolate(const std::vector<std::pair<Eigen::Matrix2d, typename K::FT> >& w_metrics)
{
#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
  std::cout << "interpolating (linear)" << w_metrics.size() << " matrices" << std::endl;
#endif
  Eigen::Matrix2d res = Eigen::Matrix2d::Zero();
  for(std::size_t i=0; i<w_metrics.size(); ++i)
    for(int j=0; j<2; ++j)
      for(int k=0; k<2; ++k)
        res(j,k) += (w_metrics[i].second) * (w_metrics[i].first)(j,k);

  return res;
}

template<typename Point>
Point transform(const Eigen::Matrix2d& f,
                const Point& p)
{
  Eigen::Vector2d ep(p.x(),p.y());
  ep = f * ep;
  Point result(ep[0],ep[1]);
  return result;
}


} //namespace Aniso
} //namespace CGAL

#endif
