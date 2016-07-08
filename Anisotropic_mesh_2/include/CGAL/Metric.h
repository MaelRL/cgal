#ifndef CGAL_ANISOTROPIC_MESH_2_METRIC_H
#define CGAL_ANISOTROPIC_MESH_2_METRIC_H

#include <CGAL/bbox.h>
#include <CGAL/Triangle_2.h>

#include <boost/algorithm/minmax_element.hpp>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K, typename KExact = K>
class Metric_base
{
public:
  typedef typename K::FT                      FT;
  typedef typename K::Point_2                 Point_2;
  typedef typename K::Vector_2                Vector_2;

public:
  Eigen::Matrix2d eigen_transformation, eigen_inverse_transformation;
  Eigen::Matrix2d full_matrix;
  mutable FT e_max, e_min;
  Vector_2 v_max, v_min;

public:
  const Eigen::Matrix2d& get_transformation() const { return eigen_transformation; }
  const Eigen::Matrix2d& get_inverse_transformation() const { return eigen_inverse_transformation; }
  const Vector_2& get_vmin() const { return v_min; }
  const Vector_2& get_vmax() const { return v_max; }
  void get_min_eigenvector(Vector_2& v) const { v = v_min; }
  void get_max_eigenvector(Vector_2& v) const { v = v_max; }
  FT get_anisotropic_ratio() const { return e_max/e_min; }
  double get_min_eigenvalue() const { return e_min; }
  double get_max_eigenvalue() const { return e_max; }
  void set_min_eigenvalue(double emin) const { e_min = emin; }
  void set_max_eigenvalue(double emax) const { e_max = emax; }

  const Eigen::Matrix2d& get_mat() const
  {
    return full_matrix;
  }

public:
  // Transform
  template <typename Object>
  Object transform(const Object &p) const
  {
    Eigen::Vector2d ep(p.x(),p.y());
    ep = eigen_transformation * ep;
    Object result(ep[0],ep[1]);
    return result;
  }

  /* transformation_exact doesn't exist anymore
#ifdef ANISO_USE_EXACT
  template<typename Object>
  Object transform(const Object &p, KExact k_exact) const {
    return p.transform(transformation_exact);
  }
#endif
  */

  // Inverse transform
  template <typename Object>
  Object inverse_transform(const Object &p) const
  {
    Eigen::Vector2d ep(p.x(),p.y());
    ep = eigen_inverse_transformation * ep;

    Object result(ep[0],ep[1]);
    return result;
  }

  template <typename Kernel>
  CGAL::Triangle_2<Kernel> inverse_transform(const CGAL::Triangle_2<Kernel>& t) const
  {
    typename Kernel::Point_2 p = t.vertex(0);
    typename Kernel::Point_2 q = t.vertex(1);

    Eigen::Vector2d ep(p.x(), p.y());
    ep = eigen_inverse_transformation * ep;
    p = Point_2(ep[0], ep[1]);

    ep = Eigen::Vector2d(q.x(), q.y());
    ep = eigen_inverse_transformation * ep;
    q = Point_2(ep[0], ep[1]);

    return CGAL::Triangle_2<Kernel>(p, q);
  }

  template <typename Kernel>
  Bbox<Kernel> inverse_transform(const Bbox<Kernel>& bb) const
  {
    double x[4],y[4];
    Eigen::Vector2d ev(bb.xmin(), bb.ymin());
    ev = eigen_inverse_transformation * ev;
    x[0] = ev[0]; y[0] = ev[1];

    ev = Eigen::Vector2d(bb.xmin(), bb.ymax());
    ev = eigen_inverse_transformation * ev;
    x[1] = ev[0]; y[1] = ev[1];

    ev = Eigen::Vector2d(bb.xmax(), bb.ymax());
    ev = eigen_inverse_transformation * ev;
    x[2] = ev[0]; y[2] = ev[1];

    ev = Eigen::Vector2d(bb.xmax(), bb.ymin());
    ev = eigen_inverse_transformation * ev;
    x[3] = ev[0]; y[3] = ev[1];

    std::pair<double*,double*> xmm = boost::minmax_element(x,x+4);
    std::pair<double*,double*> ymm = boost::minmax_element(y,y+4);
    return Bbox<Kernel>(*(xmm.first), *(ymm.first),
                        *(xmm.second), *(ymm.second));
  }

  FT compute_distortion(const Metric_base &m) const
  {
    double eigen_dis1 = (m.eigen_transformation * eigen_inverse_transformation).operatorNorm();
    double eigen_dis2 = (eigen_transformation * m.eigen_inverse_transformation).operatorNorm();

#ifdef ANISO_DEBUG_METRIC
    double max_dist = 1.001;
    if(max(eigen_dis1, eigen_dis2) > max_dist)
    {
      std::cout << "--------" << std::endl;
      std::cout << "this evs : " << get_third_eigenvalue() << " " << get_min_eigenvalue() << std::endl;
      std::cout << "m evs : " << m.get_third_eigenvalue() << " " << m.get_min_eigenvalue() << std::endl;
      std::cout << "--THIS-ET---" << std::endl;
      std::cout << eigen_transformation << std::endl;
      std::cout << "--THIS-EIT---" << std::endl;
      std::cout << eigen_inverse_transformation << std::endl;
      std::cout << "--M-ET---" << std::endl;
      std::cout << m.eigen_transformation << std::endl;
      std::cout << "--M-EIT---" << std::endl;
      std::cout << m.eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "--M-ET-*-THIS-EIT-" << std::endl;
      std::cout << m.eigen_transformation * eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "--THIS-ET-*-M-EIT" << std::endl;
      std::cout << eigen_transformation * m.eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "vals dist1,dist2 : " << eigen_dis1 << " " << eigen_dis2 << std::endl;
    }
#endif

    return max(eigen_dis1, eigen_dis2);
  }

  Eigen::Matrix2d get_inverse_mat() const
  {
    FT inv_sq_e_max = 1./(e_max*e_max);
    FT inv_sq_e_min = 1./(e_min*e_min);

    Eigen::Matrix2d eigen_m;
    eigen_m(0,0) = v_max.x();  eigen_m(0,1) = v_min.x();
    eigen_m(1,0) = v_max.y();  eigen_m(1,1) = v_min.y();

    Eigen::Matrix2d eigen_diag = Eigen::Matrix2d::Zero();
    eigen_diag(0,0) = inv_sq_e_max;
    eigen_diag(1,1) = inv_sq_e_min;

    return eigen_m * eigen_diag * eigen_m.transpose();
  }

  void construct(const Vector_2 &axis_x, //normal
                 const Vector_2 &axis_y,
                 const double& vpmax, //curvature max
                 const double& vpmin, //curvature min
                 const FT& epsilon)
  {
    e_max = (std::max)(epsilon, std::abs(vpmax)); //vpmax
    e_min = (std::max)(epsilon, std::abs(vpmin)); //vpmin

    v_max = axis_x;
    v_min = axis_y;

    Eigen::Matrix2d eigen_m;
    eigen_m(0,0) = v_max.x();  eigen_m(0,1) = v_min.x();
    eigen_m(1,0) = v_max.y();  eigen_m(1,1) = v_min.y();

    Eigen::Matrix2d eigen_diag = Eigen::Matrix2d::Zero();
    eigen_diag(0,0) = e_max;
    eigen_diag(1,1) = e_min;

    Eigen::Matrix2d eigen_mtransp = eigen_m.transpose();
    eigen_transformation = eigen_m * eigen_diag * eigen_mtransp;

    full_matrix = eigen_transformation.transpose()*eigen_transformation;

#ifdef ANISO_DEBUG_METRIC
    std::cout << "--- eigen_m ---" << std::endl;
    std::cout << eigen_m << std::endl;
    std::cout << "--- eigen_diag ---" << std::endl;
    std::cout << eigen_diag << std::endl;
    std::cout << "-- eigen_mtransp ---" << std::endl;
    std::cout << eigen_mtransp << std::endl;
    std::cout << "-- eigen_transformation ---" << std::endl;
    std::cout << eigen_transformation << std::endl;
#endif

    eigen_diag(0,0) = 1.0/eigen_diag(0,0);
    eigen_diag(1,1) = 1.0/eigen_diag(1,1);
    eigen_inverse_transformation = eigen_m * eigen_diag * eigen_mtransp;

#ifdef ANISO_DEBUG_METRIC
    std::cout << "-- eigen_inverse_transformation ---" << std::endl;
    std::cout << eigen_inverse_transformation << std::endl;
#endif
  }

public:
  friend
  std::ostream& operator<<(std::ostream& out, const Metric_base& x)
  {
    out << "M  = " << x.eigen_transformation << std::endl;
    return out;
  }

  Metric_base()
  {
    construct(Vector_2(1, 0), Vector_2(0, 1), 1, 1, 1e-6);
  }

  Metric_base(Eigen::Matrix2d eigen_transformation_)
  {
    //std::cout << "building a metric from a Matrix2d with input : " << std::endl;
    //std::cout << eigen_transformation_ << std::endl;

    construct(Vector_2(1, 0), Vector_2(0, 1), 1, 1, 1e-6);
    eigen_transformation = eigen_transformation_;

    bool invertible;
    double determinant;
    eigen_transformation.computeInverseAndDetWithCheck(eigen_inverse_transformation, determinant, invertible);
/*
    if(invertible)
      std::cout << "It is invertible, and its inverse is:" << std::endl << eigen_inverse_transformation << std::endl;
    else
      std::cout << "Not invertible and determinant was : " << determinant << std::endl;
*/
    full_matrix = eigen_transformation.transpose()*eigen_transformation;
  }

  Metric_base(const Vector_2 &axis_x,
              const Vector_2 &axis_y,
              const double& vpx,
              const double& vpy,
              const double& epsilon)
  {
#ifdef ANISO_DEBUG_METRIC
    double xy = axis_x*axis_y;
    double xx = axis_x*axis_x;
    double yy = axis_y*axis_y;

    std::cout << "produits scalaires : " << xy << std::endl;
    std::cout << "sq_norms : " << xx << " " << yy << std::endl;
    std::cout << "valeurs propres    : " << vpx << " " << vpy << std::endl;
    std::cout << std::endl;
    if(xy > 0.01)
      std::cout << "Warning : This dot product should be 0" << std::endl;
    if(std::abs(xx-1.) > 0.01 || std::abs(yy-1.) > 0.01)
      std::cout << "Warning : This sqnorm should be 1" << std::endl;
#endif
    if(std::abs(vpx) > std::abs(vpy)) // is this really useful ?
      construct(axis_x, axis_y, vpx, vpy, epsilon);
    else
      construct(axis_y, axis_x, vpy, vpx, epsilon);
  }

  Metric_base(const Metric_base& m_)
    :
      eigen_transformation(m_.eigen_transformation),
      eigen_inverse_transformation(m_.eigen_inverse_transformation),
      full_matrix(m_.full_matrix),
      e_max(m_.e_max),
      e_min(m_.e_min),
      v_max(m_.v_max),
      v_min(m_.v_min)
  { }

  ~Metric_base() { }
};

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_METRIC_H
