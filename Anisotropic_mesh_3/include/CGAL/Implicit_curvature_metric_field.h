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
class Implicit_curvature_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>                                       Base;
  typedef typename K::FT                                        FT;
  typedef typename Base::Metric                                 Metric;
  typedef typename Base::Point_3                                Point_3;
  typedef typename Base::Vector_3                               Vector_3;

  typedef Constrain_surface_3_implicit<K>                       Constrain_surface;
  typedef typename Constrain_surface::C3t3                      C3t3;

public:
  const Constrain_surface& m_pConstrain;

  virtual double global_max_sq_eigenvalue() const
  {
    if(this->m_cache_max_sq_eigenvalue)
      return this->m_max_sq_eigenvalue;

    this->m_max_sq_eigenvalue = 0.;
    this->m_min_sq_eigenvalue = DBL_MAX;
    typename C3t3::Triangulation::Finite_vertices_iterator v;
    for(v = m_pConstrain.c3t3().triangulation().finite_vertices_begin();
        v != m_pConstrain.c3t3().triangulation().finite_vertices_end();
        ++v)
    {
      if(m_pConstrain.c3t3().in_dimension(v) == 1 ||
         m_pConstrain.c3t3().in_dimension(v) == 2)
      {
        double minc, maxc;
        min_max_sq_eigenvalues(v->point(), minc, maxc);
        this->m_max_sq_eigenvalue = (std::max)(this->m_max_sq_eigenvalue, maxc);
        this->m_min_sq_eigenvalue = (std::min)(this->m_min_sq_eigenvalue, minc);
      }
    }
    this->m_cache_max_sq_eigenvalue = true;
    this->m_cache_min_sq_eigenvalue = true;
    return this->m_max_sq_eigenvalue;
  }

  virtual double global_min_sq_eigenvalue() const
  {
    if(this->m_cache_min_sq_eigenvalue)
      return this->m_min_sq_eigenvalue;

    global_max_sq_eigenvalue(); //computes both max and min
    return this->m_min_sq_eigenvalue;
  }

  void min_max_sq_eigenvalues(const Point_3& p,
                          double& minc,
                          double& maxc) const
  {
    Vector_3 n, v1, v2;
    double en, e1, e2;
    tensor_frame(p, n, v1, v2, en, e1, e2, 1e-5);
    minc = (std::min)(e1, e2);
    maxc = (std::max)(e1, e2);
  }

  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: implicit surface" << std::endl;
    fx << "epsilon: " << this->epsilon << std::endl;
  }

  virtual Vector_3 gradient(const Point_3 &p, const FT delta = 1e-5) const
  {
    FT dd = 1. / (2. * delta);
    return dd * Vector_3(
          (m_pConstrain.evaluate(p.x() + delta, p.y(), p.z()) -
           m_pConstrain.evaluate(p.x() - delta, p.y(), p.z())),
          (m_pConstrain.evaluate(p.x(), p.y() + delta, p.z()) -
           m_pConstrain.evaluate(p.x(), p.y() - delta, p.z())),
          (m_pConstrain.evaluate(p.x(), p.y(), p.z() + delta) -
           m_pConstrain.evaluate(p.x(), p.y(), p.z() - delta)));
  }

  //virtual Vector_3 normal(const Point_3 &p, const FT delta = 1e-5) const
  //{
  //  Vector_3 n = gradient(p, delta);
  //  FT inv_len = -1. / sqrt(n*n);
  //  return inv_len * n;virtual
  //}

  void orthogonal_vectors(const Vector_3& u,
                          Vector_3& v,
                          Vector_3& w) const
  {
    typename K::Line_3 support(CGAL::ORIGIN, CGAL::ORIGIN + u);
    typename K::Plane_3 p = support.perpendicular_plane(CGAL::ORIGIN);
    v = p.base1();
    w = p.base2();
  }

  int are_zeros(const double& d1, const double& d2, const double& d3,
                bool& z1, bool& z2, bool& z3) const
  {
    int count = 0;
    z1 = (std::abs(d1) < ZERO_EIGENVALUE);  if(z1) count++;
    z2 = (std::abs(d2) < ZERO_EIGENVALUE);  if(z2) count++;
    z3 = (std::abs(d3) < ZERO_EIGENVALUE);  if(z3) count++;
    return count;
  }

  bool are_equal(const Vector_3& v1,        //unit vector
                 const Vector_3& v2) const  //unit vector
  {
    return (std::abs(std::abs(v1*v2) - 1.) < 1e-10);
  }

  Vector_3 get_vector(const Eigen::Vector3d& v) const
  {
    return Vector_3(v[0], v[1], v[2]);
  }

  Vector_3 normalize(const Vector_3& v) const
  {
    return std::sqrt(1./(v*v)) * v;
  }

  Eigen::Matrix3d hessian(const Point_3 &p, const FT delta = 1e-5) const
  {
    Vector_3 xp = gradient(Point_3(p.x() + delta, p.y(), p.z()), delta);
    Vector_3 xn = gradient(Point_3(p.x() - delta, p.y(), p.z()), delta);
    Vector_3 yp = gradient(Point_3(p.x(), p.y() + delta, p.z()), delta);
    Vector_3 yn = gradient(Point_3(p.x(), p.y() - delta, p.z()), delta);
    Vector_3 zp = gradient(Point_3(p.x(), p.y(), p.z() + delta), delta);
    Vector_3 zn = gradient(Point_3(p.x(), p.y(), p.z() - delta), delta);
    FT dd = 1. / (2. * delta);
    Eigen::Matrix3d m;
    m(0,0) = dd*(xp.x() - xn.x()); m(0,1) = dd*(yp.x() - yn.x()); m(0,2) = dd*(zp.x() - zn.x());
    m(1,0) = dd*(xp.y() - xn.y()); m(1,1) = dd*(yp.y() - yn.y()); m(1,2) = dd*(zp.y() - zn.y());
    m(2,0) = dd*(xp.z() - xn.z()); m(2,1) = dd*(yp.z() - yn.z()); m(2,2) = dd*(zp.z() - zn.z());
    return m;
  }

  void tensor_frame(const Point_3 &p,
                    Vector_3 &v0, //unit normal
                    Vector_3 &v1, //unit eigenvector
                    Vector_3 &v2, //unit eigenvector
                    double& e0, //eigenvalue corresponding to v0
                    double& e1, //eigenvalue corresponding to v1
                    double& e2, //eigenvalue corresponding to v2
                    const FT delta = 1e-5) const
  {
    // method in the book
    Vector_3 gx = gradient(p, delta);
    Eigen::Vector3d grad(gx.x(), gx.y(), gx.z());
    double gl = grad.norm();
    if(gl != 0.)
      grad.normalize();
    else
      std::cout << "Warning : gradient is null vector at "<< p <<"!" << std::endl;
    Eigen::Vector3d normal = -grad;

    Eigen::Matrix3d PN = Eigen::Matrix3d::Identity() - (normal * normal.transpose());
    Eigen::Matrix3d hess = (1./gl) * hessian(p, delta) ;
    Eigen::Matrix3d m = (PN * hess * PN);

#ifdef ANISO_DEBUG_METRIC
    std::cout.precision(10);
    std::cout << "Normal : " << std::endl;
    std::cout << normal << std::endl;
    std::cout << "NN^t : " << std::endl;
    std::cout << normal * normal.transpose() << std::endl;
    std::cout << "Matrices before EVs computing : " << std::endl;
    std::cout << "PN : " << std::endl;
    std::cout << PN << std::endl;
    std::cout << "Hessian with grad_len = " << gl << std::endl;
    std::cout << hess << std::endl;
    std::cout << "M : " << std::endl;
    std::cout << m << std::endl;
#endif

    Eigen::EigenSolver<Eigen::Matrix3d> es(m, true/*compute eigenvectors and values*/);
    const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();
    const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();

    // look for smallest eigenvalue
    FT min_abs_value = DBL_MAX;
    int min_id = 0;
    for (int i = 0; i < 3; i++)
    {
      double avi = std::abs(std::real(vals[i]));
      if(avi < min_abs_value)
      {
        min_abs_value = avi;
        min_id = i;
      }
    }

    // normal is ok (computed with gradient)
    v0 = get_vector(normal);

    bool z0, z1, z2;
    double ev0 = std::abs(std::real(vals[min_id])); //should be the smallest
    double ev1 = std::abs(std::real(vals[(min_id + 1) % 3]));
    double ev2 = std::abs(std::real(vals[(min_id + 2) % 3]));
    int zeros = are_zeros(ev0, ev1, ev2, z0, z1, z2);

    if(zeros == 1) //it is ev0
    {
      if(!z0) std::cout << "Error1 : see tensor_frame, zeros==1\n";
      Vector_3 normal_test = get_eigenvector<K>(vecs.col(min_id));
      if(!are_equal(v0, normal_test))
        std::cout << "Error2 : see tensor_frame, zeros==1\n";

      e1 = ev1;
      e2 = ev2;
      v1 = get_eigenvector<K>(vecs.col((min_id + 1) % 3));
      v2 = get_eigenvector<K>(vecs.col((min_id + 2) % 3));

      // make sure the vectors form an orthogonal matrix
      // when eigenvalues are the same, vectors returned by Eigen are not
      // necessarily orthogonal
      if(std::abs(v1*v2) > ZERO_DOT_PRODUCT)
        v2 = normalize(CGAL::cross_product(v0, v1));
    }
    else if(zeros == 2)
    {
      int id = -1; // the id of the non-zero eigenvalue
      if(z0 && z1)      id = ((min_id + 2) % 3); //d2
      else if(z0 && z2) id = ((min_id + 1) % 3); //d1
      else if(z1 && z2) id = min_id;//d0
      else std::cerr << "Error1 : see tensor_frame, zeros==2\n";

      // the non-zero one
      e1 = std::abs(std::real(vals[id]));
      v1 = get_eigenvector<K>(vecs.col(id));
      // the last one : choose the vector which is not // to the normal
      Vector_3 normal_test_1 = get_eigenvector<K>(vecs.col((id + 1) % 3));
      Vector_3 normal_test_2 = get_eigenvector<K>(vecs.col((id + 2) % 3));
      if(are_equal(normal_test_1, v0))
      {
        e2 = std::abs(std::real(vals[(min_id + 2) % 3]));
        v2 = normal_test_2;
      }
      else if(are_equal(normal_test_2, v0))
      {
        e2 = std::abs(std::real(vals[(min_id + 1) % 3]));
        v2 = normal_test_1;
      }
      else
      {
        e2 = 0.;
        v2 = normalize(CGAL::cross_product(v0, v1));
        std::cerr << "Error2 : see tensor_frame, zeros==2\n";
      }
    }
    else if(zeros == 3)
    {
      orthogonal_vectors(v0, v1, v2);
      v1 = normalize(v1);
      v2 = normalize(v2);
      e1 = 0.;
      e2 = 0.;
    }
    else
    { //zeros == 0
      std::cout << "Error : see tensor_frame, zeros==0\n";
      std::cout << "p : " << p << std::endl;
      std::cout << ev0 << " " << ev1 << " " << ev2 << std::endl;
      std::cout << m << std::endl;
    }

    e0 = (std::max)(e1, e2);
  }

  virtual Metric compute_metric(const Point_3& p) const
  {
    Vector_3 vn, v1, v2;
    double en, e1, e2;
    tensor_frame(p, vn, v1, v2, en, e1, e2);
    return this->build_metric(vn, v1, v2, std::sqrt(en), std::sqrt(e1), std::sqrt(e2));

    Point_3 p_axis(0., 0., 0.);
    Point_3 p_surf;
    m_pConstrain.intersection(p_axis,
                         Point_3(p_axis + 10*Vector_3(p_axis, p)),
                         p_axis).assign(p_surf);
    tensor_frame(p_surf, vn, v1, v2, en, e1, e2);
    Metric m_surf = this->build_metric(vn, v1, v2, std::sqrt(en), std::sqrt(e1), std::sqrt(e2));

#if 1
    FT val = m_pConstrain.evaluate(p.x(), p.y(), p.z());
    FT lvlset = -0.5;
    if(val < lvlset)
      return Metric();

    FT lambda = 1. - val / lvlset;
#else
    FT d_pa_ps = CGAL::sqrt(CGAL::squared_distance(p_axis, p_surf));
    FT d_pa_p = CGAL::sqrt(CGAL::squared_distance(p_axis, p));
    FT lambda = std::min(1., d_pa_p / d_pa_ps);
#endif

    std::vector<std::pair<Eigen::Matrix3d, FT> > w_metrics;
    w_metrics.push_back(std::make_pair(Eigen::Matrix3d::Identity(), 1.-lambda));
    w_metrics.push_back(std::make_pair(m_surf.get_transformation(), lambda));
    Eigen::Matrix3d m_interpolated = logexp_interpolate<K>(w_metrics);

    get_eigen_vecs_and_vals<K>(m_interpolated, vn, v1, v2, en, e1, e2);
    Metric res = this->build_metric(vn, v1, v2, en, e1, e2);

    /*
        std::cout << "paxis : " << p_axis << std::endl;
        std::cout << "p : " << p << std::endl;
        std::cout << "psurf : " << p_surf << std::endl;
        std::cout << "dpapslambda: " << d_pa_ps << " " << d_pa_p << " " << lambda << std::endl;
        std::cout << "out: " << en << " " << e1 << " " << e2 << std::endl;
        std::cout << vn << std::endl << v1 << std::endl << v2 << std::endl;
*/

    return res;
  }

  Implicit_curvature_metric_field(const Constrain_surface& surface_,
                                  FT epsilon_ = 1e-6, FT en_factor_ = 0.999)
    :
      Metric_field<K>(epsilon_, en_factor_),
      m_pConstrain(surface_)
  { }
};

}
}

#endif
