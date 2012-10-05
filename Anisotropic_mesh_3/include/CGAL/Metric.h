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
//
// This class provides a definition of a metric by three orthogonal eigen-vectors.
// The magnitude of the vectors are the squared eigen value.
// For example, (3, 0, 0), (2, 0, 0), (1, 0, 0) denote that the space is elongated
// in (1, 0, 0) three times, in (2, 0, 0) two-times.
// Denote the three vectors V1, V2 and V3, and the eigen vectors v1, v2 and v3,
// squared eigen-values of the tensor metric M is e1, e2, e3, we have:
//		V1 = v1 * e1, V2 = v2 * e2, V3 = v3 * e3
//		R = [v1, v2, v3]
//		F_M = R diag{e1, e2, e3} R^{-1}
//      M = F_M * F_M^t

#ifndef CGAL_ANISOTROPIC_MESH_3_METRIC
#define CGAL_ANISOTROPIC_MESH_3_METRIC

#include <boost/algorithm/minmax_element.hpp>
#include <CGAL/bbox.h>
#include <iostream>
#include <fstream>
#include <utility>


#ifdef ANISO_USE_EIGEN
#include <Eigen/Dense>
#else 
#include <Klein/matrix3x3.h>
#endif

namespace CGAL{
namespace Anisotropic_mesh_3{

template<typename K, typename KExact>
class Metric_base {
public:
    typedef typename K::FT                      FT;
    typedef typename K::Vector_3                Vector_3;
#ifndef ANISO_USE_EIGEN
    typedef typename K::Aff_transformation_3    Aff_transformation_3;
    typedef typename KExact::Aff_transformation_3 Aff_transformation_exact_3;
#endif

private:
#ifndef ANISO_USE_EIGEN
    Aff_transformation_3 transformation;
    Aff_transformation_exact_3 transformation_exact;
    Aff_transformation_3 inverse_transformation;
    Aff_transformation_exact_3 inverse_transformation_exact;
    typename Klein::Kernel::Maths::Matrix3x3 transformation_matrix;
    typename Klein::Kernel::Maths::Matrix3x3 inverse_transformation_matrix;
#else 
  Eigen::Matrix3d eigen_transformation, eigen_inverse_transformation;
  double e_max, e_min, e_n;
#endif
  Vector_3 v_max, v_min, v_n;


public:
#ifndef ANISO_USE_EIGEN
  const Aff_transformation_3& get_transformation() const { return transformation; }
#else
  const Eigen::Matrix3d& get_transformation() const { return eigen_transformation; }
#endif

  const Vector_3& get_vmin() const { return v_min; }
  const Vector_3& get_vmax() const { return v_max; }
  const Vector_3& get_vn() const { return v_n; }
  
public:
  // Transform
    template <typename Object>
    Object transform(const Object &p) const {
      return transform(p, typename Kernel_traits<Object>::Kernel());
    }

    template<typename Object>
    Object transform(const Object &p, K k_) const {
#ifdef ANISO_USE_EIGEN
      Eigen::Vector3d ep(p.x(),p.y(),p.z());
      ep = eigen_transformation * ep;
      Object result(ep[0],ep[1],ep[2]);
      return result;
#else
      Object r =  p.transform(transformation);
      return r;
#endif
    }

#ifdef ANISO_USE_EXACT
    template<typename Object>
    Object transform(const Object &p, KExact k_exact) const {
      return p.transform(transformation_exact);
    }
#endif

    // Inverse transform
    template <typename Object>
    Object inverse_transform(const Object &p) const {
      return inverse_transform(p, typename Kernel_traits<Object>::Kernel());
    }


#ifdef ANISO_USE_EIGEN
  template <typename Kernel>
  Triangle_3<Kernel> inverse_transform(const Triangle_3<Kernel>& t) const {
    
    typename Kernel::Point_3 p = t.vertex(0);
    typename Kernel::Point_3 q = t.vertex(1);
    typename Kernel::Point_3 r = t.vertex(2);
    Eigen::Vector3d ep(p.x(),p.y(),p.z());
    ep = eigen_inverse_transformation * ep;
    p=Point_3<Kernel>(ep[0],ep[1],ep[2]);

    ep = Eigen::Vector3d(q.x(),q.y(),q.z());
    ep = eigen_inverse_transformation * ep;
    q=Point_3<Kernel>(ep[0],ep[1],ep[2]);


    ep = Eigen::Vector3d(r.x(),r.y(),r.z());
    ep = eigen_inverse_transformation * ep;
    r=Point_3<Kernel>(ep[0],ep[1],ep[2]);
    
    return Triangle_3<Kernel>(p,q,r);
  }


  template <typename Kernel>
  Bbox<Kernel> inverse_transform(const Bbox<Kernel>& bb) const {
    
    double x[8],y[8],z[8];
    Eigen::Vector3d ev(bb.xmin(), bb.ymin(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[0] = ev[0]; y[0] = ev[1]; z[0] = ev[2];


    ev = Eigen::Vector3d(bb.xmin(), bb.ymin(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[1] = ev[0]; y[1] = ev[1]; z[1] = ev[2];


    ev = Eigen::Vector3d(bb.xmin(), bb.ymax(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[2] = ev[0]; y[2] = ev[1]; z[2] = ev[2];


    ev = Eigen::Vector3d(bb.xmin(), bb.ymax(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[3] = ev[0]; y[3] = ev[1]; z[3] = ev[2];


    ev = Eigen::Vector3d(bb.xmax(), bb.ymin(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[4] = ev[0]; y[4] = ev[1]; z[4] = ev[2];


    ev = Eigen::Vector3d(bb.xmax(), bb.ymin(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[5] = ev[0]; y[5] = ev[1]; z[5] = ev[2];


    ev = Eigen::Vector3d(bb.xmax(), bb.ymax(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[6] = ev[0]; y[6] = ev[1]; z[6] = ev[2];


    ev = Eigen::Vector3d(bb.xmax(), bb.ymax(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[7] = ev[0]; y[7] = ev[1]; z[7] = ev[2];

    std::pair<double*,double*> xmm = boost::minmax_element(x,x+8);
    std::pair<double*,double*> ymm = boost::minmax_element(y,y+8);
    std::pair<double*,double*> zmm = boost::minmax_element(z,z+8);
    return Bbox<Kernel>(*(xmm.first), *(ymm.first), *(zmm.first),
                   *(xmm.second), *(ymm.second), *(zmm.second));
  }
#endif

    template<typename Object, typename K_obj>
    Object inverse_transform(const Object &p, K_obj k_obj) const {
      
#ifdef ANISO_USE_EIGEN
      Eigen::Vector3d ep(p.x(),p.y(),p.z());
      ep = eigen_inverse_transformation * ep;
      
      Object result(ep[0],ep[1],ep[2]);
      return result;
#else
      Object r = p.transform(inverse_transformation);
      return r;
#endif 
    }

    template<typename Object>
    Object inverse_transform(const Object &p, K k_) const {
#ifdef ANISO_USE_EIGEN      
      Eigen::Vector3d ep(p.x(),p.y(),p.z());
      ep = eigen_inverse_transformation * ep;
      
      Object result(ep[0],ep[1],ep[2]);
      return result;
#else
      Object r = p.transform(inverse_transformation);
      return r;
#endif
    }

#ifdef ANISO_USE_EXACT
    template<typename Object>
    Object inverse_transform(const Object &p, KExact k_exact) const {
        return p.transform(inverse_transformation_exact);
    }
#endif

    Metric_base &operator=(const Metric_base &m) {
#ifndef ANISO_USE_EIGEN
        transformation = m.transformation;
        inverse_transformation = m.inverse_transformation;
        transformation_matrix = m.transformation_matrix;
        inverse_transformation_matrix = m.inverse_transformation_matrix;
#ifdef ANISO_USE_EXACT
        transformation_exact = m.transformation_exact;
        inverse_transformation_exact = m.inverse_transformation_exact;
#endif
#else
        eigen_transformation = m.eigen_transformation;
        eigen_inverse_transformation = m.eigen_inverse_transformation;
#endif
        return *this;
    }
    

  FT compute_distortion(const Metric_base &m) 
  {
#ifdef ANISO_USE_EIGEN
    double eigen_dis1 = (m.eigen_transformation * eigen_inverse_transformation).operatorNorm();
    double eigen_dis2 = (eigen_transformation * m.eigen_inverse_transformation).operatorNorm();
    return max(eigen_dis1, eigen_dis2);
#else       
    Klein::Kernel::Maths::Matrix3x3 u, d;
    Klein::Kernel::Maths::eigenValue(inverse_transformation_matrix * m.transformation_matrix, u, d);
    FT dis1 = max(max(fabs(d.data_[0]), fabs(d.data_[4])), fabs(d.data_[8]));

    Klein::Kernel::Maths::eigenValue(transformation_matrix * m.inverse_transformation_matrix, u, d);
    FT dis2 = max(max(fabs(d.data_[0]), fabs(d.data_[4])), fabs(d.data_[8]));

    return max(dis1, dis2);
#endif
    }

    void construct(const Vector_3 &axis_x, //normal
                   const Vector_3 &axis_y, 
                   const Vector_3 &axis_z,
                   const double& vpn,
                   const double& vpmax, //curvature max
                   const double& vpmin, //curvature min
                   const FT& epsilon) 
    {
#ifdef ANISO_USE_EIGEN      
      e_n = std::max(epsilon, std::abs(vpn));
      e_max = std::max(epsilon, std::abs(vpmax)); //vpmax
      e_min = std::max(epsilon, std::abs(vpmin)); //vpmin
            
      v_n = axis_x;
      v_max = axis_y;
      v_min = axis_z;

      Eigen::Matrix3d eigen_m;
      eigen_m(0,0) = v_n.x();  eigen_m(0,1) = v_max.x();  eigen_m(0,2) = v_min.x();
      eigen_m(1,0) = v_n.y();  eigen_m(1,1) = v_max.y();  eigen_m(1,2) = v_min.y();
      eigen_m(2,0) = v_n.z();  eigen_m(2,1) = v_max.z();  eigen_m(2,2) = v_min.z();     

      Eigen::Matrix3d eigen_diag = Eigen::Matrix3d::Zero();
      eigen_diag(0,0) = e_n;  //should be max(e_max)
      eigen_diag(1,1) = e_max;
      eigen_diag(2,2) = e_min;
              
      Eigen::Matrix3d eigen_mtransp = eigen_m.transpose();
      eigen_transformation = eigen_m * eigen_diag * eigen_mtransp;

      eigen_diag(0,0) = 1.0/eigen_diag(0,0); 
      eigen_diag(1,1) = 1.0/eigen_diag(1,1);
      eigen_diag(2,2) = 1.0/eigen_diag(2,2);       
      eigen_inverse_transformation = eigen_m * eigen_diag * eigen_mtransp;
      
#else      
        v_max = (1./std::sqrt(axis_y*axis_y))*axis_y;
        v_min = (1./std::sqrt(axis_z*axis_z))*axis_z;

        Klein::Kernel::Maths::Matrix3x3 m(
            axis_x.x(), axis_x.y(), axis_x.z(),
            axis_y.x(), axis_y.y(), axis_y.z(),
            axis_z.x(), axis_z.y(), axis_z.z());
        Klein::Kernel::Maths::Matrix3x3 minv;
        m.inverse(minv);

        Klein::Kernel::Maths::Matrix3x3 diag(
            sqrt(axis_x.squared_length()), 0, 0,
            0, sqrt(axis_y.squared_length()), 0,
            0, 0, sqrt(axis_z.squared_length()));
        transformation_matrix = m * diag * minv;

        diag.data_[0] = 1.0 / diag.data_[0];
        diag.data_[4] = 1.0 / diag.data_[4];
        diag.data_[8] = 1.0 / diag.data_[8];
        inverse_transformation_matrix = m * diag * minv;

        build_transforms();

      //Klein::Kernel::Maths::Matrix3x3 test1 = transformation_matrix * inverse_transformation_matrix;
      //std::cout << "\nIdentity?(No Eigen)\n";
      //print_klein_matrix(test1);
      //Klein::Kernel::Maths::Matrix3x3 test2 = inverse_transformation_matrix * transformation_matrix;
      //std::cout << "  reverse :\n";
      //print_klein_matrix(test2);
#endif
    }

#ifndef ANISO_USE_EIGEN
    void print_klein_matrix(const Klein::Kernel::Maths::Matrix3x3& m) const
    {
      std::cout << m.data_[0] << "\t" << m.data_[3] << "\t" << m.data_[6] << std::endl;
      std::cout << m.data_[1] << "\t" << m.data_[4] << "\t" << m.data_[7] << std::endl;
      std::cout << m.data_[2] << "\t" << m.data_[5] << "\t" << m.data_[8] << std::endl;
    }
#endif


    void build_transforms() {
#ifndef ANISO_USE_EIGEN        
        transformation = Aff_transformation_3(
            *(transformation_matrix[0]), *(transformation_matrix[1]), *(transformation_matrix[2]),
            *(transformation_matrix[0] + 1), *(transformation_matrix[1] + 1), *(transformation_matrix[2] + 1),
            *(transformation_matrix[0] + 2), *(transformation_matrix[1] + 2), *(transformation_matrix[2] + 2));
        inverse_transformation = Aff_transformation_3(
            *(inverse_transformation_matrix[0]), *(inverse_transformation_matrix[1]), *(inverse_transformation_matrix[2]),
            *(inverse_transformation_matrix[0] + 1), *(inverse_transformation_matrix[1] + 1), *(inverse_transformation_matrix[2] + 1),
            *(inverse_transformation_matrix[0] + 2), *(inverse_transformation_matrix[1] + 2), *(inverse_transformation_matrix[2] + 2));

#ifdef ANISO_USE_EXACT
        transformation_exact = Aff_transformation_exact_3(
            *(transformation_matrix[0]), *(transformation_matrix[1]), *(transformation_matrix[2]),
            *(transformation_matrix[0] + 1), *(transformation_matrix[1] + 1), *(transformation_matrix[2] + 1),
            *(transformation_matrix[0] + 2), *(transformation_matrix[1] + 2), *(transformation_matrix[2] + 2));
        inverse_transformation_exact = Aff_transformation_exact_3(
            *(inverse_transformation_matrix[0]), *(inverse_transformation_matrix[1]), *(inverse_transformation_matrix[2]),
            *(inverse_transformation_matrix[0] + 1), *(inverse_transformation_matrix[1] + 1), *(inverse_transformation_matrix[2] + 1),
            *(inverse_transformation_matrix[0] + 2), *(inverse_transformation_matrix[1] + 2), *(inverse_transformation_matrix[2] + 2));
#endif
#endif
    }

    void get_min_eigenvector(Vector_3& v) const
    {
      v = v_min;
    }

    void get_max_eigenvector(Vector_3& v) const
    {
      v = v_max;
    }
    void get_third_eigenvector(Vector_3& v) const
    {
      v = v_n;
    }
#ifdef ANISO_USE_EIGEN
    double get_min_eigenvalue() const
    {
      return e_min;
    }
    double get_max_eigenvalue() const
    {
      return e_max;
    }
#endif

public:
#ifdef ANISO_USE_EIGEN
    Metric_base(const Metric_base &m) :
      eigen_transformation(m.eigen_transformation),
      eigen_inverse_transformation(m.eigen_inverse_transformation),
      e_max(m.e_max),
      e_min(m.e_min),
      e_n(m.e_n),
      v_max(m.v_max),
      v_min(m.v_min),
      v_n(m.v_n)
      {}
#else
    Metric_base(const Metric_base &m) :
        transformation_matrix(m.transformation_matrix),
        inverse_transformation_matrix(m.inverse_transformation_matrix)
    {
      v_max = m.v_max;
      v_min = m.v_min;
      v_n = m.v_n;
      build_transforms(); 
    }
#endif

public:
    std::ostream& operator<<(const Metric_base& x)
    {
#ifdef ANISO_USE_EIGEN
      out << "M  = " << m_eigen_transformation << std::endl;
#else
      out << "M  = " << transformation_matrix << std::endl;
#endif
      return out;
    }

    Metric_base() {
        construct(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 1, 1, 1, 1e-6);
    }
    Metric_base(const Vector_3 &axis_x,  //normal (unit)
                const Vector_3 &axis_y,  //principal direction 1 (unit) // axis_vpmax
                const Vector_3 &axis_z,  //principal direction 2 (unit) // axis_vpmin
                const double& vpn,
                const double& vpy,
                const double& vpz,
                const double& epsilon)
    {      
      //std::cout << "produits scalaires : " << axis_x*axis_z << " " 
      //                                     << axis_x*axis_y << " " 
      //                                     << axis_y*axis_z << std::endl;
      //std::cout << "norms              : " << std::sqrt(axis_x*axis_x) << " " 
      //                                     << std::sqrt(axis_y*axis_y) << " " 
      //                                     << std::sqrt(axis_z*axis_z) << std::endl;
      //std::cout << "valeurs propres    : " << vpn << " " << vpy << " " << vpz << std::endl;
      //std::cout << std::endl;
      
      if(std::abs(vpy) > std::abs(vpz))
        construct(axis_x, axis_y, axis_z, vpn, vpy, vpz, epsilon);
      else
        construct(axis_x, axis_z, axis_y, vpn, vpz, vpy, epsilon);  
    }

    ~Metric_base() { }
};
}
}

#endif
