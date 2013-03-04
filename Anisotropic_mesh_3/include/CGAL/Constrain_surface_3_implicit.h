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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_IMPLICIT_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_IMPLICIT_H

#define ZERO_DOT_PRODUCT 1e-5
#define ZERO_EIGENVALUE 1e-10

#include <CGAL/Constrain_surface_3_ex.h>

#include <Eigen/Dense>

// for Mesh_3
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/gl_draw/drawing_helper.h>

#include <boost/bind.hpp>

namespace CGAL{
  namespace Anisotropic_mesh_3{

    template<typename PointerToMemberFunction, typename PointerToObject, typename FT, typename Point>
    class Member_function_pointer_to_function_wrapper
    {
    public:
      typedef FT return_type;

      /// Constructor
      Member_function_pointer_to_function_wrapper(PointerToObject po, PointerToMemberFunction pf)
      : pf_(pf) , po_(po) {}

      // Default copy constructor and assignment operator are ok
      /// Destructor
      ~Member_function_pointer_to_function_wrapper() {}

      /// Operator ()
      return_type operator()(const Point& p) const
      {
        return ((*po_).*(pf_))(p);
      }
    private:
      /// Function to wrap
      PointerToMemberFunction pf_;
      PointerToObject po_;
    };  // end class Member_function_pointer_to_function_wrapper


    template<typename K>
    class Constrain_surface_3_implicit :
      public Constrain_surface_3_ex<K, typename Constrain_surface_3<K>::Point_container>
    {
    public:
      typedef Constrain_surface_3_implicit<K>           Self;
      typedef Constrain_surface_3_ex<K, typename Constrain_surface_3<K>::Point_container> Base;
      typedef typename Base::FT                          FT;
      typedef typename Base::Vector_3                    Vector_3;
      typedef typename Base::Point_3                     Point_3;
      typedef typename Base::Oriented_side	         Oriented_side;
       typedef typename Base::Pointset                   Pointset;

      typedef typename Constrain_surface_3<K>::Point_container Point_container;

      // for Mesh_3
      typedef FT (Self::*Function)(const Point_3& p) const;

      typedef Member_function_pointer_to_function_wrapper<Function, const Self*, FT, Point_3> Function_wrapper;
      typedef typename CGAL::Implicit_mesh_domain_3<Function_wrapper, K> Mesh_domain;

      // Triangulation
      typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
      typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

      // Criteria
      typedef typename CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    protected:
      mutable std::set<Point_3> m_poles;
      mutable C3t3 m_c3t3;

      mutable double m_max_curvature;
      mutable double m_min_curvature;
      mutable bool m_cache_max_curvature;
      mutable bool m_cache_min_curvature;

    public:
      virtual FT get_bounding_radius() const = 0;
      virtual typename CGAL::Bbox_3 get_bbox() const = 0;
      virtual std::string name() const { return std::string("Implicit"); }

      virtual FT evaluate(const FT x, const FT y, const FT z) const = 0;
      FT implicit_function(const Point_3& p) const
      {
        return evaluate(p.x(), p.y(), p.z());
      }

      std::set<Point_3>& poles() { return m_poles; }
      const std::set<Point_3>& poles() const { return m_poles; }

      virtual Oriented_side side_of_constraint(const Point_3 &p) const
      {
        FT w = evaluate(p.x(), p.y(), p.z());
        if (w < 0)
          return CGAL::ON_POSITIVE_SIDE;
        else if (w > 0)
          return CGAL::ON_NEGATIVE_SIDE;
        else
          return CGAL::ON_ORIENTED_BOUNDARY;
      }

      virtual Vector_3 gradient(const Point_3 &p, const FT delta = 1e-5) const
      {
        FT dd = 1. / (2. * delta);
        return dd * Vector_3(
          (evaluate(p.x() + delta, p.y(), p.z()) - evaluate(p.x() - delta, p.y(), p.z())),
          (evaluate(p.x(), p.y() + delta, p.z()) - evaluate(p.x(), p.y() - delta, p.z())),
          (evaluate(p.x(), p.y(), p.z() + delta) - evaluate(p.x(), p.y(), p.z() - delta)));
      }

      //virtual Vector_3 normal(const Point_3 &p, const FT delta = 1e-5) const
      //{
      //  Vector_3 n = gradient(p, delta);
      //  FT inv_len = -1. / sqrt(n*n);
      //  return inv_len * n;
      //}

public:
      FT compute_sq_approximation(const Point_3& p) const
      {
        Vector_3 g = gradient(p);
        FT v = implicit_function(p);
        if(g != CGAL::NULL_VECTOR)
          return v*v / (g*g);
        else
        {
          std::cerr << "Warning : Gradient is 0, approximation is wrong." << std::endl;
          return -1;
        }
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

      int are_zeros(const double& d1, const double& d2, const double& d3,
                    bool& z1, bool& z2, bool& z3) const
      {
        int count = 0;
        z1 = (std::abs(d1) < ZERO_EIGENVALUE);  if(z1) count++;
        z2 = (std::abs(d2) < ZERO_EIGENVALUE);  if(z2) count++;
        z3 = (std::abs(d3) < ZERO_EIGENVALUE);  if(z3) count++;
        return count;
      }

      void tensor_frame(const Point_3 &p,
                        Vector_3 &e0, //unit normal
                        Vector_3 &e1, //unit eigenvector
                        Vector_3 &e2, //unit eigenvector
                        double& v1, //eigenvalue corresponding to e1
                        double& v2, //eigenvalue corresponding to e2
                        const FT delta = 1e-5) const
      {
         tensor_frame_eigen(p, e0, e1, e2, v1, v2, delta);
      }

      void tensor_frame_eigen(const Point_3 &p,
                              Vector_3 &e0, //unit normal
                              Vector_3 &e1, //unit eigenvector
                              Vector_3 &e2, //unit eigenvector
                              double& v1, //eigenvalue corresponding to e1
                              double& v2, //eigenvalue corresponding to e2
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
          e0 = get_vector(normal);

          bool z0, z1, z2;
          double ev0 = std::abs(std::real(vals[min_id])); //should be the smallest
          double ev1 = std::abs(std::real(vals[(min_id + 1) % 3]));
          double ev2 = std::abs(std::real(vals[(min_id + 2) % 3]));
          int zeros = are_zeros(ev0, ev1, ev2, z0, z1, z2);

          if(zeros == 1) //it is ev0
          {
            if(!z0) std::cout << "Error1 : see tensor_frame, zeros==1\n";
            Vector_3 normal_test = get_eigenvector(vecs.col(min_id));
            if(!are_equal(e0, normal_test))
             std::cout << "Error2 : see tensor_frame, zeros==1\n";

            v1 = ev1;
            v2 = ev2;
            e1 = get_eigenvector(vecs.col((min_id + 1) % 3));
            e2 = get_eigenvector(vecs.col((min_id + 2) % 3));

            //make sure the vectors form an orthogonal matrix
            //when eigenvalues are the same, vectors returned by Eigen are not
            //necessarily orthogonal
            if(std::abs(e1*e2) > ZERO_DOT_PRODUCT)
              e2 = normalize(CGAL::cross_product(e0,e1));
          }
          else if(zeros == 2)
          {
            int id = -1; // the id of the non-zero eigenvalue
            if(z0 && z1)      id = ((min_id + 2) % 3); //d2
            else if(z0 && z2) id = ((min_id + 1) % 3); //d1
            else if(z1 && z2) id = min_id;//d0
            else std::cerr << "Error1 : see tensor_frame, zeros==2\n";

            //the non-zero one
            v1 = std::abs(std::real(vals[id]));
            e1 = get_eigenvector(vecs.col(id));
            //the last one : choose the vector which is not // to the normal
            Vector_3 normal_test_1 = get_eigenvector(vecs.col((id + 1) % 3));
            Vector_3 normal_test_2 = get_eigenvector(vecs.col((id + 2) % 3));
            if(are_equal(normal_test_1, e0))
            {
              v2 = std::abs(std::real(vals[(min_id + 2) % 3]));
              e2 = normal_test_2;
            }
            else if(are_equal(normal_test_2, e0))
            {
              v2 = std::abs(std::real(vals[(min_id + 1) % 3]));
              e2 = normal_test_1;
            }
            else
            {
              v2 = 0.;
              e2 = normalize(CGAL::cross_product(e0,e1));
              std::cerr << "Error2 : see tensor_frame, zeros==2\n";
            }
          }
          else if(zeros == 3)
          {
            orthogonal_vectors(e0, e1, e2);
            e1 = normalize(e1);
            e2 = normalize(e2);
            v1 = 0.;
            v2 = 0.;
          }
          else //zeros == 0
            std::cerr << "Error : see tensor_frame, zeros==0\n";
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
        Vector_3 get_eigenvector(const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType::ConstColXpr& v) const
        {
          return Vector_3(std::real(v[0]), std::real(v[1]), std::real(v[2]));
        }
        Vector_3 normalize(const Vector_3& v) const
        {
          return std::sqrt(1./(v*v)) * v;
        }
        void orthogonal_vectors(const Vector_3& u,
                                Vector_3& v,
                                Vector_3& w) const
        {
          typename K::Line_3 support(CGAL::ORIGIN, CGAL::ORIGIN + u);
          typename K::Plane_3 p = support.perpendicular_plane(CGAL::ORIGIN);
          v = p.base1();
          w = p.base2();
        }

      virtual std::set<Point_3>& compute_poles() const
      {
        Function fct = (Function)(&Self::implicit_function);
        FT r = this->get_bounding_radius();
        Function_wrapper fw(this, fct);
        Mesh_domain domain(fw, typename K::Sphere_3(CGAL::ORIGIN, r*r));
        Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.,
                               CGAL::parameters::facet_size = r * 0.05,//1,//05,
                               CGAL::parameters::facet_distance = r * 0.01);
                               // cell criteria are ignored
        // run Mesh_3
        m_c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                         CGAL::parameters::no_perturb(),
                                         CGAL::parameters::no_exude());
#ifdef ANISO_OUTPUT_MESH_FOR_POLES
        std::ofstream out("mesh_3_temp.mesh");
        m_c3t3.output_to_medit(out);
#endif
        m_poles.clear();
        compute_triangulation_poles(m_c3t3,
          std::inserter(m_poles, m_poles.end()));
        return m_poles;
      }

      virtual Pointset get_surface_points(unsigned int nb) const
      {
        std::vector<Point_3> all_points;
        typename C3t3::Triangulation::Finite_vertices_iterator v;
        for(v = m_c3t3.triangulation().finite_vertices_begin();
            v != m_c3t3.triangulation().finite_vertices_end();
            ++v)
        {
          if(m_c3t3.in_dimension(v) == 1 || m_c3t3.in_dimension(v) == 2)
            all_points.push_back(v->point());
        }
        std::random_shuffle(all_points.begin(), all_points.end());

        std::size_t n = (std::min)(std::size_t(nb), all_points.size());
        return Pointset(all_points.begin(), (all_points.begin() + n));
      }

      virtual double global_max_curvature() const
      {
        if(m_cache_max_curvature)
          return m_max_curvature;

        m_max_curvature = 0.;
        m_min_curvature = DBL_MAX;
        typename C3t3::Triangulation::Finite_vertices_iterator v;
        for(v = m_c3t3.triangulation().finite_vertices_begin();
            v != m_c3t3.triangulation().finite_vertices_end();
            ++v)
        {
          if(m_c3t3.in_dimension(v) == 1 || m_c3t3.in_dimension(v) == 2)
          {
            double minc, maxc;
            min_max_curvatures(v->point(), minc, maxc);
            m_max_curvature = (std::max)(m_max_curvature, maxc);
            m_min_curvature = (std::min)(m_min_curvature, minc);
          }
        }
        m_cache_max_curvature = true;
        m_cache_min_curvature = true;
        return m_max_curvature;
      }

      virtual double global_min_curvature() const
      {
        if(m_cache_min_curvature)
          return m_min_curvature;

        global_max_curvature(); //computes both max and min
        return m_min_curvature;
      }

      void min_max_curvatures(const Point_3& p,
                              double& minc,
                              double& maxc) const
      {
        Vector_3 n, e1, e2;
        double v1, v2;
        tensor_frame(p, n, e1, e2, v1, v2, 1e-5);
        minc = (std::min)(v1, v2);
        maxc = (std::max)(v1, v2);
      }

      void gl_draw_intermediate_mesh_3(const Plane_3& plane) const
      {
        gl_draw_c3t3<C3t3, Plane_3>(m_c3t3, plane);
      }

      virtual Constrain_surface_3_implicit* clone() const = 0; // Uses the copy constructor

      Constrain_surface_3_implicit()
        : Base(),
          m_poles(),
          m_cache_max_curvature(false),
          m_cache_min_curvature(false) {}

      virtual ~Constrain_surface_3_implicit(){}
    };


    template<typename K>
    class Constrain_surface_3_implicit_with_bounding_sphere :
      public Constrain_surface_3_implicit<K> {

    public:
      typedef Constrain_surface_3_implicit<K> Base;
      typedef typename Base::FT                                FT;
      typedef typename Base::Vector_3                          Vector_3;
      typedef typename Base::Point_3                           Point_3;
      typedef typename Base::Oriented_side		       Oriented_side;

      typedef typename std::vector<typename K::Point_3> Point_container;

    public:
      virtual FT get_bounding_radius() const = 0;
      virtual FT evaluate(const FT x, const FT y, const FT z) const = 0;

      virtual Oriented_side side_of_constraint(const Point_3 &p) const {
        FT br = get_bounding_radius();
        if (p.x() * p.x() + p.y() * p.y() + p.z() * p.z() > br * br)
          return CGAL::ON_NEGATIVE_SIDE;
        FT w = evaluate(p.x(), p.y(), p.z());
        if (w < 0)
          return CGAL::ON_POSITIVE_SIDE;
        else if (w > 0)
          return CGAL::ON_NEGATIVE_SIDE;
        else
          return CGAL::ON_ORIENTED_BOUNDARY;
      }

      Constrain_surface_3_implicit_with_bounding_sphere() { }
      ~Constrain_surface_3_implicit_with_bounding_sphere() { };
    };

  }
}

#endif
