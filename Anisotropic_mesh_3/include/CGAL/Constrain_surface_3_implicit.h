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

#include <CGAL/Constrain_surface_3.h>

#include <Eigen/Dense>

// for Mesh_3
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/gl_draw/drawing_helper.h>

#include <CGAL/helpers/metric_helper.h>
#include <CGAL/helpers/c3t3_polyhedron_builder.h>
#include <CGAL/IO/Polyhedron_iostream.h>


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


    template<typename K,
             typename Pt_container = std::vector<typename K::Point_3> >
    class Constrain_surface_3_implicit :
      public Constrain_surface_3<K, Pt_container>
    {
    public:
      typedef Constrain_surface_3_implicit<K, Pt_container> Self;
      typedef Constrain_surface_3<K, Pt_container>          Base;
      typedef typename Base::Object_3                       Object_3;
      typedef typename Base::FT                             FT;
      typedef typename Base::Vector_3                       Vector_3;
      typedef typename Base::Point_3                        Point_3;
      typedef typename Base::Plane_3                        Plane_3;
      typedef typename Base::Oriented_side                  Oriented_side;
      typedef typename Base::Colored_poly                   Colored_polyhedron;

      typedef Pt_container                                  Point_container;

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
      mutable C3t3 m_c3t3;
      FT error_bound;

    public:
      C3t3 c3t3() const { return m_c3t3; }

      virtual FT get_bounding_radius() const = 0;
      virtual typename CGAL::Bbox_3 get_bbox() const = 0;
      virtual std::string name() const { return std::string("Implicit"); }

      virtual FT evaluate(const FT x, const FT y, const FT z) const = 0;
      FT implicit_function(const Point_3& p) const
      {
        return evaluate(p.x(), p.y(), p.z());
      }

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

public:
      Object_3 intersection(const Point_3 &p0, const Point_3 &p1, const Point_3 &ref) const
      {
#ifdef ANISO_IMPLICIT_INTERSECTION_STEPS
        Point_3 lp = p0, rp = p1;
        Oriented_side lv = side_of_constraint(lp);
        Oriented_side rv = side_of_constraint(rp);
        if (lv == CGAL::ON_ORIENTED_BOUNDARY)
          return make_object(lp);
        if (rv == CGAL::ON_ORIENTED_BOUNDARY)
          return make_object(rp);
        if (lv == rv)
          return Object_3();

        Point_3 mem_p = lp;
        Vector_3 l_r(lp, rp);
        int steps = 100, count = 0;
        Vector_3 step_vec = l_r/steps;
        Point_3 step_p = lp + step_vec;

        while (count++ < steps)
        {
          Oriented_side stepv = side_of_constraint(step_p);
          if(lv != stepv)
            return intersection_dichotomy(mem_p, step_p);
          mem_p = step_p;
          step_p = step_p + step_vec;
        }
        return Object_3();
#else
        return intersection_dichotomy(p0, p1);
#endif
      }

      Object_3 intersection_dichotomy(const Point_3& p0, const Point_3& p1) const
      {
        Point_3 lp = p0, rp = p1, mp = CGAL::midpoint(p0, p1);
        Oriented_side lv = side_of_constraint(lp);
        Oriented_side rv = side_of_constraint(rp);
        if (lv == CGAL::ON_ORIENTED_BOUNDARY)
          return make_object(lp);
        if (rv == CGAL::ON_ORIENTED_BOUNDARY)
          return make_object(rp);
        if (lv == rv)
          return Object_3();

        Vector_3 r_l(rp, lp);
        FT sqdist = r_l.squared_length();
        FT r = this->get_bounding_radius();
        FT tol = r*r*error_bound*error_bound;
        while (sqdist > tol)
        {
          Oriented_side mv = side_of_constraint(mp);
          if(rv == mv)
          {
            rp = mp;
            mp = CGAL::midpoint(lp, rp);
          }
          else
          {
            lp = mp;
            mp = CGAL::midpoint(lp, rp);
          }
          sqdist *= 0.25;
        }
        return make_object(mp);
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

      void build_colored_polyhedron(Colored_polyhedron& poly) const
      {
        std::cout << "trying to convert c3t3 to polyhedron and outputing it (I)...";
        Complex_3_in_triangulation_3_polyhedron_builder<C3t3, Colored_polyhedron> builder(m_c3t3);
        poly = Colored_polyhedron();
        poly.delegate(builder);
        std::ofstream out("colored_poly.off");
        out << poly;
        std::cout << "done" << std::endl;
      }

      // Sizing field
      template<typename CSI>
      struct Aniso_sizing_field
      {
        const CSI* const csi;

        typedef typename CSI::FT                 FT;
        typedef typename CSI::Mesh_domain::Index Index;
        FT operator()(const Point_3& p, const int, const Index&) const
        {
          Vector_3 v0, v1, v2;
          double e0, e1, e2;
          csi->tensor_frame(p, v0, v1, v2, e0, e1, e2, 1e-5);
          e1 = std::max(e1, 1e-10);
          e2 = std::max(e2, 1e-10);
          FT ratio = 1./std::sqrt(std::max(e1/e2, e2/e1));
          FT ret = std::max(1e-3, ratio);

          ret *= 0.1*csi->get_bounding_radius();

          //std::cout << "ratio : " << ratio << " " << 1./ratio << std::endl;
          //std::cout << "ret : " << ret << std::endl;
          return ret;
        }

        Aniso_sizing_field(const CSI* const csi_):csi(csi_){}
      };

      virtual Point_container get_surface_points(unsigned int nb,
                                                 double facet_distance_coeff /*= 0.05*/) const
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

        if(std::size_t(nb) > all_points.size())
        {
#ifdef ANISO_VERBOSE
          std::cout << "C3T3 mesh did not contain enough points (" << all_points.size();
          std::cout << "). Generating a new c3t3...";
          std::cout << "Facet distance : " << facet_distance_coeff << std::endl;
#endif
          //not enough points on the c3t3 to fill the initial_points vector. Generating a new denser c3t3.
          Function fct = (Function)(&Self::implicit_function);
          FT r = this->get_bounding_radius();
          Function_wrapper fw(this, fct);
          Mesh_domain domain(fw, typename K::Sphere_3(CGAL::ORIGIN, r*r), error_bound);

          Aniso_sizing_field<Self> size(this);
          Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.,
                                 CGAL::parameters::facet_size = r * 0.05, //size
                                 CGAL::parameters::facet_distance = r * facet_distance_coeff,
                                 CGAL::parameters::facet_topology = MANIFOLD);
                                 // cell criteria are ignored

          // run Mesh_3
          m_c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                           CGAL::parameters::no_perturb(),
                                           CGAL::parameters::no_exude());
#ifdef ANISO_OUTPUT_MESH_FOR_POLES
        std::ofstream out("mesh_3_temp.mesh");
        m_c3t3.output_to_medit(out);
#endif
          return get_surface_points(nb, facet_distance_coeff/3.);
        }
#ifdef ANISO_VERBOSE
        std::cout << "c3t3 has enough points : " << all_points.size() << " faces : " << m_c3t3.number_of_facets_in_complex() << std::endl;
#endif
        std::random_shuffle(all_points.begin(), all_points.end());
        return Point_container(all_points.begin(), (all_points.begin() + nb));
      }

      C3t3 run_mesh_3(const FT& approx, const bool verbose = false) const
      {
        Function fct = (Function)(&Self::implicit_function);
        FT r = this->get_bounding_radius();
        Function_wrapper fw(this, fct);
        Mesh_domain domain(fw, typename K::Sphere_3(CGAL::ORIGIN, r*r), error_bound);
        Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.,
                               CGAL::parameters::facet_size = r * 0.05,
                               CGAL::parameters::facet_distance = approx,
                               CGAL::parameters::facet_topology = MANIFOLD);
                               //no cell criteria          
        C3t3 mesh = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                       CGAL::parameters::no_perturb(),
                                       CGAL::parameters::no_exude());
        if(verbose)
          std::cout << "Mesh has " 
                    << mesh.triangulation().number_of_vertices () << " vertices.\n";
        return mesh;
      }

      virtual void compute_poles(std::set<Point_3>& poles) const
      {
        compute_triangulation_poles(m_c3t3, std::inserter(poles, poles.end()), get_bbox());
      }

      void gl_draw_intermediate_mesh_3(const Plane_3& plane) const
      {
        gl_draw_c3t3<C3t3, Plane_3>(m_c3t3, plane);
      }

      virtual Constrain_surface_3_implicit* clone() const = 0; // Uses the copy constructor

      Constrain_surface_3_implicit(const FT error = 1e-5)
        :
          Base(),
          error_bound(error)
      { }

      virtual ~Constrain_surface_3_implicit(){}
    };

    template<typename K>
    class Constrain_surface_3_implicit_with_bounding_sphere :
      public Constrain_surface_3_implicit<K>
    {

    public:
      typedef Constrain_surface_3_implicit<K>                  Base;
      typedef typename Base::FT                                FT;
      typedef typename Base::Vector_3                          Vector_3;
      typedef typename Base::Point_3                           Point_3;
      typedef typename Base::Oriented_side                     Oriented_side;

      typedef typename std::vector<typename K::Point_3>        Point_container;

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
      ~Constrain_surface_3_implicit_with_bounding_sphere() { }
    };

  }
}

#endif
