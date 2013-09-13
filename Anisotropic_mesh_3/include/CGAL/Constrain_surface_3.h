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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_H

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>

#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <vector>
#include <set>

#include <CGAL/Bbox_3.h>
#include <CGAL/helpers/triangulation_helper.h>
#include <CGAL/Random.h>
#include <CGAL/internal/Intersections_3/Bbox_3_Line_3_do_intersect.h>
#include <CGAL/bbox_intersection_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/helpers/colored_polyhedron_output.h>
#include <CGAL/gl_draw/drawing_helper.h>

#include <CGAL/helpers/mpq_helper.h>
#include <CGAL/helpers/mvpq_helper.h>

#include <Eigen/Dense>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

  template <class Refs, class T, class P>
  class Metric_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
  {
    // tag
    std::size_t m_tag;
    Eigen::Matrix3d m_metric;

  public:
    Metric_vertex()  {}
    Metric_vertex(const P& pt)
      : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
        m_tag(0), m_metric(Eigen::Matrix3d::Zero()) {}

    Metric_vertex(const P& pt, const Eigen::Matrix3d m)
      : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
        m_tag(0), m_metric(m) {}

    // tag
    std::size_t& tag() {  return m_tag; }
    const std::size_t& tag() const {  return m_tag; }
    Eigen::Matrix3d& metric() {  return m_metric; }
    const Eigen::Matrix3d& metric() const {  return m_metric; }

    bool is_colored() const {return m_metric != Eigen::Matrix3d::Zero();}
  };

  template <class Refs, class T>
  class Colored_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
  {
    std::size_t m_tag; // numbering
    double m_color;
    int m_contributors;

  public:
    Colored_facet()
      : m_tag(0), m_color(0.), m_contributors(0)  {}
    const std::size_t& tag() const { return m_tag; }
    std::size_t& tag()             { return m_tag; }
    const double& color() const { return m_color; }
    double& color() { return m_color; }
    const int& contributors() const { return m_contributors; }
    int& contributors() { return m_contributors; }
  };

  struct Aniso_items : public CGAL::Polyhedron_items_3
  {
    template<class Refs, class Traits>
    struct Vertex_wrapper
    {
      typedef typename Traits::Point_3 Point;
      typedef Metric_vertex<Refs, CGAL::Tag_true, Point> Vertex;
    };

    template<class Refs, class Traits>
    struct Face_wrapper
    {
      typedef Colored_facet<Refs, CGAL::Tag_true> Face;
    };
  };

    template<typename K,
             typename Point_container = std::vector<typename K::Point_3>,
             typename Colored_polyhedron = CGAL::Polyhedron_3<K, Aniso_items> >
    class Constrain_surface_3
    {

    public:
      typedef typename K::FT                       FT;
      typedef typename K::Vector_3                 Vector_3;
      typedef typename K::Point_3                  Point_3;
      typedef typename K::Object_3                 Object_3;
      typedef typename K::Segment_3                Segment_3;
      typedef typename K::Triangle_3               Triangle_3;
      typedef typename K::Ray_3                    Ray_3;
      typedef typename K::Line_3                   Line_3;
      typedef typename K::Plane_3                  Plane_3;

      typedef std::pair<Point_3, Point_3>          Edge;
      typedef std::vector<Edge>                    EdgeList;
      typedef typename EdgeList::iterator          EdgeIterator;
      typedef typename CGAL::Oriented_side         Oriented_side;
      typedef int                                  Surface_patch_index;
      typedef int                                  Index;

      typedef Point_container                      Pointset;
      typedef Colored_polyhedron                   Colored_poly;

      typedef typename Colored_polyhedron::Face_handle      Face_handle;
      typedef typename Colored_polyhedron::Facet_iterator   Facet_iterator;
      typedef typename Colored_polyhedron::Vertex_iterator  Vertex_iterator;

      typedef CGAL::AABB_polyhedron_triangle_primitive<K, Colored_polyhedron> Primitive;
      typedef CGAL::AABB_traits<K, Primitive>                                 Traits;
      typedef CGAL::AABB_tree<Traits>                                         Tree;
      typedef typename Tree::Object_and_primitive_id                          Object_and_primitive_id;
      typedef typename Tree::Point_and_primitive_id                           Point_and_primitive_id;
      typedef typename Tree::Primitive_id                                     Primitive_id;

    public:
      EdgeList edges;
      mutable Colored_polyhedron m_colored_poly;
      mutable Colored_polyhedron m_colored_poly_mem;
      mutable Tree* m_colored_poly_tree;

    public:
      virtual FT get_bounding_radius() const = 0;
      virtual Oriented_side side_of_constraint(const Point_3 &p) const = 0;
      virtual void build_colored_polyhedron() const = 0;
      virtual typename CGAL::Bbox_3 get_bbox() const = 0;
      virtual void compute_poles(std::set<Point_3>&) const = 0;
      virtual Point_container get_surface_points(unsigned int nb, double facet_distance_coeff = 1e-3) const = 0;
      virtual std::string name() const  = 0;
      virtual FT compute_sq_approximation(const Point_3& p) const = 0;

      virtual double global_max_curvature() const = 0;
      virtual double global_min_curvature() const = 0;

      virtual void gl_draw_intermediate_mesh_3(const Plane_3& plane) const = 0;
      virtual void gl_draw_tiling(const Plane_3& plane, const int star_id) const = 0;

    protected:
      Object_3 intersection_of_ray(const Ray_3 &ray) const
      {
        return intersection(ray.source(), ray.source() +
          ray.to_vector() *
          (get_bounding_radius() * 2.0 / std::sqrt(ray.to_vector() * ray.to_vector())));
      }

      Object_3 intersection_of_line(const Line_3& line) const
      {
        if(CGAL::do_intersect(get_bbox(), line))
        {
          Point_3 p0 = line.point(0);
          Point_3 p1 = line.point(1);
          return CGAL::intersection_bl<typename K::Kernel>
            (get_bbox(), p0.x(),p0.y(),p0.z(), p1.x(),p1.y(),p1.z(),true,true);
        }
        else
          return Object_3();
      }

      Object_3 intersection_of_segment(const Segment_3 &seg) const
      {
        return intersection(seg.source(), seg.target());
      }

    public:
      virtual Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const = 0;

      Object_3 intersection(const Object_3 &obj) const
      {
        Segment_3 seg;
        Ray_3 ray;
        Line_3 line;
        if (CGAL::assign(seg, obj))
          return intersection_of_segment(seg);
        else if (CGAL::assign(ray, obj))
          return intersection_of_ray(ray);
        else if (CGAL::assign(line, obj))
          return intersection_of_line(line);
        else
          return Object_3();
      }

      Point_3 get_random_direction_end_point(const Point_3 &c) const
      {
        static CGAL::Random r(0);
        FT radius = get_bounding_radius() * 2.0;
        while (true)
        {
          FT x = r.get_double(-1.0, 1.0);
          FT y = r.get_double(-1.0, 1.0);
          FT z = r.get_double(-1.0, 1.0);
          FT w = x * x + y * y + z * z;
          if ((w > 1.0) || (w < 0.1))
            continue;
          w = sqrt(w);
          x /= w;
          y /= w;
          z /= w;
          return Point_3(x * radius + c.x(), y * radius + c.y(), z * radius + c.z());
        }
        return Point_3(0.,0.,0.); //should not be reached
      }

      virtual Point_container initial_points(const Point_container &init_points,
                                             const std::vector<Point_3> &seeds,
                                             const FT min_dist,
                                             const int nb = 8) const
      {
        Point_container points = init_points;
        FT squared_min_dist = min_dist * min_dist;

        for (int i = 0; i < (int)seeds.size(); i++)
        {
          Point_3 center = seeds[i];
          int fail_count = 0;
          while (fail_count < 12 && (int)points.size() < nb)
          {
            Point_3 p;
            Object_3 pobj = intersection(center, get_random_direction_end_point(center));
            if (CGAL::assign(p, pobj))
            {
              bool throwaway = false;
              for (int j = 0; j < (int)points.size(); j++)
              {
                Point_3 q = points[j];
                FT sql = (p.x() - q.x()) * (p.x() - q.x()) +
                         (p.y() - q.y()) * (p.y() - q.y()) +
                         (p.z() - q.z()) * (p.z() - q.z());
                if (sql < squared_min_dist)
                {
                  throwaway = true;
                  break; //for j
                }
              }// end for j
              if (throwaway)
                fail_count++;
              else
              {
                fail_count = 0;
                points.push_back(p);
              }
            } //end if assign
            else
            {
              fail_count++;
            } //end else assign
          } // end while count < 12
        } // end for i
        return points;
      }

      EdgeIterator edge_begin() { return edges.begin(); }
      EdgeIterator edge_end() { return edges.end(); }
      void edge_split(EdgeIterator edge, const Point_3 & c) { }

      void edge_split(const Point_3 &p1, const Point_3 &p2, const Point_3 & c)
      {
        EdgeIterator e = edge_begin();
        EdgeIterator eend = edge_end();
        for (; e != eend; e++)
        {
          if ((fabs(e->first.x() - p1.x()) + fabs(e->first.y() - p1.y()) + fabs(e->first.z() - p1.z()) < 1e-4) &&
            (fabs(e->second.x() - p2.x()) + fabs(e->second.y() - p2.y()) + fabs(e->second.z() - p2.z()) < 1e-4))
          {
            edge_split(e, c);
            break;
          }
        }
      }

      template<typename C3T3>
      void output_points(const C3T3& c3t3,
                         std::ofstream& fx/*.pts file*/) const
      {
        typename C3T3::Triangulation::Finite_vertices_iterator v;
        for(v = c3t3.triangulation().finite_vertices_begin();
            v != c3t3.triangulation().finite_vertices_end();
            ++v)
        {
          if(c3t3.in_dimension(v) == 1 || c3t3.in_dimension(v) == 2)
          {
            typename C3T3::Triangulation::Point p = v->point();
            fx << p.x() << " " << p.y() << " " << p.z() << std::endl;
          }
        }
        fx.close();
      }

      void number_colored_poly() const
      {
        std::size_t index = 0;
        for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v, ++index)
          v->tag() = index;

        index = 0;
        for(Facet_iterator f = m_colored_poly.facets_begin(); f != m_colored_poly.facets_end(); ++f, ++index)
          f->tag() = index;
      }

      void build_colored_poly_tree() const
      {
        std::cout << "build and accelerate tree" << std::endl;
        delete m_colored_poly_tree;
        m_colored_poly_tree = new Tree(m_colored_poly.facets_begin(), m_colored_poly.facets_end());
        m_colored_poly_tree->accelerate_distance_queries();
        std::cout << "tree has size : " << m_colored_poly_tree->size() << std::endl;
      }

      void clear_colors() const
      {
        for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
          v->metric() = Eigen::Matrix3d::Zero();

        for(Facet_iterator f = m_colored_poly.facets_begin(); f != m_colored_poly.facets_end(); ++f)
        {
          f->color() = 0.;
          f->contributors() = 0.;
        }
      }

      void color_poly_vertex()
      {
        //todo
      }

      void color_poly_facet(const Point_3& p, FT ratio) const
      {
        Point_and_primitive_id pp = m_colored_poly_tree->closest_point_and_primitive(p);
        Face_handle closest_f = pp.second; // closest primitive id
        closest_f->color() += ratio;
        closest_f->contributors()++;
      }

      void average_vertex_color_contributor() const
      {
        //todo
      }

      void average_facet_color_contributor() const
      {
        std::cout << "computing average on each facet + stats" << std::endl;

        double strongest_color = -1.;
        int colored_facets_count = 0;

        for(Facet_iterator f = m_colored_poly.facets_begin(); f != m_colored_poly.facets_end(); ++f)
        {
          if(f->contributors() > 0)
          {
            colored_facets_count++;
            f->color() /= f->contributors();
            if(f->color() > strongest_color)
              strongest_color = f->color();
          }
        }
        std::cout << colored_facets_count << "/" << m_colored_poly.size_of_facets() << " facets colored" << std::endl;

        std::ofstream out("colored_poly_ini.off");
        File_writer_OFF writer(false);
        output_colored_polyhedron(out, m_colored_poly, writer, strongest_color);
      }

      void spread_vertices_colors() const
      {
        m_colored_poly_mem = Colored_polyhedron(m_colored_poly);
        typedef Colored_modifiable_vertex_priority_queue<Colored_polyhedron>  Cvmpq;
        Cvmpq q(m_colored_poly.size_of_vertices(), typename Cvmpq::Compare(), typename Cvmpq::ID());
        q.initialize_cmvpq(m_colored_poly);
        q.color_all_vertices();
      }

      void spread_facets_colors() const
      {
        m_colored_poly_mem = Colored_polyhedron(m_colored_poly);
        typedef Colored_modifiable_priority_queue<Colored_polyhedron>  Cmpq;
        Cmpq q(m_colored_poly.size_of_facets(), typename Cmpq::Compare(), typename Cmpq::ID());
        q.initialize_cmpq(m_colored_poly);
        q.color_all_facets();
      }

      void get_vertex_metric_from_poly(const Point_3& p, Eigen::Matrix3d& m) const
      {

      }

      void get_facet_color_from_poly(const Point_3& p, FT& ratio) const
      {
        Point_and_primitive_id pp = m_colored_poly_tree->closest_point_and_primitive(p);
        Point_3 closest_point = pp.first;
        Face_handle closest_f = pp.second; // closest primitive id
        ratio = closest_f->color();
      }

      void gl_draw_colored_polyhedron() const
      {
        gl_draw_colored_poly<Colored_polyhedron>(m_colored_poly);
      }

      void gl_draw_colored_polyhedron_nospread() const
      {
        gl_draw_colored_poly<Colored_polyhedron>(m_colored_poly_mem);
      }

      Constrain_surface_3() :
      edges(),
      m_colored_poly(Colored_polyhedron ()),
      m_colored_poly_mem(Colored_polyhedron ()),
      m_colored_poly_tree(NULL)
      { }

      virtual ~Constrain_surface_3()
      {
        std::cout << "You're deleting trees!!" << std::endl;
        delete m_colored_poly_tree;
      }
    };
  } // end namespace Anisotropic_mesh_3
}

#endif //CONSTRAIN_SURFACE_3_H
