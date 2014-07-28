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

#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/gl_draw/drawing_helper.h>
#include <CGAL/helpers/aniso_polyhedron_items.h>

#include <CGAL/Polyhedron_3.h>

#include <Eigen/Dense>

namespace CGAL
{
namespace Anisotropic_mesh_3
{
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

  typedef typename CGAL::Oriented_side         Oriented_side;
  typedef int                                  Subdomain_index;
  typedef int                                  Surface_patch_index;
  typedef int                                  Index;

  typedef Point_container                      Pointset;
  typedef Colored_polyhedron                   Colored_poly;

public:
  virtual FT get_bounding_radius() const = 0;
  virtual Oriented_side side_of_constraint(const Point_3 &p) const = 0;
  virtual typename CGAL::Bbox_3 get_bbox() const = 0;
  virtual void compute_poles(std::set<Point_3>& poles) const = 0;
  virtual Point_container get_surface_points(unsigned int nb, double facet_distance_coeff = 0.01) const = 0;
  virtual std::string name() const  = 0;
  virtual FT compute_sq_approximation(const Point_3& p) const = 0;

  virtual void gl_draw_intermediate_mesh_3(const Plane_3& plane) const
  {
    std::cout << "draw intermediate mesh3: incompatible surface" << std::endl;
    return;
  }

  virtual void build_colored_polyhedron(Colored_polyhedron& poly) const
  {
    std::cout << "build polyehdron: incompatible surface" << std::endl;
    return;
  }

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
  virtual Object_3 intersection(const Point_3 &p0, const Point_3 &p1, const Point_3 &ref) const = 0;

  Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const
  {
    return intersection(p0, p1, p0); //default ref point is p0
  }

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

  Constrain_surface_3(){ }

  virtual ~Constrain_surface_3()
  { }
};

} // Anisotropic_mesh_3
} // CGAL

#endif //CONSTRAIN_SURFACE_3_H
