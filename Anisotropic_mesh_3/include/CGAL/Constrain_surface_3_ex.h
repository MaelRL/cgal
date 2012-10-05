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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_EX_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_EX_H

#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Random.h>
#include <CGAL/internal/Intersections_3/Bbox_3_Line_3_do_intersect.h>
#include <CGAL/bbox_intersection_3.h>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

    template<typename K, 
             typename Point_container = typename Constrain_surface_3<K>::Point_container>
    class Constrain_surface_3_ex : 
      public Constrain_surface_3<K> 
    {
    
    public:
      typedef Constrain_surface_3<K> Base;

      typedef typename Base::FT            FT;
      typedef typename Base::Vector_3      Vector_3;
      typedef typename Base::Point_3       Point_3;
      typedef typename Base::Object_3      Object_3;
      typedef typename Base::Segment_3     Segment_3;
      typedef typename Base::Ray_3         Ray_3;
      typedef typename Base::Line_3        Line_3;
      typedef typename Base::Oriented_side Oriented_side;
      typedef typename Base::EdgeList      EdgeList;
      typedef typename EdgeList::iterator  EdgeIterator;
      //typedef typename Base::Point_container Point_container;
      typedef typename Base::Pointset        Pointset;
      
    public:
      EdgeList edges;

    public:
      virtual FT get_bounding_radius() const = 0;
      virtual Oriented_side side_of_constraint(const Point_3 &p) const = 0;
      virtual typename CGAL::Bbox_3 get_bbox() const = 0;
      virtual std::set<Point_3>& compute_poles() const = 0;
      virtual Pointset get_surface_points(unsigned int nb) const = 0;
      virtual std::string name() const  = 0;

      virtual double global_max_curvature() const = 0;
      virtual double global_min_curvature() const = 0;
            
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
          return CGAL::intersection_bl<typename Base::Kernel>
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
      Object_3 intersection(const Point_3& p0, const Point_3& p1) const
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
        FT tol = 1e-16;
        while (sqdist > tol) 
        {
          Oriented_side mv = side_of_constraint(mp);
          if (mv == lv) 
          {
            lp = mp;
            mp = CGAL::midpoint(lp, rp);
          } 
          else if (mv == rv) 
          {
            rp = mp;
            mp = CGAL::midpoint(lp, rp);
          } 
          else 
            return make_object(mp);         
          sqdist *= 0.25;
        }
        return make_object(mp);
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

      EdgeIterator edge_begin() { return edges.begin(); }
      EdgeIterator edge_end() { return edges.end(); }
      void edge_split(EdgeIterator edge, const Point_3 & c) { }
            
      

      Constrain_surface_3_ex() : edges() { }
      ~Constrain_surface_3_ex() { };
    };
  }
}

#endif
