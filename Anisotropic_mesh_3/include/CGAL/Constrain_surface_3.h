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

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

    template<typename K>
    class Constrain_surface_3 
    {
    public:
      typedef K                                    Kernel;
      typedef typename K::FT                       FT;
      typedef typename K::Vector_3                 Vector_3;
      typedef typename K::Point_3                  Point_3;
      typedef typename K::Object_3                 Object_3;
      typedef typename K::Segment_3                Segment_3;
      typedef typename K::Ray_3                    Ray_3;
      typedef typename K::Line_3                   Line_3;
      typedef std::pair<Point_3, Point_3>          Edge;
      typedef std::vector<Edge>                    EdgeList;
      typedef typename EdgeList::iterator	   EdgeIterator;
      typedef typename CGAL::Oriented_side         Oriented_side;  

      typedef typename std::vector<Point_3>        Point_container;
      
      typedef int                                  Surface_patch_index;
      typedef int                                  Index;

      typedef Point_container                      Pointset;

      virtual Object_3 intersection(const Object_3 &obj) const = 0;
      virtual Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const = 0;
      virtual FT get_bounding_radius() const = 0;
      virtual Oriented_side side_of_constraint(const Point_3 &p) const = 0;
      virtual Point_container initial_points(const int nb = 8) const = 0;
      virtual EdgeIterator edge_begin() = 0;
      virtual EdgeIterator edge_end() = 0;
      virtual typename CGAL::Bbox_3 get_bbox() const = 0;
      virtual std::string name() const  = 0;

      virtual std::set<Point_3>& compute_poles() const = 0;
      virtual Pointset get_surface_points(unsigned int nb) const = 0;

      virtual FT compute_sq_approximation(const Point_3& p) const = 0;


      virtual double global_max_curvature() const = 0;
      virtual double global_min_curvature() const = 0;

      virtual void edge_split(EdgeIterator edge, const Point_3 & c) = 0;
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

      virtual ~Constrain_surface_3(){};
    };
  }
}

#endif //CONSTRAIN_SURFACE_3_H
