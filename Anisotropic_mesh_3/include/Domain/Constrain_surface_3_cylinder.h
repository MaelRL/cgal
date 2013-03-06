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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CYLINDER_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CYLINDER_H

#include <vector>
#include <CGAL/Bbox_3.h>
#include <CGAL/Constrain_surface_3_ex.h>
#include <CGAL/Constrain_surface_3_implicit.h>

using namespace CGAL::Anisotropic_mesh_3;



template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_cylinder :  
  public Constrain_surface_3_implicit<K> {
    //public Constrain_surface_3_ex<K, Point_container> {
public:
  typedef Constrain_surface_3_ex<K, Point_container> Base;
  typedef typename Base::FT                          FT;
  typedef typename Base::Point_3                     Point_3;
  typedef typename Base::Oriented_side               Oriented_side;

public:
  FT radius;
  FT height;

public:
  FT get_bounding_radius() const { return 1.1*(2.*radius + height); }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT r = 1.1*radius;
    FT h = 1.1*(height + radius);//1.1*height;
//    return CGAL::Bbox_3(-r, -r, 0., r, r, h);
    return CGAL::Bbox_3(-r, -r, -r, r, r, h);
  }

  FT evaluate(const FT x, const FT y, const FT z) const 
  {   
    FT ev = radius*radius - (x*x + y*y);
    //outside horizontal boundaries
    if(z > height)
      return ev - (z-height)*(z-height);
    else if(z < 0.)
      return ev - z*z;
    else //inside horizontal boundaries
      return ev;
  }

  Oriented_side side_of_constraint(const Point_3 &p) const 
  {
    FT ev = evaluate(p.x(), p.y(), p.z());
    if(ev < 0.)        return CGAL::ON_NEGATIVE_SIDE;
    else if(ev == 0.)  return CGAL::ON_ORIENTED_BOUNDARY;
    else               return CGAL::ON_POSITIVE_SIDE;

    //FT x = p.x(), y = p.y(), z = p.z();
    //FT radius2 = radius * radius;
    //FT r2 = x * x + y * y;
    //if ((z < height) && (z > 0.) && (r2 < radius2))
    //  return CGAL::ON_POSITIVE_SIDE;
    //else if ((z > height) || (z < 0.) || (r2 > radius2))
    //  return CGAL::ON_NEGATIVE_SIDE;
    //else
    //  return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Constrain_surface_3_cylinder* clone() const //covariant return types
  {
    return new Constrain_surface_3_cylinder(*this);
  }

  virtual double global_max_curvature() const
  {
    return 1. / radius;
  }
  virtual double global_min_curvature() const
  {
    return 0.;
  }
  
  Point_container initial_points(const int nb = 8) const
  {
    Point_container points;
    std::vector<Point_3> seeds;
    seeds.push_back(CGAL::ORIGIN);
    return Base::initial_points(points, seeds, 0.3);
  }

  Constrain_surface_3_cylinder(FT radius_ = 1.0, FT height_ = 10.0) 
    : radius(radius_), height(height_) 
    {}
  Constrain_surface_3_cylinder(const Constrain_surface_3_cylinder& e)
    : radius(e.radius), height(e.height)
    {}

  ~Constrain_surface_3_cylinder() { }
};

#undef	PI

#endif
