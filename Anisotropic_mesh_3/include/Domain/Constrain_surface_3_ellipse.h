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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_ELLIPSE_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_ELLIPSE_H

#include <CGAL/Constrain_surface_3_implicit.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Constrain_surface_3_ellipse : public Constrain_surface_3_implicit<K> 
{
public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT       FT;
  typedef typename Base::Point_container   Point_container;
  typedef typename Constrain_surface_3_implicit<K>::Point_3 Point_3;
public:
  FT a, b, c;
  FT bounding_radius;
public:
  FT get_a() const { return a; }
  FT get_b() const { return b; }
  FT get_c() const { return c; }

  virtual std::string name() const { return std::string("Implicit ellipse"); } 

  FT get_bounding_radius() const { return bounding_radius; }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    return CGAL::Bbox_3(-a, -b, -c, a, b, c);
  }

  FT evaluate(const FT x, const FT y, const FT z) const 
  {
    return x * x / a + y * y / b + z * z / c - 1.0;
  }

  Point_container initial_points(const int nb = 8) const 
  {
    Point_container points;
    std::vector<Point_3> seeds;
    seeds.push_back(Point_3(0, 0, 0));
    return Base::initial_points(points, seeds, 0.2);
  }

  Constrain_surface_3_ellipse* clone() const //covariant return types
  {
    return new Constrain_surface_3_ellipse(*this);
  }

  Constrain_surface_3_ellipse(const FT a_, const FT b_, const FT c_) : 
  a(a_ * a_), 
    b(b_ * b_), 
    c(c_ * c_) 
  {
    bounding_radius = (std::max)((std::max)(a_, b_), c_) * 1.1;
  }
  Constrain_surface_3_ellipse(const Constrain_surface_3_ellipse& e)
    : a(e.get_a()), b(e.get_b()), c(e.get_c())
  {}
  ~Constrain_surface_3_ellipse() { };
};

#endif
