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
  typedef Constrain_surface_3_implicit<K>  Base;
  typedef typename Base::FT                FT;
  typedef typename Base::Point_3           Point_3;
  typedef typename Base::Point_container   Point_container;
public:
  FT a, b, c;
  FT bounding_radius;
public:
  void set_a(const FT& aa) { a = aa; }
  void set_b(const FT& bb) { b = bb; }
  void set_c(const FT& cc) { c = cc; }
  FT get_a() const { return a; }
  FT get_b() const { return b; }
  FT get_c() const { return c; }

  virtual std::string name() const { return std::string("Implicit ellipse"); } 

  FT get_bounding_radius() const { return bounding_radius; }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT aa = 1.1*a;
    FT bb = 1.1*b;
    FT cc = 1.1*c;
    return CGAL::Bbox_3(-aa, -bb, -cc, aa, bb, cc);
  }

  FT evaluate(const FT x, const FT y, const FT z) const 
  {
    return x*x/(a*a) + y*y/(b*b) + z*z/(c*c) - 1.0;
  }

//  virtual double global_max_curvature() const
//  {
//    std::cout << "using bad curv values" << std::endl;
//    return 1e30; //todo
//  }
//  virtual double global_min_curvature() const
//  {
//    std::cout << "using bad curv values" << std::endl;
//    return -1e30; //todo
//  }

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

  Constrain_surface_3_ellipse(const FT a_ = 2.0, const FT b_ = 1.0, const FT c_ = 1.0) :
    a(a_),
    b(b_),
    c(c_)
  {
    bounding_radius = (std::max)((std::max)(a_, b_), c_) * 1.1;
  }

  Constrain_surface_3_ellipse(const Constrain_surface_3_ellipse& e)
      : a(e.get_a()), b(e.get_b()), c(e.get_c()), bounding_radius(e.get_bounding_radius()) { }

  ~Constrain_surface_3_ellipse() { };
};

#endif
