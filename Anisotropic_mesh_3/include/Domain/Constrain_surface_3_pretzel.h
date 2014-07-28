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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_PRETZEL_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_PRETZEL_H

#include <CGAL/Constrain_surface_3_implicit.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Constrain_surface_3_pretzel : public Constrain_surface_3_implicit<K>
{

public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT               FT;
  typedef typename Base::Point_container  Point_container;
  typedef typename Base::Point_3          Point_3;

public:
  FT c, d;

public:
  virtual std::string name() const { return std::string("Implicit pretzel"); }

  void set_c(const FT& cc) { c = cc; }
  void set_d(const FT& dd) { d = dd; }
  FT get_c() const { return c; }
  FT get_d() const { return d; }
  
  virtual typename CGAL::Bbox_3 get_bbox() const 
  {
    FT r = get_bounding_radius();
    return CGAL::Bbox_3(-r, -r, -r, r, r, r);
  }

  FT get_bounding_radius() const
  {
    FT rc = std::abs(c);
    FT rd = std::abs(d);
    if( rd < 1.)
      rd = 1./d;
    if( rc < 1.)
      rc = 1./c;
    return rc + rd;
  }

  FT evaluate(const FT x, const FT y, const FT z) const
  {
    FT sq_y = y*y;
    FT sq_z = z*z;
    FT temp = ((x-1)*(x-1) + sq_y - c*c)*((x+1)*(x+1) + sq_y - c*c);
    return temp*temp + sq_z - d;
  }

  virtual Point_container initial_points(const int nb = 8) const {
    Point_container points;
    std::vector<Point_3> seeds;
    seeds.push_back(Point_3(0, 0, 0));
    seeds.push_back(Point_3(2, 0, 0));
    seeds.push_back(Point_3(-2, 0, 0));
    return Base::initial_points(points, seeds, 0.2);
  }

  Constrain_surface_3_pretzel* clone() const //covariant return types
  {
    return new Constrain_surface_3_pretzel(*this);
  }

  Constrain_surface_3_pretzel(const FT c_ = 2, const FT d_ = 1)
    : c(c_), d(d_) {}
  Constrain_surface_3_pretzel(const Constrain_surface_3_pretzel& p)
    : c(p.c), d(p.d) {}
  ~Constrain_surface_3_pretzel() {}
};

} //Anisotropic_mesh_3
} //CGAL

#endif
