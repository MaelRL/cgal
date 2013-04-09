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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CHAIR_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CHAIR_H

#include <CGAL/Constrain_surface_3_implicit.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Constrain_surface_3_chair : public Constrain_surface_3_implicit<K> 
{

public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT               FT;
  typedef typename Base::Point_container  Point_container;
  typedef typename Base::Point_3          Point_3;

public:
  FT a, b, k;

public:
  virtual std::string name() const { return std::string("Implicit chair"); } 

  void set_a(const FT& aa) { a = aa; }
  void set_b(const FT& bb) { b = bb; }
  void set_k(const FT& kk) { k = kk; }
  FT get_a() const { return a; }
  FT get_b() const { return b; }
  FT get_k() const { return k; }
  
  virtual typename CGAL::Bbox_3 get_bbox() const 
  {
    FT r = get_bounding_radius();
    return CGAL::Bbox_3(-r, -r, -r, r, r, r);
  }

  FT get_bounding_radius() const { return 3.0; }

  FT evaluate(const FT x, const FT y, const FT z) const {
    return (x*x+y*y+z*z-a*k*k) * (x*x+y*y+z*z-a*k*k)
      - b * (((z-k)*(z-k)-2*x*x) * ((z+k)*(z+k)-2*y*y));
  }

  virtual Point_container initial_points(const int nb = 8) const {
    Point_container points;
    std::vector<Point_3> seeds;
    seeds.push_back(Point_3(0, 0, 0));
    seeds.push_back(Point_3(0, 0.5, 0.5));
    seeds.push_back(Point_3(0.5, 0, -0.5));
    seeds.push_back(Point_3(-0.5, 0, -0.5));
    seeds.push_back(Point_3(0, -0.5, 0.5));
    return Base::initial_points(points, seeds, 0.2);
  }

  Constrain_surface_3_chair* clone() const //covariant return types
  {
    return new Constrain_surface_3_chair(*this);
  }

  Constrain_surface_3_chair(const FT a_ = 0.8, const FT b_ = 0.4, const FT k_ = 1.0) 
    : a(a_), b(b_), k(k_)    {}
  Constrain_surface_3_chair(const Constrain_surface_3_chair& c)
    : a(c.a), b(c.b), k(c.k) {}
  ~Constrain_surface_3_chair() {}
};

#endif
