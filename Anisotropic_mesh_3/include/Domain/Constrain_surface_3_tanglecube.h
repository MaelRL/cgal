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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_TANGLECUBE_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_TANGLECUBE_H

#include <CGAL/Constrain_surface_3_implicit.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Constrain_surface_3_tanglecube : public Constrain_surface_3_implicit<K>
{
public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT                     FT;
  typedef typename Base::Point_container        Point_container;
  typedef typename Base::Point_3                Point_3;

public:
  FT stretch;

  void set_stretch(const FT& s) { stretch = s; }
  FT get_stretch() { return stretch; }

  virtual std::string name() const { return std::string("Implicit tanglecube"); } 

  FT get_bounding_radius() const { return stretch + 8.0; }

  typename CGAL::Bbox_3 get_bbox() const
  {
    FT r = get_bounding_radius();
    return CGAL::Bbox_3(-r, -r, -r, r, r, r);
  }

  FT evaluate(const FT x, const FT y, const FT z) const 
  {
    FT xs = x / stretch;
    return xs * xs * xs * xs - 5 * xs * xs + y * y * y * y - 5 * y * y + z * z * z * z - 5 * z * z + 11.8;
  }

  Point_container initial_points(const int nb = 8) const 
  {
    Point_container points;
    std::vector<Point_3> seeds;
    //seeds.push_back(Point_3(0, 0, 0));
    seeds.push_back(Point_3(-stretch * 1.5, 1.5, 1.5));
    seeds.push_back(Point_3(stretch * 1.5, 1.5, 1.5));
    seeds.push_back(Point_3(-stretch * 1.5, -1.5, 1.5));
    seeds.push_back(Point_3(stretch * 1.5, -1.5, 1.5));
    seeds.push_back(Point_3(-stretch * 1.5, 1.5, -1.5));
    seeds.push_back(Point_3(stretch * 1.5, 1.5, -1.5));
    seeds.push_back(Point_3(-stretch * 1.5, -1.5, -1.5));
    seeds.push_back(Point_3(stretch * 1.5, -1.5, -1.5));
    return Base::initial_points(points, seeds, 0.2 * stretch);
  }

  Constrain_surface_3_tanglecube* clone() const //covariant return types
  {
    return new Constrain_surface_3_tanglecube(*this);
  }

  Constrain_surface_3_tanglecube(const FT stretch_ = 1.0) 
    : stretch(stretch_) {}
  Constrain_surface_3_tanglecube(const Constrain_surface_3_tanglecube& t)
    : stretch(t.stretch) {}
  ~Constrain_surface_3_tanglecube() {};
};

#endif
