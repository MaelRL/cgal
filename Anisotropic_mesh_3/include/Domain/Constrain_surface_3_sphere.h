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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_SHPERE_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_SHPERE_H

#include <CGAL/Constrain_surface_3_implicit.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Constrain_surface_3_sphere : public Constrain_surface_3_implicit<K> 
{
public:
  typedef Constrain_surface_3_implicit<K> Base;
  typedef typename Base::FT               FT;
  typedef typename Base::Point_container  Point_container;
  typedef typename Base::Point_3          Point_3;
  typedef typename Base::Pointset         Pointset;

public:
  virtual std::string name() const { return std::string("Implicit sphere"); } 

  virtual FT get_bounding_radius() const { return 1.1; }

  typename CGAL::Bbox_3 get_bbox() const
  {
    return CGAL::Bbox_3(-1.1, -1.1, -1.1, 1.1, 1.1, 1.1);
  }

  FT evaluate(const FT x, const FT y, const FT z) const 
  {
    return x * x + y * y + z * z - 1.0;
  }

  Point_container initial_points(const int nb = 8) const 
  {
    Point_container points;
    std::vector<Point_3> seeds;
    seeds.push_back(Point_3(0, 0, 0));
    return Base::initial_points(points, seeds, 0.2);
  }

  virtual Pointset get_surface_points(unsigned int nb) const
  {
    Point_container ip = initial_points(nb);
    return Pointset(ip.begin(), ip.end());
  }

  virtual std::set<Point_3>& compute_poles() const
  {
    // warning : the sphere is a degenerate case
    this->m_poles.clear();
    this->m_poles.insert(CGAL::ORIGIN);
    this->m_poles.insert(Point_3(0., 0.1, 0.1)); //not a pole : to avoid degeneracies
    this->m_poles.insert(Point_3(0.1, -0.15, -0.05)); //not a pole : to avoid degeneracies
    return this->m_poles;
  }

  Constrain_surface_3_sphere* clone() const // Covariant Return Types
    // see http://www.parashift.com/c++-faq-lite/virtual-functions.html
  { 
    return new Constrain_surface_3_sphere(*this); 
  }

  Constrain_surface_3_sphere() { }
  Constrain_surface_3_sphere(const Constrain_surface_3_sphere& s) {}
  ~Constrain_surface_3_sphere() { }
};

#endif
