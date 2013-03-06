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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_PARTICLES_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_PARTICLES_H

#include "../CGAL/Constrain_surface_3_implicit.h"

#define PI 3.1415926535897932384626433832795

using namespace CGAL::Anisotropic_mesh_3;

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_particles : public Constrain_surface_3_implicit<K, Point_container>
{
public:
  typedef Constrain_surface_3_implicit<K, Point_container> Base;
  typedef typename Base::FT                                FT;

  struct Particle
  {
    FT x, y, z, q;
  };
public:
  std::vector<Particle> particles;
  FT level;
  FT bounding_radius;
public:
  void add_particle(FT x, FT y, FT z, FT q)
  {
    Particle p;
    p.x = x; p.y = y; p.z = z; p.q = q;
    particles.push_back(p);
  }

  FT get_bounding_radius() const { return bounding_radius; }

  FT evaluate(const FT x, const FT y, const FT z) const
  {
    std::vector<Particle>::const_iterator i = particles.begin();
    std::vector<Particle>::const_iterator iend = particles.end();
    FT sum = 0.0;
    for (; i != iend; i++) {
      FT dist = ((*i).x - x) * ((*i).x - x) +
                ((*i).y - y) * ((*i).y - y) +
                ((*i).z - z) * ((*i).z - z);
      sum += (*i).q / dist;
    }
    return level - sum;
  }

  Point_container initial_points() const
  {
    Point_container points;
    std::vector<Point_3> seeds;
    std::vector<Particle>::const_iterator i = particles.begin();
    std::vector<Particle>::const_iterator iend = particles.end();
    for (; i != iend; i++) {
      if ((*i).q > 0)
        seeds.push_back(Point_3((*i).x, (*i).y, (*i).z));
    }
    return Base::initial_points(points, seeds, 0.3);
  }

  Constrain_surface_3_particles(FT level_, FT bounding_radius_) :
    level(level_), bounding_radius(bounding_radius_), particles() { }
  ~Constrain_surface_3_particles() { }
};

#undef  PI

#endif
