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

#ifndef CGAL_ANISOTROPIC_MESH_3_PARTICLE_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_PARTICLE_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Particle_metric_field : public Metric_field<K> {
public:
  typedef typename K::FT                                FT;
  typedef typename Metric_field<K>::Metric              Metric;
  typedef typename Metric_field<K>::Point_3             Point_3;
  typedef typename std::vector<std::pair<Point_3, FT> > ParticleList;
  ParticleList particles;

public:
  void add_particle(const Point_3 &p, const FT &weight)
  {
    particles.push_back(std::pair<Point_3, FT>(p, weight));
  }

  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type:   particle" << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const
  {
    // compute the force
    FT sumx = 0.0, sumy = 0.0, sumz = 0.0;
    FT x = p.x(), y = p.y(), z = p.z();
    ParticleList::const_iterator p_end = particles.end();
    for (ParticleList::const_iterator pt = particles.begin(); pt != p_end; pt++)
    {
      FT px = pt->first.x(), py = pt->first.y(), pz = pt->first.z(), w = pt->second;
      FT r2 = (px - x) * (px - x) + (py - y) * (py - y) + (pz - z) * (pz - z);
      sumx += (px - x) / r2 * w;
      sumy += (py - y) / r2 * w;
      sumz += (pz - z) / r2 * w;
    }
    FT l = sqrt(sumx * sumx + sumy * sumy + sumz * sumz);
    FT x1 = sumx / l / (l + 1.0);
    FT y1 = sumy / l / (l + 1.0);
    FT z1 = sumz / l / (l + 1.0);
    FT x2, y2, z2, x3, y3, z3;
    if ((fabs(x1) < 1e-7) && (fabs(y1) < 1e-7))
    {
      x2 = 1.0; y2 = 0.0; z2 = 0.0;
      x3 = 0.0; y3 = 1.0; z3 = 0.0;
    }
    else
    {
      FT r = sqrt(x1 * x1 + y1 * y1);
      x2 = -y1 / r;
      y2 = x1 / r;
      z2 = 0.0;
      FT R = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
      x3 = -x1 * z1 / r / R;
      y3 = -y1 * z1 / r / R;
      z3 = r / R;
    }
    return Metric(Vector_3(x1, y1, z1), Vector_3(x2, y2, z2), Vector_3(x3, y3, z3));
  }

  Particle_metric_field() : particles() { }
};


#endif
