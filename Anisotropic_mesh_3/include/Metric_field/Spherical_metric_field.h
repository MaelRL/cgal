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

#ifndef CGAL_ANISOTROPIC_MESH_3_SPHERICAL_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_SPHERICAL_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;


template<typename K>
class Spherical_metric_field : public Metric_field<K> {
public:
  typedef Metric_field<K>        Base;
  typedef typename Base::FT      FT;
  typedef typename Base::Metric  Metric;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Vector_3 Vector_3;

public:
  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type:   spherical" << std::endl;
    fx << "epsilon: " << this->epsilon << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const
  {
    FT x = p.x(), y = p.y(), z = p.z();
    FT r = x * x + y * y;
    FT R = r + z * z;
    r = std::sqrt(r);
    R = std::sqrt(R);

    if (r == 0.)
      return Metric(Vector_3(1, 0, 0),
                    Vector_3(0, 1, 0), 
                    Vector_3(0, 0, (z >= 0.) ? R : (-R)),
                    1., 1., 1., this->epsilon);
    Vector_3 n(x,y,z);
    n = (1./R) * n;
    Vector_3 v1(-y/r, x/r, 0);
    Vector_3 v2 = CGAL::cross_product(n,v1);
    v2 = (1./std::sqrt(v2*v2)) * v2;

    return Metric(n, v1, v2, 1., 1., 1., this->epsilon);
  }

  Spherical_metric_field(const FT& epsilon_)
    : Metric_field<K>(epsilon_) { }
};


#endif
