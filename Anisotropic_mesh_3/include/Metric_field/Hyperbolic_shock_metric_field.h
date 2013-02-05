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

#ifndef CGAL_ANISOTROPIC_MESH_3_HYPERBOLIC_SHOCK_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_HYPERBOLIC_SHOCK_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Hyperbolic_shock_metric_field : public Metric_field<K> {
public:
  typedef Metric_field<K>                    Base;
  typedef typename Base::FT                  FT;
  typedef typename Base::Metric              Metric;
  typedef typename Base::Point_3             Point_3;
  typedef typename std::vector<std::pair<Point_3, FT> > ParticleList;

public:
  FT delta;

public:
  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type:   hyperbolic shock" << std::endl;
    fx << "epsilon: " << this->epsilon << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const
  {
    FT x = p.x();
    FT y = p.y();
    FT tanhder = tanh((2.0 * x - sin(5.0 * y)) / delta);
    tanhder = (1.0 - tanhder * tanhder) / delta;
    FT x1 = 2. * tanhder + 3.0 * x * x + y * y;
    FT y1 = -tanhder * cos(5.0 * y) * 5.0 + 2.0 * x * y;

    FT r = sqrt(x1 * x1 + y1 * y1);
    FT x2 = -y1 / r;
    FT y2 = x1 / r;
    FT l1 = std::sqrt(x1*x1 + y1*y1);
    FT l2 = std::sqrt(x2*x2 + y2*y2);

    Vector_3 v1 = (1./l1) * Vector_3(x1, y1, 0);
    Vector_3 v2 = (1./l2) * Vector_3(x2, y2, 0);
    Vector_3 v3(0, 0, 1.);
    return Metric(v1, v2, v3, l1, l2, 1., this->epsilon);
  }

  Hyperbolic_shock_metric_field(FT delta_, FT epsilon_ = 0.125)
    :Metric_field<K>(epsilon_), delta(delta_) { }
};

#endif
