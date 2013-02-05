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

#ifndef CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_METRIC_FIELD

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Metric.h>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

    template<typename K, typename KExact = K>// = CGAL::Exact_predicates_exact_constructions_kernel>
    class Metric_field 
    {
    public:
      typedef Metric_base<K, KExact>  Metric;
      typedef typename K::FT          FT;
      typedef typename K::Point_3     Point_3;
      typedef typename K::Vector_3    Vector_3;

    public:
      FT epsilon;

      virtual Metric compute_metric(const Point_3 &p) const 
      {
        return Metric(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 1., 1., 1., epsilon);
      }

      Metric uniform_metric(const Point_3& p) const
      {
        return Metric(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 1., 1., 1., epsilon);
      }

      // this function is used to report the setting of the metric
      virtual void report(typename std::ofstream &fx) const 
      {
        fx << "type: default" << std::endl;
      }

      Metric_field(FT epsilon_ = 1.0):epsilon(epsilon_) { }

    };
  }
}

#endif
