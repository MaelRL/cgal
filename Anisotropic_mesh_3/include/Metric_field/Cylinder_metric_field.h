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

#ifndef CGAL_ANISOTROPIC_MESH_3_CYLINDER_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_CYLINDER_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;


template<typename K>
class Cylinder_metric_field : public Metric_field<K> {
public:
  typedef Metric_field<K>        Base;
  typedef typename Base::FT      FT;
  typedef typename Base::Metric  Metric;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Vector_3 Vector_3;

public:
  FT lambda;
  FT radius;
  FT height;

public:
  virtual void report(typename std::ofstream &fx) const {
    fx << "type:   cylinder" << std::endl;
    fx << "lambda: " << lambda << std::endl;
    fx << "radius:  " << radius << std::endl;
    fx << "height:  " << height << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const 
  {
//    FT x = p.x(), y = p.y(), z = p.z();
//    FT l = sqrt(x * x + y * y);
//    Vector_3 g(x / l, y / l, 0); //gradient = normal
//    Vector_3 t(-y / l, x / l, 0);//tangent
//    Vector_3 q(0, 0, z*z+1.0);// lambda * (z * z + 1.0));//vertical
////    return Metric(g, t, q);
//    
//    double ng = std::sqrt(g.squared_length());
//    double nt = std::sqrt(t.squared_length());
//    double nq = std::sqrt(q.squared_length());
//
//    return Metric(1./ng * g,
//                  1./nt * t,
//                  1./nq * q,
//                  1./radius,
//                  nt,
//                  nq,
//                  lambda/*epsilon*/);


    Vector_3 n(p.x(), p.y(), 0.);//normal
    Vector_3 v1(0., 0., 1.);//vertical
    Vector_3 v2 = CGAL::cross_product(n, v1);//tangent

    double nn = std::sqrt(n.squared_length());
    double nv1 = std::sqrt(v1.squared_length());
    double nv2 = std::sqrt(v2.squared_length());

    return Metric(1./nn * n, 
                  1./nv1 * v1, 
                  1./nv2 * v2,
      1./radius,//normal -> e_n = max curvature
      0.,//1./(height*height), 
      1./radius, 
      lambda/*epsilon*/);
  }
  Cylinder_metric_field(FT radius_ = 1.0, 
                        FT height_ = 10.,
                        FT lambda_ = 1.0) 
    : radius(radius_), 
      height(height_),
      lambda(lambda_) { }
};


#endif
