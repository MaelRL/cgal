// Copyright (c) 2012, 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the Licenxse, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Claudia WERNER

#ifndef CGAL_DELAUNAY_MESH_SPHERE_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_SPHERE_TRAITS_2_H

#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>

#include <CGAL/triangulation_assertions.h>

namespace CGAL {

//projects the computed points (e.g. midpoint, circumcenter) on the sphere
template <class K, class Predicate_>
class Project_on_sphere_adaptor
{
public:
  typedef Predicate_ Predicate;
  typedef typename K::Point_3 Point;
  Project_on_sphere_adaptor(double radius,Point sphere):_radius(radius), _sphere(sphere) { }

  double _radius;
  Point _sphere;

  typedef typename Predicate::result_type                         result_type;

  result_type operator()(const Point& p0, const Point& p1)
  {
    Point p = Predicate()(p0, p1);
    return project(p);
  }

  result_type operator()(const Point& p0, const Point& p1, const Point& p2)
  {
    Point p = Predicate()(p0, p1, p2);
    return project(p);
  }

private:
  Point project(const Point& p) { return project(p.x(), p.y(), p.z()); }
  Point project(double x, double y, double z)
  {
    double dist = x*x + y*y + z*z;
    if(dist == 0)
      return Point(0,0,0);

    double scale = _radius / sqrt(dist);
    return Point(x*scale, y*scale, z*scale);
  }
};

template < class Base >
class Delaunay_mesh_sphere_traits_2
  : public Base
{
public:
  typedef Delaunay_mesh_sphere_traits_2<Base>                 Self;

  typedef typename Base::Point_2                              Point_2;
  typedef typename Base::Kernel                               Kernel;
  typedef typename Kernel::Vector_3                           Vector_2;
  typedef typename Kernel::Compute_area_3                     Compute_area_2;
  typedef typename Kernel::Angle_3                            Angle_2;
  typedef typename Kernel::Compare_distance_3                 Compare_distance_2;
  typedef typename Kernel::Construct_vector_3                 Construct_vector_2;
  typedef typename Kernel::Construct_scaled_vector_3          Construct_scaled_vector_2;
  typedef typename Kernel::Construct_translated_point_3       Construct_translated_point_2;
  typedef typename Kernel::Construct_circumcenter_3           Circumcenter_2;

  typedef Project_on_sphere_adaptor<Kernel, typename Kernel::Construct_midpoint_3> Construct_midpoint_2;
  typedef Project_on_sphere_adaptor<Kernel, typename Kernel::Construct_circumcenter_3> Construct_circumcenter_2;
  typedef boost::true_type                                    requires_test;

  double _radius;
  double _minDistSquared;

  Delaunay_mesh_sphere_traits_2(double radius = 1);

  double radius() { return _radius; }
  void set_radius(double a)
  {
    _radius = a;
    _minDistSquared = (_radius * pow(2, -23)) * (_radius * pow(2, -23));
  }

  Compare_distance_2
  compare_distance_2_object() const
  { return Compare_distance_2(); }

  Compute_area_2
  compute_area_2_object() const
  { return Compute_area_2(); }

  Angle_2
  angle_2_object() const
  { return Angle_2(); }

  Construct_vector_2
  construct_vector_2_object() const
  { return Construct_vector_2(); }

  Construct_scaled_vector_2
  construct_scaled_vector_2_object() const
  { return Construct_scaled_vector_2(); }

  Construct_translated_point_2
  construct_translated_point_2_object() const
  { return Construct_translated_point_2(); }

  Construct_midpoint_2
  construct_midpoint_2_object() const
  { return Construct_midpoint_2(_radius, Base::_sphere); }

  Construct_circumcenter_2
  construct_circumcenter_2_object() const
  { return Construct_circumcenter_2(_radius, Base::_sphere); }

protected:
  Point_2 _sphere;
};

template <class R>
  Delaunay_mesh_sphere_traits_2<R>::
  Delaunay_mesh_sphere_traits_2(double radius)
    : _radius(radius)
  { }

} // namespace CGAL

#endif // CGAL_Reg_TRIANGULATION_SPHERE_TRAITS_2_H
