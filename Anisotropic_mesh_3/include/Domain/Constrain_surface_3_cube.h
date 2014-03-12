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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CUBE_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CUBE_H

#include <CGAL/Constrain_surface_3.h>

#define LINEAR_BLEND(p0,p1,t)  Point_3((p0).x() * (1 - t) + (p1).x() * t, \
                                       (p0).y() * (1 - t) + (p1).y() * t, \
                                       (p0).z() * (1 - t) + (p1).z() * t)

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_cube : public Constrain_surface_3<K, Point_container>
{
public:
  typedef typename K::Point_3           Point_3;
  typedef typename K::Object_3          Object_3;
  typedef typename K::Segment_3         Segment_3;
  typedef typename K::Ray_3             Ray_3;
  typedef typename K::FT                FT;
  typedef typename K::Vector_3          Vector_3;
  typedef typename CGAL::Oriented_side  Oriented_side;

public:
  FT radius;

protected:
  Object_3 intersection_of_ray(const Ray_3 &ray) const
  {
    return intersection(ray.source(), ray.source() +
      ray.to_vector() * ((radius * 3.5) / std::sqrt(ray.to_vector() * ray.to_vector())));
  }

  Object_3 intersection_of_segment(const Segment_3 &seg) const
  {
    return intersection(seg.source(), seg.target());
  }

public:
  Object_3 intersection(const Point_3 &p0, const Point_3 &p1, const Point_3& ref) const
  {
    /*
    FT min_t = 2.0;
    if (p1.x() - p0.x() != 0.0)
    {
      FT t = (-radius - p0.x()) / (p1.x() - p0.x());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.y() >= -radius) && (p.y() <= radius) &&
          (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
      t = (radius - p0.x()) / (p1.x() - p0.x());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.y() >= -radius) && (p.y() <= radius) &&
          (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
    }
    if (p1.y() - p0.y() != 0.0)
    {
      FT t = (-radius - p0.y()) / (p1.y() - p0.y());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.x() >= -radius) && (p.x() <= radius) &&
          (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
        t = (radius - p0.y()) / (p1.y() - p0.y());
        if ((t >= 0) && (t <= 1.0))
        {
          Point_3 p = LINEAR_BLEND(p0, p1, t);
          if ((p.x() >= -radius) && (p.x() <= radius) &&
           (p.z() >= -radius) && (p.z() <= radius))
           if (t < min_t)
             min_t = t;
        }
    }
    if (p1.z() - p0.z() != 0.0)
    {
      FT t = (-radius - p0.z()) / (p1.z() - p0.z());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.x() >= -radius) && (p.x() <= radius) &&
            (p.y() >= -radius) && (p.y() <= radius))
          if (t < min_t)
            min_t = t;
      }
      t = (radius - p0.z()) / (p1.z() - p0.z());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.x() >= -radius) && (p.x() <= radius) &&
            (p.y() >= -radius) && (p.y() <= radius))
          if (t < min_t)
            min_t = t;
      }
    }
    if (min_t > 1.0)
      return Object_3();
    else
      return make_object(LINEAR_BLEND(p0, p1, min_t));
    }*/

    Point_3 lp = p0, rp = p1, mp = CGAL::midpoint(p0, p1);
    Oriented_side lv = side_of_constraint(lp);
    Oriented_side rv = side_of_constraint(rp);
    if (lv == CGAL::ON_ORIENTED_BOUNDARY)
      return make_object(lp);
    if (rv == CGAL::ON_ORIENTED_BOUNDARY)
      return make_object(rp);
    if (lv == rv)
      return Object_3();

    Vector_3 r_l(rp, lp);
    FT sqdist = r_l.squared_length();
    FT error_bound = 1e-5;
    FT tol = radius*radius*error_bound*error_bound;
    while (sqdist > tol)
    {
      Oriented_side mv = side_of_constraint(mp);
      if(rv == mv)
      {
        rp = mp;
        mp = CGAL::midpoint(lp, rp);
      }
      else
      {
        lp = mp;
        mp = CGAL::midpoint(lp, rp);
      }
      sqdist *= 0.25;
    }
    return make_object(mp);
  }

  Object_3 intersection(const Object_3 &obj) const
  {
    Segment_3 seg;
    Ray_3 ray;
    if (CGAL::assign(seg, obj))
      return intersection_of_segment(seg);
    else if (CGAL::assign(ray, obj))
      return intersection_of_ray(ray);
    else
      return Object_3();
  }

  FT get_bounding_radius() const { return radius * 2.0; }
  std::string name() const { return std::string("Cube"); }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT rr = 1.1*radius;
    return CGAL::Bbox_3(-rr, -rr, -rr, rr, rr, rr);
  }

  Oriented_side side_of_constraint(const Point_3 &p) const
  {
    if ((p.x() > radius) || (p.x() < -radius) ||
        (p.y() > radius) || (p.y() < -radius) ||
        (p.z() > radius) || (p.z() < -radius))
      return CGAL::ON_NEGATIVE_SIDE;
    else if ((p.x() > -radius) && (p.x() < radius) &&
             (p.y() > -radius) && (p.y() < radius) &&
             (p.z() > -radius) && (p.z() < radius))
      return CGAL::ON_POSITIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  void compute_poles(std::set<Point_3>& poles) const
  {
    poles.insert(Point_3(0.,0.,0.));
  }

  Point_container initial_points() const
  {
    Point_container points;
    points.push_back(Point_3(-radius, -radius, -radius));
    points.push_back(Point_3(+radius, -radius, -radius));
    points.push_back(Point_3(-radius, +radius, -radius));
    points.push_back(Point_3(+radius, +radius, -radius));
    points.push_back(Point_3(-radius, -radius, +radius));
    points.push_back(Point_3(+radius, -radius, +radius));
    points.push_back(Point_3(-radius, +radius, +radius));
    points.push_back(Point_3(+radius, +radius, +radius));
    points.push_back(Point_3(0, 0, 0));

    points.push_back(Point_3(-radius, 0, 0));
    points.push_back(Point_3(+radius, 0, 0));
    points.push_back(Point_3(0, -radius, 0));
    points.push_back(Point_3(0, +radius, 0));
    points.push_back(Point_3(0, 0, -radius));
    points.push_back(Point_3(0, 0, +radius));

    points.push_back(Point_3(0, -radius, -radius));
    points.push_back(Point_3(0, -radius, +radius));
    points.push_back(Point_3(0, +radius, -radius));
    points.push_back(Point_3(0, +radius, +radius));
    points.push_back(Point_3(-radius, 0, -radius));
    points.push_back(Point_3(-radius, 0, +radius));
    points.push_back(Point_3(+radius, 0, -radius));
    points.push_back(Point_3(+radius, 0, +radius));
    points.push_back(Point_3(-radius, -radius, 0));
    points.push_back(Point_3(-radius, +radius, 0));
    points.push_back(Point_3(+radius, -radius, 0));
    points.push_back(Point_3(+radius, +radius, 0));

    return points;
  }

  Point_container get_surface_points(unsigned int, double) const
  {
    return initial_points();
  }

  Point_3 project(const Point_3 &p) const
  {
    Point_3 q = p;
    if (q.x() < -radius) q = Point_3(-radius, q.y(), q.z());
    if (q.y() < -radius) q = Point_3(q.x(), -radius, q.z());
    if (q.z() < -radius) q = Point_3(q.x(), q.y(), -radius);
    if (q.x() > radius) q = Point_3(radius, q.y(), q.z());
    if (q.y() > radius) q = Point_3(q.x(), radius, q.z());
    if (q.z() > radius) q = Point_3(q.x(), q.y(), radius);
    FT snap_dist = DBL_MAX, sd;
    int snap_type = 0;
    if ((sd = q.x() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 0; }
    if ((sd = q.y() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 1; }
    if ((sd = q.z() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 2; }
    if ((sd = radius - q.x()) < snap_dist) { snap_dist = sd; snap_type = 3; }
    if ((sd = radius - q.y()) < snap_dist) { snap_dist = sd; snap_type = 4; }
    if ((sd = radius - q.z()) < snap_dist) { snap_dist = sd; snap_type = 5; }
    switch (snap_type) {
      case 0:  q = Point_3(-radius, q.y(), q.z()); break;
      case 1:  q = Point_3(q.x(), -radius, q.z()); break;
      case 2:  q = Point_3(q.x(), q.y(), -radius); break;
      case 3:  q = Point_3(radius, q.y(), q.z()); break;
      case 4:  q = Point_3(q.x(), radius, q.z()); break;
      case 5:  q = Point_3(q.x(), q.y(), radius); break;
    }
  return q;
  }

  FT compute_sq_approximation(const Point_3 &p) const
  {
    Point_3 pp = project(p);
    return CGAL::sqrt(CGAL::squared_distance(p,pp));
  }

  double global_max_curvature() const { return 0.;}
  double global_min_curvature() const { return 0.;}

  Constrain_surface_3_cube(const FT &radius_) : radius(radius_) { }
  ~Constrain_surface_3_cube() { }
};

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_flat_cube : public Constrain_surface_3<K, Point_container>
{ //this needs changes, the flat coeff has to be a member instead of hard coded and
  //there should be a coeff for every direction
public:
  typedef typename K::Point_3           Point_3;
  typedef typename K::Object_3          Object_3;
  typedef typename K::Segment_3         Segment_3;
  typedef typename K::Ray_3             Ray_3;
  typedef typename K::FT                FT;
  typedef typename CGAL::Oriented_side  Oriented_side;

public:
  FT radius;

protected:
  Object_3 intersection_of_ray(const Ray_3 &ray) const
  {
    return intersection(ray.source(), ray.source() +
      ray.to_vector() * ((radius * 3.5) / std::sqrt(ray.to_vector() * ray.to_vector())));
  }

  Object_3 intersection_of_segment(const Segment_3 &seg) const
  {
    return intersection(seg.source(), seg.target());
  }

  Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const
  {
    FT min_t = 2.0;
    if (p1.x() - p0.x() != 0.0)
    {
      FT t = (-radius - p0.x()) / (p1.x() - p0.x());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.y() >= -radius) && (p.y() <= radius) &&
          (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
      t = (radius - p0.x()) / (p1.x() - p0.x());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.y() >= -radius) && (p.y() <= radius) &&
          (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
    }
    if (p1.y() - p0.y() != 0.0)
    {
      FT t = (-radius - p0.y()) / (p1.y() - p0.y());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.x() >= -radius) && (p.x() <= radius) &&
            (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
      t = (radius - p0.y()) / (p1.y() - p0.y());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.x() >= -radius) && (p.x() <= radius) &&
          (p.z() >= -radius) && (p.z() <= radius))
          if (t < min_t)
            min_t = t;
      }
    }
    if (p1.z() - p0.z() != 0.0)
    {
      FT t = (-radius - p0.z()) / (p1.z() - p0.z());
      if ((t >= 0) && (t <= 1.0))
      {
        Point_3 p = LINEAR_BLEND(p0, p1, t);
        if ((p.x() >= -radius) && (p.x() <= radius) &&
          (p.y() >= -radius) && (p.y() <= radius))
          if (t < min_t)
            min_t = t;
      }
        t = (radius - p0.z()) / (p1.z() - p0.z());
        if ((t >= 0) && (t <= 1.0))
        {
          Point_3 p = LINEAR_BLEND(p0, p1, t);
          if ((p.x() >= -radius) && (p.x() <= radius) &&
              (p.y() >= -radius) && (p.y() <= radius))
            if (t < min_t)
              min_t = t;
        }
    }
    if (min_t > 1.0)
      return Object_3();
    else
      return make_object(LINEAR_BLEND(p0, p1, min_t));

    /*
    Point_3 pl = p0, pr = p1;
    Oriented_side vl = side_of_constraint(pl);
    Oriented_side vr = side_of_constraint(pr);
    if (vl == CGAL::ON_ORIENTED_BOUNDARY)
      return make_object(pl);
    if (vr == CGAL::ON_ORIENTED_BOUNDARY)
      return make_object(pr);
    FT dist = CGAL::squared_distance(pl, pr);
    while (true) {
      Point_3 pm = LINEAR_BLEND(pl, pr, 0.5);
      Oriented_side vm = side_of_constraint(pm);
      if (vm == CGAL::ON_ORIENTED_BOUNDARY)
          return make_object(project(pm));
      if (vm == vl)
      {
          pl = pm;
          vl = vm;
      }
      else
      {
          pr = pm;
          vr = vm;
      }
      dist /= 4.0;
      if (dist < 1e-10)
          return make_object(project(pm));
    }*/
  }

public:
  Object_3 intersection(const Object_3 &obj) const
  {
    Segment_3 seg;
    Ray_3 ray;
    if (CGAL::assign(seg, obj))
      return intersection_of_segment(seg);
    else if (CGAL::assign(ray, obj))
      return intersection_of_ray(ray);
    else
      return Object_3();
  }

  FT get_bounding_radius() const { return radius * 2.0; } //wrong if stretch coeff <0.5
  std::string name() const { return std::string("Flat cube"); }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT rr = 1.1*radius;
    return CGAL::Bbox_3(-rr, -rr, -rr, rr, rr, rr);
  }

  Oriented_side side_of_constraint(const Point_3 &p) const
  {
    if ((p.x() > radius) || (p.x() < -radius) ||
        (p.y() > radius) || (p.y() < -radius) ||
        (p.z() > radius / 5.0) || (p.z() < -radius / 5.0))
      return CGAL::ON_NEGATIVE_SIDE;
    else if ((p.x() > -radius) && (p.x() < radius) &&
             (p.y() > -radius) && (p.y() < radius) &&
             (p.z() > -radius / 5.0) && (p.z() < radius / 5.0))
      return CGAL::ON_POSITIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Point_container initial_points() const
  {
    Point_container points;
    points.push_back(Point_3(-radius, -radius, -radius / 5.0));
    points.push_back(Point_3(+radius, -radius, -radius / 5.0));
    points.push_back(Point_3(-radius, +radius, -radius / 5.0));
    points.push_back(Point_3(+radius, +radius, -radius / 5.0));
    points.push_back(Point_3(-radius, -radius, +radius / 5.0));
    points.push_back(Point_3(+radius, -radius, +radius / 5.0));
    points.push_back(Point_3(-radius, +radius, +radius / 5.0));
    points.push_back(Point_3(+radius, +radius, +radius / 5.0));
    points.push_back(Point_3(0, 0, 0));

    points.push_back(Point_3(-radius, 0, 0));
    points.push_back(Point_3(+radius, 0, 0));
    points.push_back(Point_3(0, -radius, 0));
    points.push_back(Point_3(0, +radius, 0));
    points.push_back(Point_3(0, 0, -radius / 5.0));
    points.push_back(Point_3(0, 0, +radius / 5.0));

    points.push_back(Point_3(0, -radius, -radius / 5.0));
    points.push_back(Point_3(0, -radius, +radius / 5.0));
    points.push_back(Point_3(0, +radius, -radius / 5.0));
    points.push_back(Point_3(0, +radius, +radius / 5.0));
    points.push_back(Point_3(-radius, 0, -radius / 5.0));
    points.push_back(Point_3(-radius, 0, +radius / 5.0));
    points.push_back(Point_3(+radius, 0, -radius / 5.0));
    points.push_back(Point_3(+radius, 0, +radius / 5.0));
    points.push_back(Point_3(-radius, -radius, 0));
    points.push_back(Point_3(-radius, +radius, 0));
    points.push_back(Point_3(+radius, -radius, 0));
    points.push_back(Point_3(+radius, +radius, 0));

    return points;
  }

  Point_container get_surface_points(unsigned int, double)
  {
    return initial_points();
  }

  Point_3 project(const Point_3 &p) const
  {
    Point_3 q = p;
    if (q.x() < -radius) q = Point_3(-radius, q.y(), q.z());
    if (q.y() < -radius) q = Point_3(q.x(), -radius, q.z());
    if (q.z() < -radius / 5.0) q = Point_3(q.x(), q.y(), -radius / 5.0);
    if (q.x() > radius) q = Point_3(radius, q.y(), q.z());
    if (q.y() > radius) q = Point_3(q.x(), radius, q.z());
    if (q.z() > radius / 5.0) q = Point_3(q.x(), q.y(), radius / 5.0);
    FT snap_dist = DBL_MAX, sd;
    int snap_type = 0;
    if ((sd = q.x() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 0; }
    if ((sd = q.y() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 1; }
    if ((sd = q.z() - (-radius / 5.0)) < snap_dist) { snap_dist = sd; snap_type = 2; }
    if ((sd = radius - q.x()) < snap_dist) { snap_dist = sd; snap_type = 3; }
    if ((sd = radius - q.y()) < snap_dist) { snap_dist = sd; snap_type = 4; }
    if ((sd = radius / 5.0 - q.z()) < snap_dist) { snap_dist = sd; snap_type = 5; }
    switch (snap_type)
    {
    case 0:	q = Point_3(-radius, q.y(), q.z()); break;
    case 1:	q = Point_3(q.x(), -radius, q.z()); break;
    case 2:	q = Point_3(q.x(), q.y(), -radius / 5.0); break;
    case 3:	q = Point_3(radius, q.y(), q.z()); break;
    case 4:	q = Point_3(q.x(), radius, q.z()); break;
    case 5:	q = Point_3(q.x(), q.y(), radius / 5.0); break;
    }
    return q;
  }

  Constrain_surface_3_flat_cube(const FT &radius_) : radius(radius_) { }
  ~Constrain_surface_3_flat_cube() { }
};


template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_free_cube : public Constrain_surface_3<K, Point_container>
{
public:
  typedef typename K::Point_3                   Point_3;
  typedef typename K::Object_3                  Object_3;
  typedef typename K::Segment_3                 Segment_3;
  typedef typename K::Ray_3                     Ray_3;
  typedef typename K::FT                        FT;
  typedef typename CGAL::Oriented_side          Oriented_side;

public:
  FT xmin, ymin, zmin;
  FT xmax, ymax, zmax;

protected:
  Object_3 intersection_of_ray(const Ray_3 &ray) const
  {
    throw "not implemented";
  }

  Object_3 intersection_of_segment(const Segment_3 &seg) const
  {
    throw "not implemented";
  }

  Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const
  {
    throw "not implemented";
  }

public:
  Object_3 intersection(const Object_3 &obj) const
  {
    Segment_3 seg;
    Ray_3 ray;
    if (CGAL::assign(seg, obj))
      return intersection_of_segment(seg);
    else if (CGAL::assign(ray, obj))
      return intersection_of_ray(ray);
    else
      return Object_3();
  }

  FT get_bounding_radius() const { return xmax + ymax + zmax; }
  virtual std::string name() const { return std::string("Free cube"); }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT rx = 0.1*(xmax-xmin);
    FT ry = 0.1*(ymax-ymin);
    FT rz = 0.1*(zmax-zmin);
    return CGAL::Bbox_3(xmin-rx, ymin-ry, zmin-rz, xmax+rx, ymax+ry, zmax+rz);
  }

  Oriented_side side_of_constraint(const Point_3 &p) const
  {
    if ((p.x() > xmax) || (p.x() < xmin) ||
        (p.y() > ymax) || (p.y() < ymin) ||
        (p.z() > zmax) || (p.z() < zmin))
      return CGAL::ON_NEGATIVE_SIDE;
    else if ((p.x() > xmin) && (p.x() < xmax) &&
        (p.y() > ymin) && (p.y() < ymax) &&
        (p.z() > zmin) && (p.z() < zmax))
      return CGAL::ON_POSITIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Point_container initial_points() const
  {
    Point_container points;
    int split = 3;
    for (int i = 0; i <= split; i++)
      for (int j = 0; j <= split; j++)
        for (int k = 0; k <= split; k++)
          points.push_back(Point_3( (FT)i * xmax / (FT)split + xmin,
                                    (FT)j * ymax / (FT)split + ymin,
                                    (FT)k * zmax / (FT)split + zmin));
    return points;
  }

  Point_container get_surface_points(unsigned int, double)
  {
    return initial_points();
  }

  Point_3 project(const Point_3 &p) const
  {
    throw "not implemented";
  }

  Constrain_surface_3_free_cube(const FT &xmin_, const FT &ymin_, const FT &zmin_,
                                const FT &xmax_, const FT &ymax_, const FT &zmax_)
    :
    xmin(xmin_), ymin(ymin_), zmin(zmin_),
    xmax(xmax_), ymax(ymax_), zmax(zmax_)
  { }
  ~Constrain_surface_3_free_cube() { }
};

} // Anisotropic_mesh_3
} // CGAL

#undef LINEAR_BLEND

#endif
