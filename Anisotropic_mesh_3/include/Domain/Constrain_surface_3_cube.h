#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CUBE_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CUBE_H

#include <CGAL/Constrain_surface_3.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_cube :
  public Constrain_surface_3<K, Point_container>
{
public:
  typedef typename K::Point_3           Point_3;
  typedef typename K::Object_3          Object_3;
  typedef typename K::Segment_3         Segment_3;
  typedef typename K::Ray_3             Ray_3;
  typedef typename K::FT                FT;
  typedef typename K::Vector_3          Vector_3;
  typedef CGAL::Oriented_side           Oriented_side;

public:
  FT hside; // half side length

protected:
  Object_3 intersection_of_ray(const Ray_3 &ray) const
  {
    return intersection(ray.source(), ray.source() +
      ray.to_vector() * ((hside * 3.5) / std::sqrt(ray.to_vector() * ray.to_vector())));
  }

  Object_3 intersection_of_segment(const Segment_3 &seg) const
  {
    return intersection(seg.source(), seg.target());
  }

public:
  Object_3 intersection(const Point_3 &p0, const Point_3 &p1, const Point_3& ref) const
  {
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
    FT tol = hside*hside*error_bound*error_bound;
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

  FT get_bounding_radius() const { return hside * 2.0; } // r = hside*sqrt(2)
  std::string name() const { return std::string("Cube"); }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT rr = 1.1*hside;
    return CGAL::Bbox_3(-rr, -rr, -rr, rr, rr, rr);
  }

  Oriented_side side_of_constraint(const Point_3 &p) const
  {
    if ((p.x() > hside) || (p.x() < -hside) ||
        (p.y() > hside) || (p.y() < -hside) ||
        (p.z() > hside) || (p.z() < -hside))
      return CGAL::ON_NEGATIVE_SIDE;
    else if ((p.x() > -hside) && (p.x() < hside) &&
             (p.y() > -hside) && (p.y() < hside) &&
             (p.z() > -hside) && (p.z() < hside))
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
    points.push_back(Point_3(-hside, -hside, -hside));
    points.push_back(Point_3(-hside, +hside, -hside));
    points.push_back(Point_3(-hside, -hside, +hside));
    points.push_back(Point_3(-hside, +hside, +hside));
    points.push_back(Point_3(+hside, -hside, -hside));
    points.push_back(Point_3(+hside, +hside, -hside));
    points.push_back(Point_3(+hside, -hside, +hside));
    points.push_back(Point_3(+hside, +hside, +hside));
    points.push_back(Point_3(-hside/2., -hside, -hside));
    points.push_back(Point_3(-hside/4., +hside, -hside));
    points.push_back(Point_3(-hside/2., -hside, +hside));
    points.push_back(Point_3(-hside/2., +hside, +hside));
    points.push_back(Point_3(+hside/16., -hside, -hside));
    points.push_back(Point_3(+hside/2., +hside, -hside));
    points.push_back(Point_3(+hside/4., -hside, +hside));
    points.push_back(Point_3(+hside/2., +hside, +hside));

    points.push_back(Point_3(-hside, 0, 0));
    points.push_back(Point_3(+hside, 0, 0));
    points.push_back(Point_3(+hside/2., 0, 0));
    points.push_back(Point_3(-hside/2., 0, 0));
    points.push_back(Point_3(-hside/4., 0, -hside));
    points.push_back(Point_3(-hside/32., 0, +hside));
    points.push_back(Point_3(-hside, -hside, 0));
    points.push_back(Point_3(-hside, +hside, 0));
    points.push_back(Point_3(-hside/2., -hside, 0));
    points.push_back(Point_3(-hside/8., +hside, 0));
    points.push_back(Point_3(+hside, -hside, 0));
    points.push_back(Point_3(+hside, +hside, 0));
    points.push_back(Point_3(+hside/2., -hside, 0));
    points.push_back(Point_3(+hside/2., +hside, 0));

    return points;
  }

  Point_container get_surface_points(unsigned int, double) const
  {
    return initial_points();
  }

  Point_3 project(const Point_3 &p) const
  {
    FT x = p.x(), y = p.y(), z = p.z();

    //p is outside the cube
    bool b1 = (p.x() > hside && (x = hside));
    bool b2 = (p.y() > hside && (y = hside));
    bool b3 = (p.z() > hside && (z = hside));
    bool b4 = (p.x() < -hside && (x = -hside));
    bool b5 = (p.y() < -hside && (y = -hside));
    bool b6 = (p.z() < -hside && (z = -hside));

    if(b1 || b2 || b3 || b4 || b5 || b6)
      return Point_3(x,y,z);

    //point is inside the cube
    std::vector<Point_3> possible_projs;

    //simple comparisons to reduce the number of distances to be compared to 3
    if(p.x() < 0.)
      possible_projs.push_back(Point_3(-hside, p.y(), p.z()));
    else
      possible_projs.push_back(Point_3(hside, p.y(), p.z()));

    if(p.y() < 0.)
      possible_projs.push_back(Point_3(p.x(), -hside, p.z()));
    else
      possible_projs.push_back(Point_3(p.x(), hside, p.z()));

    if(p.z() < 0.)
      possible_projs.push_back(Point_3(p.x(), p.y(), -hside));
    else
      possible_projs.push_back(Point_3(p.x(), p.y(), hside));

    int id = -1;
    FT min_sq_dist = 1e16;

    for(int i=0; i<3; ++i)
    {
      FT sq_d = CGAL::squared_distance(p, possible_projs[i]);
      if(sq_d < min_sq_dist)
      {
        id = i;
        min_sq_dist = sq_d;
      }
    }

    return possible_projs[id];
  }

  FT compute_sq_approximation(const Point_3 &p) const
  {
    Point_3 pp = project(p);
    return CGAL::sqrt(CGAL::squared_distance(p,pp));
  }

  Constrain_surface_3_cube(const FT &hside_) : hside(hside_) { }
  ~Constrain_surface_3_cube() { }
};

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_free_cube :
  public Constrain_surface_3<K, Point_container>
{
public:
  typedef typename K::FT                        FT;
  typedef typename K::Point_3                   Point_3;
  typedef typename K::Object_3                  Object_3;
  typedef typename K::Segment_3                 Segment_3;
  typedef typename K::Ray_3                     Ray_3;
  typedef typename K::Vector_3                  Vector_3;
  typedef CGAL::Oriented_side                   Oriented_side;

public:
  FT xmin, ymin, zmin;
  FT xmax, ymax, zmax;

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

  Object_3 intersection(const Point_3 &p0, const Point_3 &p1, const Point_3&) const
  {
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
    FT br = get_bounding_radius();
    FT tol = br*br*error_bound*error_bound;
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

  FT get_bounding_radius() const
  {
    return (xmax-xmin) + (ymax-ymin) + (zmax-zmin);
  }
  virtual std::string name() const { return std::string("Free cube"); }

  virtual typename CGAL::Bbox_3 get_bbox() const
  {
    FT rx = 0.05*(xmax-xmin);
    FT ry = 0.05*(ymax-ymin);
    FT rz = 0.05*(zmax-zmin);
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
    points.push_back(Point_3(xmin, ymin, zmin));
    points.push_back(Point_3(xmin, ymin, zmax));
    points.push_back(Point_3(xmin, ymax, zmin));
    points.push_back(Point_3(xmin, ymax, zmax));
    points.push_back(Point_3(xmax, ymin, zmin));
    points.push_back(Point_3(xmax, ymin, zmax));
    points.push_back(Point_3(xmax, ymax, zmin));
    points.push_back(Point_3(xmax, ymax, zmax));

    points.push_back(Point_3(0.5*(xmin+xmax), ymin, zmin));
    points.push_back(Point_3(0.5*(xmin+xmax), ymin, zmax));
    points.push_back(Point_3(0.5*(xmin+xmax), ymax, zmin));
    points.push_back(Point_3(0.5*(xmin+xmax), ymax, zmax));
    points.push_back(Point_3(0.1*xmin+0.75*xmax, ymin, zmin));
    points.push_back(Point_3(0.1*xmin+0.75*xmax, ymin, zmax));
    points.push_back(Point_3(0.1*xmin+0.75*xmax, ymax, zmin));
    points.push_back(Point_3(0.1*xmin+0.75*xmax, ymax, zmax));

    return points;
  }

  Point_container get_surface_points(unsigned int, double) const
  {
    return initial_points();
  }

  void compute_poles(std::set<Point_3>& poles) const
  {
    FT dx = xmax - xmin;
    FT dy = ymax - ymin;
    FT dz = zmax - zmin;

    poles.insert(Point_3(xmin + 0.60 * dx, ymin + 0.81 * dy, zmin + 0.37 * dz));
    poles.insert(Point_3(xmin + 0.58 * dx, ymin + 0.59 * dy, zmin + 0.41 * dz));
    poles.insert(Point_3(xmin + 0.76 * dx, ymin + 0.87 * dy, zmin + 0.52 * dz));
    poles.insert(Point_3(xmin + 0.10 * dx, ymin + 0.11 * dy, zmin + 0.23 * dz));
    poles.insert(Point_3(xmin + 0.22 * dx, ymin + 0.16 * dy, zmin + 0.37 * dz));
    poles.insert(Point_3(xmin + 0.88 * dx, ymin + 0.72 * dy, zmin + 0.11 * dz));
    poles.insert(Point_3(xmin + 0.66 * dx, ymin + 0.11 * dy, zmin + 0.75 * dz));
    poles.insert(Point_3(xmin + 0.16 * dx, ymin + 0.82 * dy, zmin + 0.95 * dz));
    poles.insert(Point_3(xmin + 0.95 * dx, ymin + 0.86 * dy, zmin + 0.45 * dz));
    poles.insert(Point_3(xmin + 0.29 * dx, ymin + 0.60 * dy, zmin + 0.40 * dz));
  }

  Point_3 project(const Point_3 &p) const
  {
    FT x = p.x(), y = p.y(), z = p.z();

    //p is outside the cube
    bool b1 = (p.x() > xmax && (x = xmax));
    bool b2 = (p.y() > ymax && (y = ymax));
    bool b3 = (p.z() > zmax && (z = zmax));
    bool b4 = (p.x() < xmin && (x = xmin));
    bool b5 = (p.y() < ymin && (y = ymin));
    bool b6 = (p.z() < zmin && (z = zmin));

    if(b1 || b2 || b3 || b4 || b5 || b6)
      return Point_3(x,y,z);

    //point is inside the cube
    std::vector<Point_3> possible_projs;

    //simple comparisons to reduce the number of distances to be compared to 3
    if(p.x() < 0.5*(xmin+xmax))
      possible_projs.push_back(Point_3(xmin, p.y(), p.z()));
    else
      possible_projs.push_back(Point_3(xmax, p.y(), p.z()));

    if(p.y() < 0.5*(ymin+ymax))
      possible_projs.push_back(Point_3(p.x(), ymin, p.z()));
    else
      possible_projs.push_back(Point_3(p.x(), ymax, p.z()));

    if(p.z() < 0.5*(zmin+zmax))
      possible_projs.push_back(Point_3(p.x(), p.y(), zmin));
    else
      possible_projs.push_back(Point_3(p.x(), p.y(), zmax));

    int id = -1;
    FT min_sq_dist = 1e16;

    for(int i=0; i<3; ++i)
    {
      FT sq_d = CGAL::squared_distance(p, possible_projs[i]);
      if(sq_d < min_sq_dist)
      {
        id = i;
        min_sq_dist = sq_d;
      }
    }

    return possible_projs[id];
  }

  FT compute_sq_approximation(const Point_3& p) const
  {
    Point_3 pp = project(p);
    return CGAL::sqrt(CGAL::squared_distance(p,pp));
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

#endif
