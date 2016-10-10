#ifndef CGAL_ANISOTROPIC_MESH_2_RECTANGLE_DOMAIN_H
#define CGAL_ANISOTROPIC_MESH_2_RECTANGLE_DOMAIN_H

#include <CGAL/Domain_2.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K,
         typename Point_container = std::vector<typename K::Point_2> >
class Rectangle_domain :
    public Domain_2<K, Point_container>
{
public:
  typedef typename K::Point_2           Point_2;
  typedef typename K::Object_2          Object_2;
  typedef typename K::Segment_2         Segment_2;
  typedef typename K::Ray_2             Ray_2;
  typedef typename K::FT                FT;
  typedef typename K::Vector_2          Vector_2;
  typedef CGAL::Oriented_side           Oriented_side;

public:
  FT hside_x; // half side length along Ox
  FT hside_y; // half side length along Oy
  FT hside;
  Point_2 center;

public:
  FT get_bounding_radius() const { return hside * 2.0; } // r = hside*sqrt(2)
  std::string name() const { return std::string("Rectangle"); }

  virtual typename CGAL::Bbox_2 get_bbox() const
  {
    FT rx = 1.1*hside_x;
    FT ry = 1.1*hside_y;
    FT x = center.x();
    FT y = center.y();
    return CGAL::Bbox_2(x-rx, y-ry, x+rx, y+ry);
  }

  Oriented_side side_of_constraint(const Point_2 &p) const
  {
    FT x = center.x();
    FT y = center.y();
    if ((p.x() > x + hside_x) || (p.x() < x - hside_x) ||
        (p.y() > y + hside_y) || (p.y() < y - hside_y))
      return CGAL::ON_NEGATIVE_SIDE;
    else if ((p.x() > x - hside_x) && (p.x() < x + hside_x) &&
             (p.y() > y - hside_y) && (p.y() < y + hside_y))
      return CGAL::ON_POSITIVE_SIDE;
    else
      return CGAL::ON_ORIENTED_BOUNDARY;
  }

  Point_container initial_points() const
  {
    FT x = center.x();
    FT y = center.y();

    Point_container points;
    points.push_back(Point_2(x-hside_x, y-hside_y));
    points.push_back(Point_2(x+hside_x, y-hside_y));
    points.push_back(Point_2(x+hside_x, y+hside_y));
    points.push_back(Point_2(x-hside_x, y+hside_y));

    points.push_back(Point_2(x-hside_x, y));
    points.push_back(Point_2(x, y-hside_y));
    points.push_back(Point_2(x+hside_x, y));
    points.push_back(Point_2(x, y+hside_y));

    // to introduce some asymmetry
    points.push_back(Point_2(x-hside_x/3., y+hside_y/2.));

    //some points really far away on each side
//    bbox_vertices.push_back(Point_2(2.*xmin-xmax, ymin+0.5*(ymax-ymin)));
//    bbox_vertices.push_back(Point_2(2.*xmax-xmin, ymin+0.5*(ymax-ymin)));
//    bbox_vertices.push_back(Point_2(xmin+0.5*(xmax-xmin), 2.*ymin-ymax));
//    bbox_vertices.push_back(Point_2(xmin+0.5*(xmax-xmin), 2.*ymax-ymin));

    //a net on the border just to have a clean border with Tangential_Complex
//    int n = 10;
//    for(int i=1; i<n; ++i)
//    {
//      FT xstep = (xmax-xmin)/(double) n;
//      FT ystep = (ymax-ymin)/(double) n;

//      bbox_vertices.push_back(Point_2(xmin, ymin+i*ystep));
//      bbox_vertices.push_back(Point_2(xmax, ymin+i*ystep));
//      bbox_vertices.push_back(Point_2(xmin+i*xstep, ymin));
//      bbox_vertices.push_back(Point_2(xmin+i*xstep, ymax));
//    }

    return points;
  }

  Point_container get_boundary_points(unsigned int, double) const
  {
    return initial_points();
  }

  Rectangle_domain(const FT hside_x_,
                   const FT hside_y_,
                   const Point_2 center_ = CGAL::ORIGIN)
    : hside_x(hside_x_), hside_y(hside_y_), center(center_)
  {
    hside = (std::max)(hside_x, hside_y);
  }

  ~Rectangle_domain() { }
};

} // Anisotropic_mesh_2
} // CGAL

#endif
