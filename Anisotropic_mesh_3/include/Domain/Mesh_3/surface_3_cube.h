#ifndef CGAL_ANISOTROPIC_MESH_3_SURFACE_3_CUBE_H
#define CGAL_ANISOTROPIC_MESH_3_SURFACE_3_CUBE_H

#include <list>
#include <vector>

namespace CGAL
{

template<typename K>
typename K::FT cube_function(const typename K::Point_3& p)
{
  if( p.x() > 0. && p.x() < 3. &&
      p.y() > 0. && p.y() < 3. &&
      p.z() > 0. && p.z() < 3. )
    return -1.;
  return 1.;
}

template<typename K>
void create_polylines(std::list<std::vector<typename K::Point_3> >& polylines)
{
  typedef typename K::FT                      FT;
  typedef typename K::Point_3                 Point_3;
  typedef std::vector<Point_3>                Polyline_3;

  FT side_x = 3.; // fixme don't hardcode stuff like that...
  FT side_y = 3.;
  FT side_z = 3.;

  // bot 4 edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, 0, 0));
    polyline.push_back(Point_3(side_x, 0, 0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, 0, 0));
    polyline.push_back(Point_3(side_x, side_y, 0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, side_y, 0));
    polyline.push_back(Point_3(0, side_y, 0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, side_y, 0));
    polyline.push_back(Point_3(0, 0, 0));
    polylines.push_back(polyline);
  }

  // vertical edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, 0, 0));
    polyline.push_back(Point_3(0, 0, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, 0, 0));
    polyline.push_back(Point_3(side_x, 0, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, side_y, 0));
    polyline.push_back(Point_3(0, side_y, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, side_y, 0));
    polyline.push_back(Point_3(side_x, side_y, side_z));
    polylines.push_back(polyline);
  }

  // top 4 edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, 0, side_z));
    polyline.push_back(Point_3(side_x, 0, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, 0, side_z));
    polyline.push_back(Point_3(side_x, side_y, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,side_y,side_z));
    polyline.push_back(Point_3(0,side_y,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,side_y,side_z));
    polyline.push_back(Point_3(0,0,side_z));
    polylines.push_back(polyline);
  }
}

template<typename Domain>
Domain const * cube_domain()
{
  typedef typename Domain::R                     Traits;
  typedef typename Traits::Point_3               Point_3;
  typedef typename Traits::Sphere_3              Sphere_3;

  typename Traits::FT radius = 5.*5.;
  return new Domain(cube_function<Traits>,
                    Sphere_3(Point_3(1., 1., 1.), radius),
                    1e-6);
}

template<typename Domain>
Domain const * cube_domain_with_features()
{
  typedef typename Domain::R                     Traits;
  typedef typename Traits::Point_3               Point_3;
  typedef typename Traits::Sphere_3              Sphere_3;

  typedef std::vector<Point_3>                   Polyline_3;
  typedef std::list<Polyline_3>                  Polylines;

  typename Traits::FT radius = 5.*5.;
  Domain* dom = new Domain(cube_function<Traits>,
                           Sphere_3(Point_3(1., 1., 1.), radius),
                           1e-6);

  // Polylines
  Polylines polylines;
  create_polylines<Traits>(polylines);
  dom->add_features(polylines.begin(), polylines.end());

  return dom;
}

} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_SURFACE_3_CUBE_H
