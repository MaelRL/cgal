#ifndef CGAL_ANISO_2_BBOX_3_H
#define CGAL_ANISO_2_BBOX_3_H

#include <CGAL/Bbox_3.h>

namespace CGAL
{
  template<typename K>
  class Aniso_bbox_3 : public CGAL::Bbox_3
  {
  public:
    typedef typename K::FT                                  FT;
    typedef typename K::Point_3                             Point_3;
    typedef typename K::Iso_cuboid_3                        Iso_cuboid;
    typedef typename K::Aff_transformation_3                Aff_transformation;

  public:
    Aniso_bbox_3() {}

    Aniso_bbox_3(const typename CGAL::Bbox_3& bb)
      : CGAL::Bbox_3(bb) { }

    Aniso_bbox_3(double x_min, double y_min, double z_min,
                 double x_max, double y_max, double z_max)
      : CGAL::Bbox_3(x_min, y_min, z_min, x_max, y_max, z_max) { }

  public:
    CGAL::Bbox_3 bbox()
    {
      return *this;
    }

    CGAL::Bbox_3 transform(const Aff_transformation& t) const
    {
      FT x_min = DBL_MAX;   FT x_max = -DBL_MAX;
      FT y_min = DBL_MAX;   FT y_max = -DBL_MAX;
      FT z_min = DBL_MAX;   FT z_max = -DBL_MAX;

      Iso_cuboid cuboid(Point_3(xmin(), ymin(), zmin()),
                        Point_3(xmax(), ymax(), zmax()));
      for(unsigned int i = 0; i < 8; i++)
      {
        Point_3 pi = t.transform(cuboid[i]);
        x_min = (std::min)(x_min, pi.x());    x_max = (std::max)(x_max, pi.x());
        y_min = (std::min)(y_min, pi.y());    y_max = (std::max)(y_max, pi.y());
        z_min = (std::min)(z_min, pi.z());    z_max = (std::max)(z_max, pi.z());
      }

      Iso_cuboid t_cuboid(Point_3(x_min, y_min, z_min),
                          Point_3(x_max, y_max, z_max));
      return t_cuboid.bbox();
    }

  };  //end class Bounding_box
} // end namespace CGAL

#endif // CGAL_ANISO_2_BBOX_3_H
