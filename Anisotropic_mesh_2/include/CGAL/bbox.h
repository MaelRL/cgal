#ifndef CGAL_ANISO_2_BBOX_2_H
#define CGAL_ANISO_2_BBOX_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Iso_rectangle_2.h>

namespace CGAL
{

template<typename K>
class Bbox : public CGAL::Bbox_2
{
public:
  typedef typename K::FT                                  FT;
  typedef typename K::Point_2                             Point_2;
  typedef typename K::Iso_rectangle_2                     Iso_rectangle;
  typedef typename K::Aff_transformation_2                Aff_transformation;

public:
  Bbox() {}

  Bbox(const typename CGAL::Bbox_2& bb)
    : CGAL::Bbox_2(bb) { }

  Bbox(double x_min, double y_min, double x_max, double y_max)
    : CGAL::Bbox_2(x_min, y_min, x_max, y_max) { }

public:
  CGAL::Bbox_2 bbox() { return *this; }

  CGAL::Bbox_2 transform(const Aff_transformation& t) const
  {
    FT x_min = DBL_MAX;   FT x_max = -DBL_MAX;
    FT y_min = DBL_MAX;   FT y_max = -DBL_MAX;

    Iso_rectangle rectangle(Point_2(xmin(), ymin()),
                            Point_2(xmax(), ymax()));
    for(unsigned int i = 0; i < 4; i++)
    {
      Point_2 pi = t.transform(rectangle[i]);
      x_min = (std::min)(x_min, pi.x());    x_max = (std::max)(x_max, pi.x());
      y_min = (std::min)(y_min, pi.y());    y_max = (std::max)(y_max, pi.y());
    }

    Iso_rectangle t_rectangle(Point_2(x_min, y_min),
                              Point_2(x_max, y_max));
    return t_rectangle.bbox();
  }

};  //end class Bounding_box

} // end namespace CGAL

#endif // CGAL_ANISO_2_BBOX_2_H
