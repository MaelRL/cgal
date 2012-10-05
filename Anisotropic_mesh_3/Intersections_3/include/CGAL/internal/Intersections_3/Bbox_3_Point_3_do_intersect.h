#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_POINT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_POINT_3_DO_INTERSECT_H


#include <CGAL/Bbox_3.h>
#include <CGAL/Point_3.h>


namespace CGAL {
  
namespace internal {

  template <class K>
  bool do_intersect(const CGAL::Bbox_3& bbox,
                    const CGAL::Point_3<K>& p,
                    const K& k)
  {
    return bbox.xmin() <= p.x() &&   p.x() <= bbox.xmax()
      && bbox.ymin() <= p.y() &&   p.y() <= bbox.ymax()
      && bbox.zmin() <= p.z() &&   p.z() <= bbox.zmax();
  }
  template <class K>
  bool do_intersect(const CGAL::Point_3<K>& p,
                    const CGAL::Bbox_3& bbox,
                    const K& k)
  {
    return do_intersect(bbox, p, k);
  }

} //namespace internal

template<class K>
bool do_intersect(const CGAL::Bbox_3& bbox,
                  const CGAL::Point_3<K>& p)
{
  CGAL_PROFILER("[intersections Bbox-Pt]");
  return typename K::Do_intersect_3()(bbox,p);
}

template<class K>
bool do_intersect(const CGAL::Point_3<K>& p,
                  const CGAL::Bbox_3& bbox)
{
  CGAL_PROFILER("[intersections Bbox-Pt]");
  return typename K::Do_intersect_3()(bbox,p);
}

} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_POINT_3_DO_INTERSECT_H
