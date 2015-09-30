#ifndef CGAL_TETRAHEDRON_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_TETRAHEDRON_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/Triangle_3_Segment_3_do_intersect.h>

namespace CGAL
{

template <class K>
class Tetrahedron_3;

namespace internal
{

template <class K>
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& t,
             const typename K::Segment_3& s,
             const K & k)
{
  typedef typename K::Triangle_3        Triangle;

  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (s) );
  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (t) );

  if(do_intersect(s, Triangle(t[0], t[1], t[2]), k)) return true;
  if(do_intersect(s, Triangle(t[0], t[1], t[3]), k)) return true;
  if(do_intersect(s, Triangle(t[0], t[2], t[3]), k)) return true;
  if(do_intersect(s, Triangle(t[1], t[2], t[3]), k)) return true;

  CGAL_kernel_assertion(k.bounded_side_3_object()(t, s[0]) ==
                        k.bounded_side_3_object()(t, s[1]));

  return k.has_on_bounded_side_3_object()(t, s[0]);
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Segment_3& s,
             const typename K::Tetrahedron_3& t,
             const K & k)
{
  return do_intersect<K>(t, s, k);
}

} // namespace internal

CGAL_DO_INTERSECT_FUNCTION(Tetrahedron_3, Segment_3, 3)

} // namespace CGAL

#endif // CGAL_TETRAHEDRON_3_SEGMENT_3_DO_INTERSECT_H
