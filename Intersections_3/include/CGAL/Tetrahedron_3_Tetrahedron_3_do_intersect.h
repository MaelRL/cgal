#ifndef CGAL_TETRAHEDRON_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_TETRAHEDRON_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Triangle_3_Triangle_3_do_intersect.h>

namespace CGAL
{

template <class K>
class Tetrahedron_3;

namespace internal
{

// This code is not optimized:
template <class K>
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet_1,
             const typename K::Tetrahedron_3& tet_2,
             const K & k)
{
  typedef typename K::Triangle_3 Triangle;

  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (tet_1) );
  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (tet_2) );

  for(std::size_t i=0; i<4; ++i)
  {
    Triangle tr_1(tet_1[i], tet_1[(i+1)%4], tet_1[(i+2)%4]);
    for(std::size_t j=0; j<4; ++j)
    {
      Triangle tr_2(tet_2[j], tet_2[(j+1)%4], tet_2[(j+2)%4]);

      if(do_intersect(tr_1, tr_2, k))
        return true;
    }
  }

  CGAL_kernel_assertion(k.bounded_side_3_object()(tet_1, tet_2[0]) ==
                        k.bounded_side_3_object()(tet_1, tet_2[1]));
  CGAL_kernel_assertion(k.bounded_side_3_object()(tet_1, tet_2[0]) ==
                        k.bounded_side_3_object()(tet_1, tet_2[2]));
  CGAL_kernel_assertion(k.bounded_side_3_object()(tet_1, tet_2[0]) ==
                        k.bounded_side_3_object()(tet_1, tet_2[3]));

  CGAL_kernel_assertion(k.bounded_side_3_object()(tet_2, tet_1[0]) ==
                        k.bounded_side_3_object()(tet_2, tet_1[1]));
  CGAL_kernel_assertion(k.bounded_side_3_object()(tet_2, tet_1[0]) ==
                        k.bounded_side_3_object()(tet_2, tet_1[2]));
  CGAL_kernel_assertion(k.bounded_side_3_object()(tet_2, tet_1[0]) ==
                        k.bounded_side_3_object()(tet_2, tet_1[3]));

  return (k.has_on_bounded_side_3_object()(tet_1, tet_2[0]) ||
          k.has_on_bounded_side_3_object()(tet_2, tet_1[0]));
}

} // namespace internal

CGAL_DO_INTERSECT_FUNCTION_SELF(Tetrahedron_3, 3)

} // namespace CGAL

#endif // CGAL_TETRAHEDRON_3_TETRAHEDRON_3_DO_INTERSECT_H
