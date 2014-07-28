#ifndef NEAREST_POINT_TETRAHEDRON_3_H_
#define NEAREST_POINT_TETRAHEDRON_3_H_

#include <CGAL/internal/AABB_tree/nearest_point_triangle_3.h>

#include <CGAL/kernel_basic.h>
#include <CGAL/enum.h>

namespace CGAL
{
namespace internal
{
template <class K>
inline
typename K::Point_3
nearest_point_3(const typename K::Point_3& origin,
                const typename K::Point_3& p1,
                const typename K::Point_3& p2,
                const typename K::Point_3& p3,
                const typename K::Point_3& p4,
                const K& k)
{
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_3 sq_distance =
      k.compute_squared_distance_3_object();

  const FT dist_origin_p1 = sq_distance(origin, p1);
  const FT dist_origin_p2 = sq_distance(origin, p2);
  const FT dist_origin_p3 = sq_distance(origin, p3);
  const FT dist_origin_p4 = sq_distance(origin, p4);

  if( dist_origin_p2 >= dist_origin_p1 &&
      dist_origin_p3 >= dist_origin_p1 &&
      dist_origin_p4 >= dist_origin_p1)
  {
    return p1;
  }
  if( dist_origin_p3 >= dist_origin_p2 &&
      dist_origin_p4 >= dist_origin_p2)
  {
    return p2;
  }
  if( dist_origin_p4 >= dist_origin_p3 )
  {
    return p3;
  }

  return p4;
}

template <class K>
inline
bool
is_inside_tetrahedron_3_aux(const typename K::Point_3& p0,
                            const typename K::Point_3& p1,
                            const typename K::Point_3& p2,
                            const typename K::Point_3& origin,
                            const typename K::Point_3& gp,
                            typename K::Point_3& result,
                            typename K::FT& sq_d_orig_result,
                            const typename K::Point_3& bound,
                            const K& k)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Plane_3 Plane_3;

  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();

  Triangle_3 tr(p0, p1, p2);
  Plane_3 plane = tr.supporting_plane();

  if(plane.oriented_side(gp) != plane.oriented_side(origin))
  {
    //point is outside of tet && the triangle is visible
    //so we get the shortest dist to triangle and compare
    Point_3 closest_point = CGAL::internal::nearest_point_3(origin, tr, bound, k);
    FT sq_d_orig_clpoint = sq_distance(origin, closest_point);

    if(sq_d_orig_clpoint < sq_d_orig_result)
    {
      result = closest_point;
      sq_d_orig_result = sq_d_orig_clpoint;
    }
    return false;
  }
  return true;
}

template <class K>
inline
bool
is_inside_tetrahedron_3(const typename K::Point_3& origin,
                        const typename K::Tetrahedron_3& tet,
                        typename K::Point_3& result,
                        const typename K::Point_3& bound,
                        const K& k)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on =
      k.construct_vertex_3_object();

  const Point_3& tet0 = vertex_on(tet,0);
  const Point_3& tet1 = vertex_on(tet,1);
  const Point_3& tet2 = vertex_on(tet,2);
  const Point_3& tet3 = vertex_on(tet,3);

  Point_3 gp = CGAL::barycenter(tet0, 0.25, tet1, 0.25, tet2, 0.25, tet3, 0.25);
  FT sq_d_orig_result = 1e30;

  bool b0 = is_inside_tetrahedron_3_aux(tet1, tet3, tet2, origin, gp, result, sq_d_orig_result, bound, k);
  bool b1 = is_inside_tetrahedron_3_aux(tet0, tet2, tet3, origin, gp, result, sq_d_orig_result, bound, k);
  bool b2 = is_inside_tetrahedron_3_aux(tet0, tet3, tet1, origin, gp, result, sq_d_orig_result, bound, k);
  bool b3 = is_inside_tetrahedron_3_aux(tet0, tet1, tet2, origin, gp, result, sq_d_orig_result, bound, k);

  return (b0 && b1 && b2 && b3);
}

template <class K>
typename K::Point_3
nearest_point_3(const typename K::Point_3& origin,
                const typename K::Tetrahedron_3& tetrahedron,
                const typename K::Point_3& bound,
                const K& k)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;

  typename K::Compute_squared_distance_3 sq_distance =
      k.compute_squared_distance_3_object();
  typename K::Compare_squared_distance_3 compare_sq_distance =
      k.compare_squared_distance_3_object();

  // Distance from origin to bound
  const FT bound_sq_dist = sq_distance(origin, bound);

  Point_3 closest_point;
  bool inside = is_inside_tetrahedron_3(origin, tetrahedron, closest_point, bound, k);

  if(inside)
    return origin;
  else if(compare_sq_distance(origin, closest_point, bound_sq_dist) == CGAL::LARGER)
    return bound;
  else
    return closest_point;
}

} // end namespace internal


template <class K>
inline
Point_3<K>
nearest_point_3(const Point_3<K>& origin,
                const Tetrahedron_3<K>& tetrahedron,
                const Point_3<K>& bound)
{
  return internal::nearest_point_3(origin, tetrahedron, bound, K());
}

} // end namespace CGAL


#endif // NEAREST_POINT_TETRAHEDRON_3_H_
