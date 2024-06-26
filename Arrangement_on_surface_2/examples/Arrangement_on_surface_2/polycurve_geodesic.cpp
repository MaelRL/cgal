//! \file examples/Arrangement_on_surface_2/polycurve_geodesic.cpp
// Constructing an arrangement of polygeodesics.

#include <vector>
#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_on_surface_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Segment_traits_2 = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Poly_traits_2 = CGAL::Arr_polyline_traits_2<Segment_traits_2>;

using Point_2 = Poly_traits_2::Point_2;
using Poly_curve_2 = Poly_traits_2::Curve_2;
using X_poly_curve_2 = Poly_traits_2::X_monotone_curve_2;
using Topol_poly_traits_2 =
  CGAL::Arr_spherical_topology_traits_2<Poly_traits_2>;
using Poly_arr =
  CGAL::Arrangement_on_surface_2<Poly_traits_2, Topol_poly_traits_2>;

using Seg_curve_2 = Segment_traits_2::Curve_2;
using X_seg_curve_2 = Segment_traits_2::X_monotone_curve_2;
using Topol_segment_traits_2 =
  CGAL::Arr_spherical_topology_traits_2<Segment_traits_2>;
using Segment_arr =
  CGAL::Arrangement_on_surface_2<Segment_traits_2, Topol_segment_traits_2>;

int main() {
  Segment_traits_2 seg_traits;
  auto ctr_p = seg_traits.construct_point_2_object();
  auto ctr_seg = seg_traits.construct_x_monotone_curve_2_object();

  Point_2 p1 = ctr_p(0, 1, -1);
  Point_2 p2 = ctr_p(-11, 7, -7);
  Point_2 p3 = ctr_p(-1, 0, 0);
  Point_2 p4 = ctr_p(-11, 7, 7);
  Point_2 p5 = ctr_p(-1, 1, 1);

  Segment_arr seg_arr(&seg_traits);
  X_seg_curve_2 seg_cv1 = ctr_seg(p1, p2);
  X_seg_curve_2 seg_cv2 = ctr_seg(p2, p3);
  X_seg_curve_2 seg_cv3 = ctr_seg(p3, p4);
  X_seg_curve_2 seg_cv4 = ctr_seg(p4, p5);

  insert(seg_arr, seg_cv1);
  insert(seg_arr, seg_cv2);
  insert(seg_arr, seg_cv3);
  insert(seg_arr, seg_cv4);
  std::cout << "# seg. vertives: " << seg_arr.number_of_vertices() << std::endl;

  std::list<Point_2> points;
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  points.push_back(p4);
  points.push_back(p5);

  Poly_traits_2 poly_traits;
  auto ctr = poly_traits.construct_x_monotone_curve_2_object();
  Poly_arr poly_arr(&poly_traits);
  insert(poly_arr, ctr(seg_cv1));
  insert(poly_arr, ctr(seg_cv2));
  insert(poly_arr, ctr(seg_cv3));
  insert(poly_arr, ctr(seg_cv4));
  std::cout << "# poly vertives: " << poly_arr.number_of_vertices() << std::endl;

  return 0;
}
