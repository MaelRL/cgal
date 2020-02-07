// Copyright(c) 2018, 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL(www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAPPING_SNAP_FACES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAPPING_SNAP_FACES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/circulator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/concurrent_vector.h>
# include <tbb/parallel_for.h>
#endif

#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <type_traits>
#include <utility>
#include <unordered_set>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename TriangleMesh, typename VPM>
struct Close_face_detector
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor                    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor                  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor                      face_descriptor;

  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 2, face_descriptor>            Box;

  typedef typename boost::property_traits<VPM>::value_type                                 Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                                      Kernel;

  typedef CGAL::Exact_predicates_exact_constructions_kernel                                Exact_kernel;
  typedef Exact_kernel::FT                                                                 EFT;
  typedef Exact_kernel::Point_2                                                            EPoint_2;
  typedef Exact_kernel::Segment_2                                                          ESegment_2;
  typedef Exact_kernel::Triangle_2                                                         ETriangle_2;
  typedef Exact_kernel::Point_3                                                            EPoint_3;
  typedef Exact_kernel::Vector_3                                                           EVector_3;
  typedef Exact_kernel::Line_3                                                             ELine_3;
  typedef Exact_kernel::Plane_3                                                            EPlane_3;

  typedef Exact_kernel::Intersect_3                                                        EIntersect_3;

  typedef CGAL::dynamic_vertex_property_t<EPoint_2>                                        Vertex_point_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_point_tag>::const_type         PPM;

  typedef CGAL::Cartesian_converter<Kernel, Exact_kernel>                                  To_exact;
  typedef CGAL::Cartesian_converter<Exact_kernel, Kernel>                                  To_input;

public:
  CGAL::Bbox_2 bbox_of_projected_face(const face_descriptor f,
                                      const TriangleMesh& tm,
                                      const PPM ppm)
  {
    const halfedge_descriptor h = halfedge(f, tm);

    return get(ppm, target(h, tm)).bbox() +
           get(ppm, target(next(h, tm), tm)).bbox() +
           get(ppm, target(prev(h, tm), tm)).bbox();
  }

  struct Point_collector
    : public boost::static_visitor<>
  {
    Point_collector(std::vector<EPoint_2>& ipoints) : m_ipoints(ipoints) { }

    void operator()(const EPoint_2& p)
    {
      m_ipoints.push_back(p);
    }

    void operator()(const ESegment_2& s)
    {
      m_ipoints.push_back(s.source());
      m_ipoints.push_back(s.target());
    }

    void operator()(const ETriangle_2& t)
    {
      m_ipoints.push_back(t[0]);
      m_ipoints.push_back(t[1]);
      m_ipoints.push_back(t[2]);
    }

    void operator()(const std::vector<EPoint_2>& pts)
    {
      m_ipoints.insert(std::end(m_ipoints), std::begin(pts), std::end(pts));
    }

  private:
    std::vector<EPoint_2>& m_ipoints;
  };

  struct Face_inter_info
  {
    Face_inter_info(const face_descriptor f)
      : face(f), closest_face(NULL), min_distance(std::numeric_limits<double>::max())
    { }

    void add_contact(face_descriptor other_face,
                     const double min_d,
                     EFT /* max_d */)
    {
      if(min_d < min_distance)
      {
        min_distance = min_d;
        closest_face = other_face;
      }
    }

  private:
    const face_descriptor face; // the face considered
    face_descriptor closest_face; // face that realizes the min line-distance
    double min_distance; // min line-distance to another face
  };

  struct Box_box_intersection
  {
    Box_box_intersection(const EPlane_3& plane,
                         const EVector_3& direction,
                         const TriangleMesh& tm,
                         const PPM ppm,
                         const VPM vpm,
                         std::vector<Face_inter_info>& result)
      :
        m_result(result),
        m_projection_plane(plane),
        m_direction(direction),
        m_tm(tm),
        m_vpm(vpm),
        m_ppm(ppm)
    { }

private:
    ETriangle_2 projected_triangle(const face_descriptor f,
                                   const PPM ppm)
    {
      const halfedge_descriptor h = halfedge(f, m_tm);

      return ETriangle_2(get(ppm, target(h, m_tm)),
                         get(ppm, target(next(h, m_tm), m_tm)),
                         get(ppm, target(prev(h, m_tm), m_tm)));
    }

    EPlane_3 plane(const face_descriptor f)
    {
      const halfedge_descriptor h = halfedge(f, m_tm);
      To_exact to_exact;

      return EPlane_3(to_exact(get(m_vpm, target(h, m_tm))),
                      to_exact(get(m_vpm, target(next(h, m_tm), m_tm))),
                      to_exact(get(m_vpm, target(prev(h, m_tm), m_tm))));
    }

    double signed_line_distance_of_faces(const face_descriptor f1,
                                         const face_descriptor f2,
                                         const EPoint_3& pt,
                                         const EVector_3& direction)
    {
      const ELine_3 line(pt, direction);
      const EPlane_3 plane1 = plane(f1);
      const EPlane_3 plane2 = plane(f2);

      CGAL::cpp11::result_of<EIntersect_3(EPlane_3, ELine_3)>::type
          inter1 = CGAL::intersection(plane1, line),
          inter2 = CGAL::intersection(plane2, line);

      const EPoint_3* pt1_ptr = boost::get<EPoint_3>(&*inter1);
      const EPoint_3* pt2_ptr = boost::get<EPoint_3>(&*inter2);

      CGAL_assertion(pt1_ptr != NULL && pt2_ptr != NULL);

      CGAL::Interval_nt_advanced::Protector protector;

      return (*pt2_ptr - *pt1_ptr) * direction > 0 ?
            CGAL::sqrt(CGAL::approx(CGAL::squared_distance(*pt1_ptr, *pt2_ptr))).inf():
            -CGAL::sqrt(CGAL::approx(CGAL::squared_distance(*pt1_ptr, *pt2_ptr))).sup();
    }

public:
    template <class Box>
    void operator()(const Box* b1, const Box* b2)
    {
      const ETriangle_2 t1 = projected_triangle(b1->info(), m_ppm);
      const ETriangle_2 t2 = projected_triangle(b2->info(), m_ppm);

      if(CGAL::do_intersect(t1, t2))
      {
        std::vector<EPoint_2> ipoints;
        Point_collector collector(ipoints);
        boost::apply_visitor(collector, *CGAL::intersection(t1, t2));

        std::size_t nb_ipts = ipoints.size();

        double min_d = signed_line_distance_of_faces(b1->info(), b2->info(),
                                                     m_projection_plane.to_3d(ipoints[0]),
                                                     m_direction);
        double max_d = min_d;

        for(std::size_t i=1; i<nb_ipts; ++i)
        {
          double distance = signed_line_distance_of_faces(b1->info(), b2->info(),
                                                          m_projection_plane.to_3d(ipoints[i]),
                                                          m_direction);
          if(distance < min_d)
            min_d = distance;

          if(distance > max_d)
            max_d = distance;
        }

        m_result[b1->info()].add_contact(b2->info(), min_d, max_d);
        m_result[b2->info()].add_contact(b1->info(), min_d, max_d);
      }
    }

  private:
    std::vector<Face_inter_info>& m_result;

    const EPlane_3& m_projection_plane;
    const EVector_3& m_direction;
    const TriangleMesh& m_tm;
    const VPM m_vpm;
    const PPM m_ppm;
  };

  void operator()(const face_descriptor f,
                  const double sq_d,
                  const TriangleMesh& tm,
                  const VPM vpm)
  {
    // Contains the precomputed projections of 3D points onto a plane
    PPM ppm = get(Vertex_point_tag(), tm);

    To_exact to_exact;

    const EVector_3 direction = to_exact(compute_face_normal(f, tm, parameters::vertex_point_map(vpm)));
    EPlane_3 projection_plane(EPoint_3(0,0,0), direction);

    for(vertex_descriptor v : vertices(tm))
      put(ppm, v, projection_plane.to_2d(to_exact(get(vpm, v))));

    Box f_box(bbox_of_projected_face(f, tm, ppm), f);

    // Collect faces from the fixed polyhedron and set a unique id to each face
    std::vector<Box> boxes;
    for(face_descriptor tf : faces(tm))
      boxes.push_back(Box(bbox_of_projected_face(tf, tm, ppm), tf));

    std::vector<const Box*> f_box_ptr;
    f_box_ptr.push_back(&f_box);

    std::vector<const Box*> boxes_ptrs;
    for(const Box& b : boxes)
      boxes_ptrs.push_back(&b);

    std::vector<Face_inter_info> result;
    CGAL::box_intersection_d(f_box_ptr.begin(), f_box_ptr.end(),
                             boxes_ptrs.begin(), boxes_ptrs.end(),
                             Box_box_intersection(projection_plane, direction, tm, ppm, vpm, result));
  }
};

// Obviously terrible complexity for now
template <typename TriangleMesh, typename VPM>
void detect_face_proximity(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                           const double sq_d,
                           const TriangleMesh& tm,
                           const VPM vpm)
{
  Close_face_detector<TriangleMesh, VPM> proximity_detector;
  proximity_detector(f, sq_d, tm, vpm);
}

template <typename TriangleMesh, typename ToleranceMap, typename VPM>
void mark_close_faces(TriangleMesh& tm,
                      const ToleranceMap tolerance_map,
                      const VPM vpm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  for(const face_descriptor f : faces(tm))
  {
    const halfedge_descriptor h = halfedge(f, tm);

    double tol = get(tolerance_map, target(next(h, tm), tm));
    tol = (std::min)(tol, get(tolerance_map, source(h, tm)));
    tol = (std::min)(tol, get(tolerance_map, target(h, tm)));

    detect_face_proximity(f, CGAL::square(tol), tm, vpm);
  }
}

template <typename FaceRange, typename TriangleMesh>
void project_borders(const FaceRange& faces,
                     TriangleMesh& tm)
{

}

template <typename TriangleMesh, typename ToleranceMap, typename VPM>
void mark_faces_to_be_removed(TriangleMesh& tm,
                              const ToleranceMap tolerance_map,
                              const VPM vpm)
{
  bool stable = false;
  while(!stable)
  {
    mark_close_faces(tm, tolerance_map, vpm);
    project_borders();
  }
}

template <typename TriangleMesh>
void remove_marked_faces(TriangleMesh& tm)
{

}

template <typename TriangleMesh>
void close_gap_with_snap(TriangleMesh& tm)
{

}

template <typename TriangleMesh>
void close_gap_with_hole_filling(TriangleMesh& tm)
{

}

template <typename TriangleMesh>
void close_gap(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
               TriangleMesh& tm)
{
  close_gap_with_snap(tm);
  close_gap_with_hole_filling(tm);
}

template <typename TriangleMesh>
void close_gaps(TriangleMesh& tm)
{
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h;
  for(;;) // holes that matter, not all of them
    close_gap(h, tm);
}

} // namespace internal

namespace experimental {

template <typename TriangleMesh, typename ToleranceMap, typename NamedParameters>
std::size_t snap_faces(TriangleMesh& tm,
                       const ToleranceMap tolerance_map,
                       const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type       VPM;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(boost::vertex_point, tm));

  internal::mark_faces_to_be_removed(tm, tolerance_map, vpm);

  internal::remove_marked_faces(tm);

  internal::close_gaps(tm);
}

template <typename TriangleMesh, typename ToleranceMap>
std::size_t snap_faces(TriangleMesh& tm,
                       const ToleranceMap tolerance_map)
{
  return snap_faces(tm, tolerance_map, CGAL::parameters::all_default());
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAPPING_SNAP_FACES_H
