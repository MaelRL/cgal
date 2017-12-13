// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Maxime Gimeno and Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Side_of_triangle_mesh.h>

namespace CGAL {
namespace internal {

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class VertexPointMap1,
         class VertexPointMap2>
struct Intersect_faces
{

  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;

  // members
  const TM& m_tmesh1;
  const VertexPointMap1 m_vpmap1;
  const TM& m_tmesh2;
  const VertexPointMap2 m_vpmap2;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;





  Intersect_faces(const TM& mesh1, const TM& mesh2,
                   OutputIterator it,
                   VertexPointMap1 vpmap1, VertexPointMap2 vpmap2,
                   const Kernel& kernel)
    :
      m_tmesh1(mesh1),
      m_vpmap1(vpmap1),
      m_tmesh2(mesh2),
      m_vpmap2(vpmap2),
      m_iterator(it),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tmesh1);
    halfedge_descriptor g = halfedge(c->info(), m_tmesh2);


    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap1, target(h,m_tmesh1)),
                                    get(m_vpmap1, target(next(h,m_tmesh1),m_tmesh1)),
                                    get(m_vpmap1, source(h,m_tmesh1)));

    Triangle t2 = triangle_functor( get(m_vpmap2, target(g,m_tmesh2)),
                                    get(m_vpmap2, target(next(g,m_tmesh2),m_tmesh2)),
                                    get(m_vpmap2, source(g,m_tmesh2)));
    if(do_intersect_3_functor(t1, t2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_faces

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class Polyline,
         class VertexPointMap>
struct Intersect_face_polyline
{
  // wrapper to check whether anything is inserted to output iterator

  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point;


  // members
  const TM& m_tmesh;
  const std::vector<face_descriptor>& faces;
  const VertexPointMap m_vpmap;
  const Polyline& polyline;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_face_polyline(const TM& mesh,
                           const std::vector<face_descriptor>& faces,
                           const Polyline& polyline,
                           OutputIterator it,
                           VertexPointMap vpmap,
                           const Kernel& kernel)
    :
      m_tmesh(mesh),
      faces(faces),
      m_vpmap(vpmap),
      polyline(polyline),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(faces[b->info()], m_tmesh);


    // check for geometric intersection
    Triangle t = triangle_functor( get(m_vpmap, target(h,m_tmesh)),
                                   get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                   get(m_vpmap, source(h,m_tmesh)));

    Segment s = segment_functor(polyline[c->info()], polyline[c->info() + 1]);
    if(do_intersect_3_functor(t, s)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_face_polyline

template<class TM,
         class Kernel,
         class Box,
         class PolylineRange,
         class OutputIterator,
         class VertexPointMap>
struct Intersect_face_polylines
{
  // wrapper to check whether anything is inserted to output iterator

  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point;


  // members
  const TM& m_tmesh;
  const std::vector<face_descriptor>& faces;
  const VertexPointMap m_vpmap;
  const PolylineRange& polylines;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_face_polylines(const TM& mesh,
                           const std::vector<face_descriptor>& faces,
                           const PolylineRange& polylines,
                           OutputIterator it,
                           VertexPointMap vpmap,
                           const Kernel& kernel)
    :
      m_tmesh(mesh),
      faces(faces),
      m_vpmap(vpmap),
      polylines(polylines),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(faces[b->info().second], m_tmesh);


    // check for geometric intersection
    Triangle t = triangle_functor( get(m_vpmap, target(h,m_tmesh)),
                                   get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                   get(m_vpmap, source(h,m_tmesh)));

    Segment s = segment_functor(polylines[c->info().first][c->info().second], polylines[c->info().first][c->info().second + 1]);
    if(do_intersect_3_functor(t, s)){
      *m_iterator++ = std::make_pair(b->info().second, c->info());
    }
  } // end operator ()
}; // end struct Intersect_face_polylines


template<class Polyline,
         class Kernel,
         class Box,
         class OutputIterator>
struct Intersect_polylines
{


  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Point_3 Point;


  // members
  const Polyline& polyline1;
  const Polyline& polyline2;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_polylines(const Polyline& polyline1,
                      const Polyline& polyline2,
                      OutputIterator it,
                      const Kernel& kernel)
    :
      polyline1(polyline1),
      polyline2(polyline2),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {


    // check for geometric intersection

    Segment s1 = segment_functor(polyline1[b->info()], polyline1[b->info() + 1]);
    Segment s2 = segment_functor(polyline2[c->info()], polyline2[c->info() + 1]);
    if(do_intersect_3_functor(s1, s2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_polylines

template<class PolylineRange,
         class Kernel,
         class Box,
         class OutputIterator>
struct Intersect_polyline_ranges
{


  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Point_3 Point;


  // members
  const PolylineRange& polyline1;
  const PolylineRange& polyline2;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_polyline_ranges(const PolylineRange& polyline1,
                            const PolylineRange& polyline2,
                            OutputIterator it,
                            const Kernel& kernel)
    :
      polyline1(polyline1),
      polyline2(polyline2),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {


    // check for geometric intersection

    Segment s1 = segment_functor(polyline1[b->info().first][b->info().second], polyline1[b->info().first][b->info().second + 1]);
    Segment s2 = segment_functor(polyline2[c->info().first][c->info().second], polyline2[c->info().first][c->info().second + 1]);

    if(do_intersect_3_functor(s1, s2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_polyline_ranges

struct Throw_at_first_output {
  class Throw_at_first_output_exception: public std::exception
  { };

  template<class T>
  void operator()(const T& /* t */) const {
    throw Throw_at_first_output_exception();
  }
};

// Note this is not officially documented
/*
 * reports all the pairs of faces intersecting between two triangulated surface meshes.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(tm1)`
 * \pre `CGAL::is_triangle_mesh(tm2)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `RandomAccessRange`.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param face_range1 the range of faces of `tm1` to check for intersections.
 * \param face_range2 the range of faces of `tm2` to check for intersections.
 * \param tm1 the first triangulated surface mesh.
 * \param tm2 the second triangulated surface mesh.
 * \param out output iterator to be filled with all pairs of faces that intersect.
 *  First and second element in the pairs correspond to faces of `tm1` and `tm2` respectively
 * \param np1 optional sequence of \ref namedparameters for `tm1`, among the ones listed below
 * \param np2 optional sequence of \ref namedparameters for `tm2`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh,
          class FaceRange,
          class OutputIterator,
          class NamedParameters1,
          class NamedParameters2>
OutputIterator
compute_face_face_intersection(const FaceRange& face_range1,
                               const FaceRange& face_range2,
                               const TriangleMesh& tm1,
                               const TriangleMesh& tm2,
                               OutputIterator out,
                               const NamedParameters1& np1,
                               const NamedParameters2& np2)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm1));
  CGAL_precondition(CGAL::is_triangle_mesh(tm2));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor> Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range1), boost::end(face_range1) )
        );
  boxes2.reserve(
        std::distance( boost::begin(face_range2), boost::end(face_range2) )
        );

  typedef typename GetVertexPointMap<TM, NamedParameters1>::const_type VertexPointMap1;
  typedef typename GetVertexPointMap<TM, NamedParameters2>::const_type VertexPointMap2;

  VertexPointMap1 vpmap1 = boost::choose_param(get_param(np1, internal_np::vertex_point),
                                              get_const_property_map(boost::vertex_point, tm1));
  VertexPointMap2 vpmap2 = boost::choose_param(get_param(np2, internal_np::vertex_point),
                                              get_const_property_map(boost::vertex_point, tm2));
  CGAL_static_assertion(
      (boost::is_same<
       typename boost::property_traits<VertexPointMap1>::value_type,
       typename boost::property_traits<VertexPointMap2>::value_type
       >::value) );
  BOOST_FOREACH(face_descriptor f, face_range1)
  {
    boxes1.push_back(Box(Polygon_mesh_processing::face_bbox(f, tm1), f));
  }

  BOOST_FOREACH(face_descriptor f, face_range2)
  {
    boxes2.push_back(Box(Polygon_mesh_processing::face_bbox(f, tm2), f));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr(boost::make_counting_iterator<const Box*>(&boxes1[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes1[0]+boxes1.size()));
  std::vector<const Box*> box2_ptr(boost::make_counting_iterator<const Box*>(&boxes2[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes2[0]+boxes2.size()));


  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters1>::type GeomTraits;
  GeomTraits gt = boost::choose_param(get_param(np1, internal_np::geom_traits), GeomTraits());

  CGAL::internal::Intersect_faces<TM,
      GeomTraits,
      Box,
      OutputIterator,
      VertexPointMap1,
      VertexPointMap2>
      Intersect_faces(tm1, tm2,
                       out,
                       vpmap1, vpmap2,
                       gt);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           Intersect_faces,cutoff);
  return Intersect_faces.m_iterator;
}

// Note this is not officially documented
/*
 * reports all the pairs of segments and faces intersecting between
 * a triangulated surface mesh and a polyline.
 * \attention If a polyline vertex intersects a face, the intersection will
 * be reported twice (even more if it is on a vertex, edge, or point).
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(tm)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `RandomAccessRange`.
 * \tparam Polyline a `RandomAccessRange` of points. The point type of the range must be
 * the same as the value type of the vertex point map.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`.This `OutputIterator` will hold the position of the
 *  elements in their respective range. In the case of the polyline, this position is the index of
 * the segment intersecting the face (which is the index  of the first point of the
 * segment following the range order.)
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param face_range the range of faces of `tm` to check for intersections.
 * \param polyline the polyline to check for intersections.
 * \param tm the triangulated surface mesh to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np optional sequence of \ref namedparameters for `tm`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tm`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh,
          class FaceRange,
          class Polyline,
          class OutputIterator,
          class NamedParameters>
OutputIterator
compute_face_polyline_intersection( const FaceRange& face_range,
               const Polyline& polyline,
               const TriangleMesh& tm,
               OutputIterator out,
               const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;

  VertexPointMap vpmap = boost::choose_param(get_param(np, internal_np::vertex_point),
                                             get_const_property_map(boost::vertex_point, tm));
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;
  CGAL_static_assertion(
        (boost::is_same<Point,
        typename boost::range_value<Polyline>::type>::value));

  std::vector<face_descriptor> faces;
  faces.reserve(std::distance( boost::begin(face_range), boost::end(face_range) ));

  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t> Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range), boost::end(face_range) )
        );

  boxes2.reserve(
        std::distance( boost::begin(polyline), boost::end(polyline) ) - 1
        );


  BOOST_FOREACH(face_descriptor f, face_range)
  {
    faces.push_back(f);
    boxes1.push_back(Box(Polygon_mesh_processing::face_bbox(f, tm), faces.size()-1));
  }

  for(std::size_t i =0; i< polyline.size()-1; ++i)
  {
    Point p1 = polyline[i];
    Point p2 = polyline[i+1];
    boxes2.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  // generate box pointers

  std::vector<const Box*> box1_ptr(boost::make_counting_iterator<const Box*>(&boxes1[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes1[0]+boxes1.size()));
  std::vector<const Box*> box2_ptr(boost::make_counting_iterator<const Box*>(&boxes2[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes2[0]+boxes2.size()));

  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;
  GeomTraits gt = boost::choose_param(get_param(np, internal_np::geom_traits), GeomTraits());

  CGAL::internal::Intersect_face_polyline<TM,
      GeomTraits,
      Box,
      OutputIterator,
      Polyline,
      VertexPointMap>
      Intersect_face_polyline(tm,
                               faces,
                               polyline,
                               out,
                               vpmap,
                               gt);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           Intersect_face_polyline, cutoff);
  return Intersect_face_polyline.m_iterator;
}

// Note this is not officially documented
/*
 * reports all the pairs of segments and faces intersecting between
 * a triangulated surface mesh and a range of polylines.
 * \attention If a polyline vertex intersects a face, the intersection will
 * be reported twice (even more if it is on a vertex, edge, or point).
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `RandomAccessRange`.
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points. The point type of the range must be
 * the same as the value type of the vertex point map.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::pair<std::size_t, std::size_t> >`.
 * Each pair holds the index of the face and a pair containing the index of the polyline in the range and the index of
 * the first point of the segment in the polyline.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param face_range the range of `mesh` faces to check for intersections.
 * \param polyline_range the range of polylines to check for intersections.
 * \param mesh the triangulated surface mesh to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np optional sequence of \ref namedparameters for `mesh`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh,
          class FaceRange,
          class PolylineRange,
          class OutputIterator,
          class NamedParameters>
OutputIterator
compute_face_polylines_intersection(const FaceRange& face_range,
                                    const PolylineRange& polyline_range,
                                    const TriangleMesh& mesh,
                                    OutputIterator out,
                                    const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;

  VertexPointMap vpmap = boost::choose_param(get_param(np, internal_np::vertex_point),
                                             get_const_property_map(boost::vertex_point, mesh));
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;
  typedef typename boost::range_value<PolylineRange>::type Polyline;
  CGAL_static_assertion(
        (boost::is_same<Point,
        typename boost::range_value<Polyline>::type>::value));

  std::vector<face_descriptor> faces;
  faces.reserve(std::distance( boost::begin(face_range), boost::end(face_range) ));

  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::pair<std::size_t, std::size_t> > Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range), boost::end(face_range) )
        );

  int polylines_size = 0;
  BOOST_FOREACH(Polyline poly, polyline_range)
  {
    polylines_size += std::distance( boost::begin(poly), boost::end(poly) ) -1;
  }
  boxes2.reserve( polylines_size );

  BOOST_FOREACH(face_descriptor f, face_range)
  {
    faces.push_back(f);
    boxes1.push_back(Box(Polygon_mesh_processing::face_bbox(f, mesh), std::make_pair(0, faces.size()-1)));
  }
  int range_size = std::distance( boost::begin(polyline_range), boost::end(polyline_range) );
  for(int j = 0; j < range_size; ++j)
  {
    Polyline poly = polyline_range[j];
    int size = std::distance( boost::begin(poly), boost::end(poly) );
    for(int i =0; i< size - 1; ++i)
    {
      Point p1 = poly[i];
      Point p2 = poly[i+1];
      boxes2.push_back(Box(p1.bbox() + p2.bbox(), std::make_pair(j, i)));
    }
  }

  // generate box pointers

  std::vector<const Box*> box1_ptr(boost::make_counting_iterator<const Box*>(&boxes1[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes1[0]+boxes1.size()));
  std::vector<const Box*> box2_ptr(boost::make_counting_iterator<const Box*>(&boxes2[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes2[0]+boxes2.size()));

  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;
  GeomTraits gt = boost::choose_param(get_param(np, internal_np::geom_traits), GeomTraits());

  CGAL::internal::Intersect_face_polylines<TM,
      GeomTraits,
      Box,
      PolylineRange,
      OutputIterator,
      VertexPointMap>
      Intersect_face_polyline(mesh,
                               faces,
                               polyline_range,
                               out,
                               vpmap,
                               gt);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           Intersect_face_polyline, cutoff);
  return Intersect_face_polyline.m_iterator;
}

// Note this is not officially documented
/*
 * detects and records intersections between two polylines.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * \attention If a polyline vertex intersects another polyline, the intersection will
 * be reported twice (even more if it is on a vertex).
 * \tparam Polyline a `RandomAccessRange` of points.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 *  elements in their respective range. This position is the index of the segment that holds the
 * intersection, so it is the index of the first point of the segment following the range order.
 * \tparam Kernel a model of `Kernel`
 *
 * \param polyline1 the first polyline to check for intersections.
 * \param polyline2 the second polyline to check for intersections.
 * \param out output iterator to be filled with all pairs of segments that intersect
 * \param K an instance of `Kernel`
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template < class Polyline,
           class OutputIterator,
           class Kernel>
OutputIterator
compute_polyline_polyline_intersection(const Polyline& polyline1,
                                       const Polyline& polyline2,
                                       OutputIterator out,
                                       const Kernel& K)
{
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t> Box;
  typedef typename Kernel::Point_3 Point;
  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(polyline1), boost::end(polyline1) ) - 1
        );

  boxes2.reserve(
        std::distance( boost::begin(polyline2), boost::end(polyline2) ) - 1
        );

  for(std::size_t i =0; i< polyline1.size()-1; ++i)
  {
    const Point& p1 = polyline1[i];
    const Point& p2 = polyline1[i+1];
    boxes1.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  for(std::size_t i =0; i< polyline2.size()-1; ++i)
  {
    const Point& p1 = polyline2[i];
    const Point& p2 = polyline2[i+1];
    boxes2.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr(boost::make_counting_iterator<const Box*>(&boxes1[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes1[0]+boxes1.size()));
  std::vector<const Box*> box2_ptr(boost::make_counting_iterator<const Box*>(&boxes2[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes2[0]+boxes2.size()));


  // compute intersections filtered out by boxes

  CGAL::internal::Intersect_polylines<Polyline,
      Kernel,
      Box,
      OutputIterator>
      intersect_polylines(polyline1,
                          polyline2,
                          out,
                          K);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           intersect_polylines, cutoff);
  return intersect_polylines.m_iterator;
}

// Note this is not officially documented
/*
 * detects and records intersections between two ranges of polylines.
 *  \attention If a polyline vertex intersects another polyline, the intersection will
 * be reported twice (even more if it is on a vertex).
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t> >`.
 * Each pair holds the index of the face and a pair containing the index of the polyline in the range and the index of
 * the first point of the segment in the polyline.
 * \tparam Kernel a model of `Kernel`
 *
 * \param polylines1 the first range of polylines to check for intersections.
 * \param polylines2 the second range of polylines to check for intersections.
 * \param out output iterator to be filled with all pairs of segments that intersect
 * \param K an instance of `Kernel`
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template < class PolylineRange,
           class OutputIterator,
           class Kernel>
OutputIterator
compute_polylines_polylines_intersection(const PolylineRange& polylines1,
                                         const PolylineRange& polylines2,
                                         OutputIterator out,
                                         const Kernel& K)
{
  //info.first is the index of the polyline in the range, info.second is the index of the point in the polyline
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::pair<std::size_t, std::size_t> > Box;
  typedef typename Kernel::Point_3 Point;
  typedef typename boost::range_value<PolylineRange>::type Polyline;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  int polylines_size = 0;
  BOOST_FOREACH(Polyline poly, polylines1)
  {
    polylines_size += std::distance( boost::begin(poly), boost::end(poly) ) -1;
  }
  boxes1.reserve( polylines_size );
  polylines_size = 0;
  BOOST_FOREACH(Polyline poly, polylines2)
  {
    polylines_size += std::distance( boost::begin(poly), boost::end(poly) ) -1;
  }
  boxes2.reserve(polylines_size);

  int range_size = std::distance( boost::begin(polylines1), boost::end(polylines1) );
  for(int j = 0; j < range_size; ++j)
  {
    Polyline poly = polylines1[j];
    int size = std::distance( boost::begin(poly), boost::end(poly) );
    for(int i =0; i< size - 1; ++i)
    {
      const Point& p1 = poly[i];
      const Point& p2 = poly[i+1];
      boxes1.push_back(Box(p1.bbox() + p2.bbox(), std::make_pair(j, i)));
    }
  }

  range_size = std::distance( boost::begin(polylines2), boost::end(polylines2) );
  for(int j = 0; j < range_size; ++j)
  {
    Polyline poly = polylines2[j];
    int size = std::distance( boost::begin(poly), boost::end(poly) );
    for(int i =0; i< size - 1; ++i)
    {
      const Point& p1 = poly[i];
      const Point& p2 = poly[i+1];
      boxes2.push_back(Box(p1.bbox() + p2.bbox(), std::make_pair(j, i)));
    }
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr(boost::make_counting_iterator<const Box*>(&boxes1[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes1[0]+boxes1.size()));
  std::vector<const Box*> box2_ptr(boost::make_counting_iterator<const Box*>(&boxes2[0]),
                                   boost::make_counting_iterator<const Box*>(&boxes2[0]+boxes2.size()));


  // compute intersections filtered out by boxes

  CGAL::internal::Intersect_polyline_ranges<PolylineRange,
      Kernel,
      Box,
      OutputIterator>
      intersect_polylines(polylines1,
                          polylines2,
                          out,
                          K);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           intersect_polylines, cutoff);
  return intersect_polylines.m_iterator;
}

// Note this is not officially documented
/*
 * reports all the pairs of faces intersecting between two triangulated surface meshes.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * @pre `CGAL::is_triangle_mesh(mesh1)`
 * @pre `CGAL::is_triangle_mesh(mesh2)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param mesh1 the first triangulated surface mesh to check for intersections
 * \param mesh2 the second triangulated surface mesh to check for intersections
 * \param out output iterator to be filled with all pairs of faces that intersect
 * \param np1 optional sequence of \ref namedparameters for `mesh1`, among the ones listed below
 * \param np2 optional sequence of \ref namedparameters for `mesh2`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh,
          class OutputIterator,
          class NamedParameters1,
          class NamedParameters2>
OutputIterator
compute_face_face_intersection(const TriangleMesh& mesh1,
                               const TriangleMesh& mesh2,
                               OutputIterator out,
                               const NamedParameters1& np1,
                               const NamedParameters2& np2)
{
  return compute_face_face_intersection(faces(mesh1), faces(mesh2),
                                        mesh1, mesh2, out, np1, np2);
}

// Note this is not officially documented
/*
 * detects and records intersections between a triangulated surface mesh
 * and a polyline.
 *  \attention If a polyline vertex intersects a face or another polyline, the intersection will
 * be reported twice (even more if it is on a vertex, edge, or point).
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam Polyline a `RandomAccessRange` of points. The point type of the range must be the
 * same as the value type of the vertex point map.
 * \cgalDescribePolylineType
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 *  elements in their respective range. In the case of the polyline, this position is the index
 * of the segment that holds the intersection, so it is the index of the first point of the
 * segment following the range order.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param mesh the triangulated surface mesh to check for intersections.
 * \param polyline the polyline to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np optional sequence of \ref namedparameters for `mesh`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh,
          class Polyline,
          class OutputIterator,
          class NamedParameters>
OutputIterator
compute_face_polyline_intersection(const TriangleMesh& mesh,
              const Polyline& polyline,
              OutputIterator out,
              const NamedParameters& np)
{
  return compute_face_polyline_intersection(faces(mesh), polyline, mesh, out, np);
}


// functions to check for overlap of meshes
template <class GT, class TriangleMesh, class VPM>
void get_one_point_per_cc(TriangleMesh& mesh,
                          const VPM& vpm,
                          std::vector<typename GT::Point_3>& points_of_interest)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  boost::unordered_map<face_descriptor, int> fid_map;
  int id = 0;
  BOOST_FOREACH(face_descriptor fd, faces(mesh))
  {
    fid_map.insert(std::make_pair(fd,id++));
  }
  boost::associative_property_map< boost::unordered_map<face_descriptor, int> >
      fid_pmap(fid_map);
  std::map<face_descriptor, int> fcc_map;

  int nb_cc = Polygon_mesh_processing::connected_components(mesh,
                                                            boost::make_assoc_property_map<std::map<face_descriptor, int> >(fcc_map),
                                                            Polygon_mesh_processing::parameters::face_index_map(fid_pmap));
  std::vector<bool> is_cc_treated(nb_cc, false);
  points_of_interest.resize(nb_cc);
  int cc_treated = 0;
  BOOST_FOREACH(face_descriptor fd, faces(mesh))
  {
    int cc=fcc_map[fd];
    if(!is_cc_treated[cc])
    {
      points_of_interest[cc]=get(vpm, target(halfedge(fd, mesh),mesh));
      is_cc_treated[cc] = true;
      if(++cc_treated == nb_cc)
        break;
    }
  }
}

//this assumes the meshes does not intersect
template <class TriangleMesh, class VPM, class GT, class AABB_tree>
bool is_mesh2_in_mesh1_impl(const AABB_tree& tree1,
                            const std::vector<typename GT::Point_3>& points_of_interest2,
                            const GT& gt)
{
  //for each CC, take a point on it and test bounded side
  Side_of_triangle_mesh<TriangleMesh, GT, VPM> sotm(tree1, gt);
  BOOST_FOREACH(const typename GT::Point_3& p, points_of_interest2)
  {
    if(sotm(p) == CGAL::ON_BOUNDED_SIDE) // sufficient as we know meshes do not intersect
    {
      return true;
    }
  }
  return false;
}

template <class TriangleMesh, class VPM, class GT>
bool is_mesh2_in_mesh1(const TriangleMesh& mesh_1,
                       const TriangleMesh& mesh_2,
                       const VPM& vpm1,
                       const VPM& vpm2,
                       const GT& gt)
{
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef CGAL::AABB_traits<GT, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABBTree;

  AABBTree tree1(faces(mesh_1).begin(), faces(mesh_1).end(), mesh_1, vpm1);
  std::vector<typename GT::Point_3> points_of_interest2;
  get_one_point_per_cc<GT>(mesh_2, vpm2, points_of_interest2);

  return is_mesh2_in_mesh1_impl<TriangleMesh, VPM>(tree1, points_of_interest2, gt);
}


}// namespace internal

namespace Polygon_mesh_processing{
/**
 * \ingroup PMP_predicates_grp
 * returns `true` if any pair of segments from `polylines1` and `polylines2` intersect, and `false` otherwise.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points.
 *         The point type must be from a \cgal Kernel.
 * \cgalDescribePolylineType
 *
 * @param polylines1 the first range of polylines to check for intersections.
 * @param polylines2 the second range of polylines to check for intersections.
 *
 */
template <class PolylineRange>
bool do_intersect(const PolylineRange& polylines1,
                  const PolylineRange& polylines2
#ifndef DOXYGEN_RUNNING
                  , const typename boost::enable_if<
                    typename boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<PolylineRange>,
                        boost::range_value<PolylineRange>,
                        boost::false_type >::type
                    >::type
                   >::type* = 0//end enable_if
#endif
    )
{
  typedef typename boost::range_value<PolylineRange>::type Polyline;
  typedef typename boost::range_value<Polyline>::type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;
  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    CGAL::internal::compute_polylines_polylines_intersection(polylines1, polylines2, OutputIterator(), K());
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_predicates_grp
 * returns `true` if any pair of segments from `polyline1` and `polyline2` intersect, and `false` otherwise.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \tparam Polyline a `RandomAccessRange` of points.
 *         The point type must be from a \cgal Kernel.
 * \cgalDescribePolylineType
 *
 * @param polyline1 the first polyline to check for intersections.
 * @param polyline2 the second polyline to check for intersections.
 *
 */
template <class Polyline>
bool do_intersect(const Polyline& polyline1,
                  const Polyline& polyline2
#ifndef DOXYGEN_RUNNING
                , const typename boost::enable_if<
                    typename boost::has_range_const_iterator<Polyline>::type
                  >::type* = 0,
                  const typename boost::disable_if<
                    typename boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<Polyline>,
                        boost::range_value<Polyline>,
                        boost::false_type
                      >::type
                    >::type
                  >::type* = 0//end enable_if
#endif
                 )
{
  typedef typename boost::range_value<Polyline>::type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;
  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    CGAL::internal::compute_polyline_polyline_intersection(polyline1, polyline2, OutputIterator(), K());
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_predicates_grp
 * returns `true` if any pair of faces from `mesh1` and `mesh2` intersect, and `false` otherwise.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * If `test_overlap` is `true`, bounded volume overlaps are tested as well. In that case, the meshes must be closed.
 *
 * @pre `CGAL::is_triangle_mesh(mesh1)`
 * @pre `CGAL::is_triangle_mesh(mesh2)`
 * @pre `!test_overlap || CGAL::is_closed(mesh1)`
 * @pre `!test_overlap || CGAL::is_closed(mesh2)`
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref namedparameters for `mesh1`
 * @tparam NamedParameters2 a sequence of \ref namedparameters for `mesh2`
 *
 * @param mesh1 the first triangulated surface mesh to check for intersections
 * @param mesh2 the second triangulated surface mesh to check for intersections
 * @param np1 optional sequence of \ref namedparameters for `mesh1`, among the ones listed below
 * @param np2 optional sequence of \ref namedparameters for `mesh2`, among the ones listed below
 * @param test_overlap if `true` tests also for mesh inclusions. Note that the inclusion tests do not depend on the orientation of the meshes.
 *                     if `false`, only the intersection of surface triangles are tested. Default is `false`.
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh1` (mesh2`).
 *   The two property map types must be the same.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool do_intersect(const TriangleMesh& mesh1,
                  const TriangleMesh& mesh2,
                  bool test_overlap,
                  const NamedParameters1& np1,
                  const NamedParameters2& np2)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh1));
  CGAL_precondition(CGAL::is_triangle_mesh(mesh2));
  CGAL_precondition(!test_overlap || CGAL::is_closed(mesh1));
  CGAL_precondition(!test_overlap || CGAL::is_closed(mesh2));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    CGAL::internal::compute_face_face_intersection(mesh1,mesh2, OutputIterator(), np1, np2);
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  if (test_overlap)
  {
    typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::const_type VertexPointMap1;
    typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::const_type VertexPointMap2;
    VertexPointMap1 vpm1 = boost::choose_param(get_param(np1, internal_np::vertex_point),
                                                get_const_property_map(boost::vertex_point, mesh1));
    VertexPointMap2 vpm2 = boost::choose_param(get_param(np2, internal_np::vertex_point),
                                                get_const_property_map(boost::vertex_point, mesh2));
    typedef typename GetGeomTraits<TriangleMesh, NamedParameters1>::type GeomTraits;
    GeomTraits gt = boost::choose_param(get_param(np1, internal_np::geom_traits), GeomTraits());

    return CGAL::internal::is_mesh2_in_mesh1(mesh1, mesh2, vpm1, vpm2, gt) ||
           CGAL::internal::is_mesh2_in_mesh1(mesh2, mesh1, vpm2, vpm1, gt);
  }
  return false;
}

//convenient overload
template <class TriangleMesh>
bool do_intersect(const TriangleMesh& mesh1,
                  const TriangleMesh& mesh2,
                  bool test_overlap = false,
                  const typename boost::disable_if<
                                    typename boost::has_range_const_iterator<TriangleMesh>::type
                  >::type* = 0)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh1));
  CGAL_precondition(CGAL::is_triangle_mesh(mesh2));
  return do_intersect(mesh1, mesh2, test_overlap, parameters::all_default(), parameters::all_default());
}

/**
 * \ingroup PMP_predicates_grp
 * returns `true` if any pair of face and segment from `mesh` and `polylines` intersect, and `false` otherwise.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points. The point type of the range must be the
 *  same as the value type of the vertex point map.
 * \cgalDescribePolylineType
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param mesh the triangulated surface mesh to check for intersections
 * @param polylines the range of polylines to check for intersections.
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class PolylineRange,
          class NamedParameters>
bool do_intersect(const TriangleMesh& mesh,
                  const PolylineRange& polylines,
                  const NamedParameters& np
#ifndef DOXYGEN_RUNNING
                , const typename boost::enable_if<
                    typename boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<PolylineRange>,
                        boost::range_value<PolylineRange>,
                        boost::false_type
                      >::type
                    >::type
                  >::type* = 0//end enable_if
#endif
                 )
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    CGAL::internal::compute_face_polylines_intersection(faces(mesh), polylines, mesh, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_predicates_grp
 * returns `true` if any pair of face and segment from `mesh` and `polyline` intersect, and `false` otherwise.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam Polyline a `RandomAccessRange` of points. The point type of the range must be the
 *  same as the value type of the vertex point map.
 * \cgalDescribePolylineType
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param mesh the triangulated surface mesh to check for intersections
 * @param polyline the polyline to check for intersections.
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class Polyline,
          class NamedParameters>
bool do_intersect(const TriangleMesh& mesh,
                  const Polyline& polyline,
                  const NamedParameters& np
#ifndef DOXYGEN_RUNNING
                , const typename boost::disable_if<
                    typename boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<Polyline>,
                        boost::range_value<Polyline>,
                        boost::false_type
                      >::type
                    >::type
                  >::type* = 0
#endif
                 )
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    CGAL::internal::compute_face_polyline_intersection(mesh,polyline, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

template <class TriangleMesh,
          class PolylineRange>
bool do_intersect(const TriangleMesh& mesh,
                  const PolylineRange& polylines,
                  const typename boost::enable_if<
                    typename boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<PolylineRange>,
                        boost::range_value<PolylineRange>,
                        boost::false_type
                      >::type
                    >::type
                  >::type* = 0)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  return do_intersect(mesh, polylines, parameters::all_default());
}


template <class TriangleMesh,
          class Polyline>
bool do_intersect(const TriangleMesh& mesh,
                  const Polyline& polyline,
                  const typename boost::disable_if<
                    typename boost::has_range_const_iterator<TriangleMesh>::type
                  >::type* = 0,
                  const typename boost::disable_if<
                    typename boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<Polyline>,
                        boost::range_value<Polyline>,
                        boost::false_type
                      >::type
                    >::type
                  >::type* = 0)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  return do_intersect(mesh, polyline, parameters::all_default());
}
}//end PMP
namespace internal{

template<class TriangleMeshRange,
         typename OutputIterator,
         class NamedParametersRange>
struct Mesh_callback
{
  typedef typename boost::range_value<TriangleMeshRange>::type TriangleMesh;
  typedef typename boost::range_value<NamedParametersRange>::type NamedParameter;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameter>::type GT;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameter>::const_type VPM;
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef CGAL::AABB_traits<GT, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABBTree;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  // fill them in the operator and test for inclusion with
  // Side_of_triangle_mesh.
  const TriangleMeshRange& meshes;
  OutputIterator m_iterator;
  const bool report_overlap;
  const NamedParametersRange& nps;
  std::vector<AABBTree*> trees;

  std::vector<std::vector<CGAL::Point_3<GT> > > points_of_interest;

  Mesh_callback(const TriangleMeshRange& meshes,
                OutputIterator iterator,
                const bool report_overlap,
                const NamedParametersRange& nps)
    : meshes(meshes), m_iterator(iterator),
      report_overlap(report_overlap), nps(nps)
  {
    int size = std::distance(meshes.begin(), meshes.end());
    trees = std::vector<AABBTree*>(size, NULL);
    points_of_interest.resize(size);
  }

  ~Mesh_callback()
  {
    BOOST_FOREACH(AABBTree* tree, trees)
    {
      delete tree;
    }
  }

  template<class TriangleMesh,
           class VPM,
           class GT>
  bool is_mesh2_in_mesh1(const TriangleMesh& mesh_1,
                         const TriangleMesh& mesh_2,
                         const int mesh_id_1,
                         const int mesh_id_2,
                         const VPM& vpm1,
                         const VPM& vpm2,
                         const GT& gt)
  {
    //test if mesh2 is included in mesh1

    //get AABB_tree for mesh1
    if(!trees[mesh_id_1])
    {
      trees[mesh_id_1] = new AABBTree(faces(mesh_1).begin(),
                                      faces(mesh_1).end(),
                                      mesh_1, vpm1);
    }
    //get a face-index map for mesh2
    if(points_of_interest[mesh_id_2].size() == 0)
      get_one_point_per_cc<GT>(mesh_2, vpm2, points_of_interest[mesh_id_2]);

    //test if mesh_2 is included in mesh_1:
    return is_mesh2_in_mesh1_impl<TriangleMesh, VPM, GT>(
      *trees[mesh_id_1],  points_of_interest[mesh_id_2], gt);
  }

  template<class Mesh_box>
  void operator()(const Mesh_box* b1, const Mesh_box* b2)
  {
    std::size_t mesh_id_1 = std::distance(meshes.begin(), b1->info());
    std::size_t mesh_id_2 = std::distance(meshes.begin(), b2->info());


    VPM vpm1 = boost::choose_param(get_param(*(nps.begin() + mesh_id_1), internal_np::vertex_point),
                                   get_const_property_map(CGAL::vertex_point, *b1->info()));

    VPM vpm2 = boost::choose_param(get_param(*(nps.begin() + mesh_id_2), internal_np::vertex_point),
                                   get_const_property_map(CGAL::vertex_point, *b2->info()));

    GT gt = boost::choose_param(get_param(*nps.begin(), internal_np::geom_traits), GT());

    //surfacic test
    if(Polygon_mesh_processing::do_intersect(*b1->info(),
                                             *b2->info(),
                                             false, // we do not put report_overlap here since
                                                    // we want to cache aabb-trees and reference points
                                             Polygon_mesh_processing::parameters::vertex_point_map(vpm1)
                                             .geom_traits(gt),
                                             Polygon_mesh_processing::parameters::vertex_point_map(vpm2)
                                             .geom_traits(gt)))
    {
      *m_iterator++ = std::make_pair(mesh_id_1, mesh_id_2);
    }
    //volumic test
    else if(report_overlap)
    {
      if(!CGAL::do_overlap(b1->bbox(), b2->bbox()))
        return;
      if(is_mesh2_in_mesh1(*b1->info(), *b2->info(), mesh_id_1, mesh_id_2, vpm1, vpm2, gt))
        *m_iterator++ = std::make_pair(mesh_id_1, mesh_id_2);
      else if(is_mesh2_in_mesh1(*b2->info(), *b1->info(), mesh_id_2, mesh_id_1, vpm2, vpm1, gt))
        *m_iterator++ = std::make_pair(mesh_id_2, mesh_id_1);
    }
  }
};
}//end internal

namespace Polygon_mesh_processing{
/*!
 * \ingroup PMP_predicates_grp
 * detects and reports all the pairs of meshes intersecting in a range of triangulated surface meshes.
 * A pair of meshes intersecting is put in the output iterator `out` as a `std::pair<std::size_t, std::size_t>`,
 * each index refering to the index of the triangle mesh in the input range.
 * If `report_overlap` is `true`, bounded volume overlaps are tested as well. In that case, the meshes must be closed.
 *
 * \tparam TriangleMeshRange a model of `Bidirectional Range` of triangulated surface meshes model of `FaceListGraph`.
 * \tparam OutputIterator an output iterator in which `std::pair<std::size_t, std::size_t>` can be put.
 * \tparam NamedParametersRange a range of named parameters.
 *
 * \param range the range of triangulated surface meshes to be checked for intersections.
 * \param out output iterator used to collect pairs of intersecting meshes.
 * \param report_overlap if `true` reports also mesh inclusions. Note that the inclusion tests do not depend on the orientation of the meshes.
 *                       if `false`, only the intersection of surface triangles are tested.
 *                       Default is `false`.
 * \param nps an optional range of `vertex_point_map` namedparameters containing the `VertexPointMap` of each mesh in `range`, in the same order.
 * The first named parameter of the range may contain a custom `geom_traits` for the whole range of meshes.
 * if this parameter is omitted, then an internal property map for
 *   `CGAL::vertex_point_t` should be available for every mesh in the range.
 * All the vertex point maps must be of the same type.
 * The default value for `geom_traits` is `CGAL::Kernel_traits<Point>::Kernel`, where `Point` is the
 * value type of the verter point map.
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 */

template <class TriangleMeshRange,
          class OutputIterator,
          class NamedParametersRange>
OutputIterator intersecting_meshes(const TriangleMeshRange& range,
                                         OutputIterator out,
                                   const bool report_overlap,
                                         NamedParametersRange nps)
{
  typedef typename TriangleMeshRange::const_iterator TriangleMeshIterator;

  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, TriangleMeshIterator> Mesh_box;
  std::vector<Mesh_box> boxes;
  boxes.reserve(std::distance(range.begin(), range.end()));

  for(TriangleMeshIterator it = range.begin(); it != range.end(); ++it)
  {
    boxes.push_back( Mesh_box(Polygon_mesh_processing::bbox(*it), it) );
  }

  std::vector<Mesh_box*> boxes_ptr(
        boost::make_counting_iterator(&boxes[0]),
    boost::make_counting_iterator(&boxes[0]+boxes.size()));

  //get all the pairs of meshes intersecting (no strict inclusion test)
  std::ptrdiff_t cutoff = 2000;
  CGAL::internal::Mesh_callback<TriangleMeshRange, OutputIterator, NamedParametersRange> callback(range, out, report_overlap, nps);
  CGAL::box_self_intersection_d(boxes_ptr.begin(), boxes_ptr.end(),
                                callback, cutoff);
  return callback.m_iterator;
}

template <class TriangleMeshRange, class OutputIterator>
OutputIterator intersecting_meshes(const TriangleMeshRange& range,
                                          OutputIterator out,
                                          bool report_overlap = false)
{
  std::vector<pmp_bgl_named_params<bool, internal_np::all_default_t> >nps;
  nps.reserve(std::distance(range.begin(), range.end()));
  for(int i = 0; i<std::distance(range.begin(), range.end()); ++i)
  {
    nps.push_back(parameters::all_default);
  }
  return intersecting_meshes(range, out, report_overlap, nps);
}

/**
 * \ingroup PMP_corefinement_grp
 * computes the intersection of triangles of `tm1` and `tm2`. The output is a
 * set of polylines with all vertices but endpoints being of degree 2.
 *
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref namedparameters
 * @tparam NamedParameters2 a sequence of \ref namedparameters
 * @tparam OutputIterator an output iterator in which `std::vector` of points
 *                        can be put. The point type is the one from the
 *                        vertex property map
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param polyline_output output iterator of polylines. Each polyline will be
 *        given as a vector of points
 * @param throw_on_self_intersection if `true`, for each input triangle mesh,
 *        the set of triangles closed to the intersection of `tm1` and `tm2` will be
 *        checked for self-intersection and `CGAL::Corefinement::Self_intersection_exception`
 *        will be thrown if at least one is found.
 * @param np1 optional sequence of \ref namedparameters among the ones listed below
 * @param np2 optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *    a property map with the points associated to the vertices of `tm1`
 *    (`tm2`). The two property map types must be the same.
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
template <class OutputIterator,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2 >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const NamedParameters1& np1,
                     const NamedParameters2& np2,
                     const bool throw_on_self_intersection=false)
{
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::const_type Vpm;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters2>::const_type Vpm2;
  CGAL_USE_TYPE(Vpm2);
  CGAL_assertion_code(
    static const bool same_vpm = (boost::is_same<Vpm,Vpm2>::value);)
  CGAL_static_assertion(same_vpm);

  Vpm vpm1 = choose_const_pmap(get_param(np1, internal_np::vertex_point),
                               tm1,
                               vertex_point);
  Vpm vpm2 = choose_const_pmap(get_param(np2, internal_np::vertex_point),
                               tm2,
                               vertex_point);

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh,Vpm>
    functor(tm1, tm2, vpm1, vpm2);
  functor(polyline_output, throw_on_self_intersection, true);
  return polyline_output;
}

template <class OutputIterator,
          class TriangleMesh >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const bool throw_on_self_intersection=false)
{
  return surface_intersection(tm1, tm2, polyline_output,
    CGAL::Polygon_mesh_processing::parameters::all_default(),
    CGAL::Polygon_mesh_processing::parameters::all_default(),
    throw_on_self_intersection
  );
}

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
