// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_FLATNESS_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_FLATNESS_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/utility.h>

#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

template <typename TriangleMesh, typename VPM>
void output_mesh(const TriangleMesh& tmesh, const VPM& vpm, std::ofstream& out)
{
  if(!tmesh.is_valid())
    return;

  out.precision(17);
  write_off(out, tmesh, CGAL::parameters::vertex_point_map(vpm));
  out.close();
}

template <typename TriangleMesh>
void output_mesh(const TriangleMesh& tmesh, std::ofstream& out)
{
  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type CVPM;

  CVPM cvpm = get_const_property_map(CGAL::vertex_point, tmesh);
  return output_mesh(tmesh, cvpm, out);
}

template <typename TriangleMesh>
void output_mesh(const TriangleMesh& tmesh, const std::string& filename)
{
  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type CVPM;

  std::ofstream out(filename);
  CVPM cvpm = get_const_property_map(CGAL::vertex_point, tmesh);
  output_mesh(tmesh, cvpm, out);
}

template <typename TriangleMesh, typename VPM>
void output_mesh(const TriangleMesh& tmesh, const VPM& vpm, std::stringstream& oss)
{
  std::cout << "[[ ---------- Dumping mesh: " << oss.str() << " ---------- ]]" << std::endl;

  std::ofstream out(oss.str().c_str());
  output_mesh(tmesh, vpm, out);
}

template <typename TriangleMesh>
void output_mesh(const TriangleMesh& tmesh, std::stringstream& oss)
{
  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type CVPM;

  CVPM cvpm = get_const_property_map(CGAL::vertex_point, tmesh);
  return output_mesh(tmesh, cvpm, oss);
}

template <typename TriangleMesh, typename AABBTree, typename VPM, typename GT>
double lfs_at_face(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                   const AABBTree& tree,
                   const TriangleMesh& tmesh,
                   const VPM& vpm,
                   const GT& gt)
{
  typedef typename GT::Point_3                                                           Point;
  typedef typename GT::Segment_3                                                         Segment;
  typedef typename GT::Vector_3                                                          Vector;
  typedef typename GT::Triangle_3                                                        Triangle;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor                halfedge_descriptor;

  typedef typename AABBTree::Primitive_id                                                Primitive_id;

  halfedge_descriptor h = halfedge(f, tmesh);

  // compute face centroid
  const Point p0 = get(vpm, source(h, tmesh));
  const Point p1 = get(vpm, target(h, tmesh));
  const Point p2 = get(vpm, target(next(h, tmesh), tmesh));
  const Point c = CGAL::centroid(p0, p1, p2);

  // make tiny segment
#ifdef CGAL_PMP_REPAIR_FLATNESS_COMPUTE_LFS_APPROXIMATION
  Vector n = Polygon_mesh_processing::compute_face_normal(f, tmesh);
#else
  Vector n = 0.015 * Polygon_mesh_processing::compute_face_normal(f, tmesh); // @fixme hardcoded
#endif

  Segment s(c - n, c + n);

  // intersect segment with mesh
  std::list<Primitive_id> primitives;
  tree.all_intersected_primitives(s, std::back_inserter(primitives));
  std::cout << primitives.size() << " intersected primitives" << std::endl;
  std::cout << tree.size() << " tree!!!" << std::endl;

  double lfs = std::numeric_limits<double>::max();

  for(const auto & pf : primitives)
  {
    if(pf == f)
      continue;

#ifdef CGAL_PMP_REPAIR_FLATNESS_COMPUTE_LFS_APPROXIMATION
    halfedge_descriptor ph = halfedge(pf, tmesh);
    const Point op0 = get(vpm, source(ph, tmesh));
    const Point op1 = get(vpm, target(ph, tmesh));
    const Point op2 = get(vpm, target(next(ph, tmesh), tmesh));
    const Triangle tr(op0, op1, op2);

    CGAL::Object obj = gt.intersect_3_object()(tr, s);

    Point inter_pt;
    if(gt.assign_3_object()(inter_pt, obj))
    {
      std::cout << "c/inter " << c << " " << inter_pt << std::endl;

      lfs = (std::min)(lfs, CGAL::approximate_sqrt(CGAL::squared_distance(c, inter_pt)));
    }
#else
    lfs = 0.;
#endif
  }

  std::cout << "face: " << f << " lfs: " << lfs << std::endl;
  return lfs;
}

#ifdef CGAL_PMP_REPAIR_FLATNESS_INITIALIZE_FROM_DIHEDRAL_ANGLES
template <typename TriangleMesh, typename MarkPmap, typename AABBTree>
void grow_selection(const typename boost::graph_traits<TriangleMesh>::face_descriptor seed_f,
                    const MarkPmap& marks,
                    const AABBTree& tree,
                    const TriangleMesh& tmesh)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor               halfedge_descriptor;

  do
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_face(seed_f, tmesh))
    {
      // @todo
    }
  }
  while();
}
#endif

template <typename TriangleMesh, typename GT>
void iterative_mesh_snap(TriangleMesh& tmesh,
                         const GT& gt)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor               halfedge_descriptor;

  typedef typename GT::FT                                                               FT;
  typedef typename TriangleMesh::template Property_map<vertex_descriptor, FT>           VertexFTMap;

  FT tolerance(1e-5); // @fixme hardcoded
  FT max_tolerance(0.1);
  int snap_id = 0;

  while(tolerance < max_tolerance)
  {
    std::cout << "Snapping with tolerance " << tolerance << std::endl;

    VertexFTMap tol_map = tmesh.template add_property_map<vertex_descriptor, FT>("v:tol", 0).first;

    // First, try to resolve everything within boundary cycles
    std::vector<halfedge_descriptor> boundaries;
    PMP::extract_boundary_cycles(tmesh, std::back_inserter(boundaries));
    std::cout << "Snapping " << boundaries.size() << " border(s)" << std::endl;

    for(halfedge_descriptor hd : boundaries)
    {
      std::cout << "boundary starting at : " << edge(hd, tmesh) << std::endl;

      std::vector<halfedge_descriptor> border;
      halfedge_descriptor done = hd;
      do
      {
        border.push_back(hd);
//        std::cout << tmesh.point(source(hd, tmesh)) << " " << tmesh.point(target(hd, tmesh)) << std::endl;
//        std::cout << "\t dist: " << CGAL::squared_distance(tmesh.point(source(hd, tmesh)),
//                                                           tmesh.point(target(hd, tmesh))) << std::endl;
        hd = next(hd, tmesh);
      }
      while(hd != done);

      CGAL_assertion(border.size() == PMP::internal::border_length(hd, tmesh));

      std::cout << "border of length: " << border.size() << std::endl;
      CGAL_assertion(border.size() >= 3);

      // Uses 'tolerance' at each vertex but bounds it by a % of the smallest incident edge length
      PMP::internal::assign_tolerance_with_local_edge_length_bound(border, tol_map, tolerance, tmesh);
      PMP::experimental::snap_vertex_range_onto_vertex_range_non_conforming(border, tmesh, border, tmesh, tol_map);
      PMP::stitch_boundary_cycle(hd, tmesh);
    }

    // Below this is basically 'PMP::border_halfedges()', but checking is_border on the whole mesh,
    // and not border of the face range
    std::vector<halfedge_descriptor> border_vertices;
    PMP::border_halfedges(tmesh, std::back_inserter(border_vertices));
    std::cout << border_vertices.size() << " total border halfedges" << std::endl;

    for(halfedge_descriptor hd : border_vertices)
    {
      CGAL_assertion(is_border(hd, tmesh));
//      std::cout << tmesh.point(source(hd, tmesh)) << " " << tmesh.point(target(hd, tmesh)) << std::endl;
    }

    std::cout << std::endl << "Snapping ALL borders" << std::endl;

    // re-assign tolerance values since we have moved vertices
    PMP::internal::assign_tolerance_with_local_edge_length_bound(border_vertices, tol_map, tolerance, tmesh);
    PMP::experimental::snap_vertex_range_onto_vertex_range_non_conforming(border_vertices, tmesh,
                                                                          border_vertices, tmesh,
                                                                          tol_map,
                                                                          parameters::all_default(),
                                                                          parameters::all_default());
    CGAL_assertion(is_valid_polygon_mesh(tmesh));

    PMP::remove_degenerate_faces(tmesh);
    PMP::stitch_borders(tmesh);

    std::stringstream oss;
    oss << "results/snap_intermediate_" << snap_id++ << "_" << tolerance << ".off" << std::ends;
    output_mesh(tmesh, oss);

    tolerance *= 10;
  }
}

template <typename HalfedgeRange, typename TriangleMesh, typename GT>
bool fill_hole_with_polyline(const HalfedgeRange& border_loop,
                             TriangleMesh& tmesh,
                             TriangleMesh& hole_mesh,
                             const GT& gt)
{
  typedef typename GT::Point_3                                                                  Point;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor                         vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor                       halfedge_descriptor;

  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::type       VPM;
  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type CVPM;

  CVPM cvpm = get_const_property_map(CGAL::vertex_point, tmesh);

  std::vector<Point> hole_points, third_points;
  hole_points.reserve(border_loop.size());
  third_points.reserve(border_loop.size());

  for(halfedge_descriptor hd : border_loop)
  {
    vertex_descriptor vd = source(hd, tmesh);
    hole_points.push_back(get(cvpm, vd));
    third_points.push_back(get(cvpm, target(next(opposite(hd, tmesh), tmesh), tmesh)));
  }
  CGAL_assertion(hole_points.size() >= 3);

//  std::cout << "fill_hole_with_polyline with hole_points.size() = " << hole_points.size() << std::endl;

  std::vector<CGAL::Triple<int, int, int> > patch;
  Polygon_mesh_processing::triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch));

  std::vector<vertex_descriptor> hole_vertices;
  hole_vertices.reserve(hole_points.size());
  VPM vpm = get_property_map(CGAL::vertex_point, hole_mesh);
  for(const Point& pt : hole_points)
  {
    vertex_descriptor vd = add_vertex(hole_mesh);
    hole_vertices.push_back(vd);
    put(vpm, vd, pt);
  }

  for(const CGAL::Triple<int, int, int>& face : patch)
  {
    std::vector<vertex_descriptor> vertex_range {{ hole_vertices[face.first],
                                                   hole_vertices[face.second],
                                                   hole_vertices[face.third] }};
    Euler::add_face(vertex_range, hole_mesh);
  }

  Polygon_mesh_processing::stitch_borders(hole_mesh);

  return tmesh.is_valid();
}

template <typename TriangleMesh, typename GT>
bool fill_hole(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor bhd,
               TriangleMesh& tmesh,
               const GT& gt)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor                       halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor                           face_descriptor;

  typedef typename GT::Point_3                                                                  Point;

  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::type       VPM;
  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type CVPM;

  std::cout << "Fill hole " << PMP::number_of_borders(tmesh) << " left" << std::endl;

  bool is_healthy_hole = true;
  std::set<Point> hole_points;
  std::set<Point> non_manifold_points; // @todo optimize by caching non-manifold vertices

  CVPM cvpm = get_const_property_map(CGAL::vertex_point, tmesh);
  std::size_t hole_size = 0;

  std::vector<face_descriptor> local_new_faces;

  halfedge_descriptor iter_bhd = bhd, done = bhd;
  do
  {
    const Point& hole_pt = get(cvpm, target(iter_bhd, tmesh));
    if(!hole_points.insert(hole_pt).second)
    {
      is_healthy_hole = false;
      non_manifold_points.insert(hole_pt);
      std::cout << "non-manifold vertex at: " << hole_pt << std::endl;
    }

    ++hole_size;
    iter_bhd = next(iter_bhd, tmesh);
  }
  while(iter_bhd != done);

  CGAL_assertion(hole_size >= 3); // otherwise it should have been stitched before

  std::cout << "hole of size: " << hole_points.size() << " (health: " << is_healthy_hole << ")" << std::endl;

  if(is_healthy_hole)
  {
    // can simply fill normally
    PMP::triangulate_hole(tmesh, bhd, std::back_inserter(local_new_faces));
  }
  else
  {
    return false; // @tmp

    // Non-manifold vertex on the hole: split the hole into multiple holes
    std::cout << non_manifold_points.size() << " non-manifold vertices on the border" << std::endl;

    std::deque<halfedge_descriptor> hole_border;
    iter_bhd = bhd;
    do
    {
      hole_border.push_back(iter_bhd);
      iter_bhd = next(iter_bhd, tmesh);
    }
    while (iter_bhd != done);

    std::vector<TriangleMesh> hole_fillers;
    hole_fillers.reserve(non_manifold_points.size());

    std::stack<std::deque<halfedge_descriptor> > border_loops;
    border_loops.push(hole_border);

    while(!border_loops.empty())
    {
      std::cout << "stack has size: " << border_loops.size() << std::endl;

      const std::deque<halfedge_descriptor> border_loop = border_loops.top(); // @todo without copy
      border_loops.pop();

      std::cout << "border loop has size: " << border_loop.size() << std::endl;

      CGAL_assertion(border_loop.size() > 1);
      if(border_loop.size() == 2)
      {
        std::cerr << "Something is wrong... This loop should have been stitched" << std::endl;
        PMP::stitch_borders(tmesh); // just to see the error given by stitch()
        std::cin.get();
//        CGAL_assertion(false);
        return false;
      }

      // Find any non-manifold vertex
      typedef typename std::deque<halfedge_descriptor>::const_iterator           VHIT;

      VHIT halfedge_with_non_manifold_source_it;
      VHIT vhit = border_loop.cbegin(), last = std::prev(border_loop.cend()), vhend = border_loop.cend();
      for(; vhit!=vhend; ++vhit)
      {
        halfedge_descriptor hd = *vhit;
        if(non_manifold_points.find(get(cvpm, target(hd, tmesh))) == non_manifold_points.end())
          continue;

        std::cout << "found a non-manifold vertex at " << get(cvpm, target(hd, tmesh)) << std::endl;

        // move forward one so that source(vit, tmesh) is this non-manifold vertex
        if(vhit == last)
          halfedge_with_non_manifold_source_it = border_loop.cbegin();
        else
          halfedge_with_non_manifold_source_it = ++vhit;

        break;
      }

      // must have found a non-manifold vertex, even if we are in a simple loop
      CGAL_assertion(halfedge_with_non_manifold_source_it != VHIT());
      halfedge_descriptor hd_with_non_manifold_source = *halfedge_with_non_manifold_source_it;
      std::cout << "hd_with_non_manifold_source: " << hd_with_non_manifold_source << std::endl;

      vhit = halfedge_with_non_manifold_source_it;

      // Now, walk the border loop and split into left and right part
      bool walking_right_loop = true;
      int non_manifold_vertex_counter = 0;
      std::deque<halfedge_descriptor> left_border_loop;
      std::deque<halfedge_descriptor> right_border_loop;
      do
      {
        halfedge_descriptor hd = *vhit;
        std::cout << "at " << hd << std::endl;
        std::cout << "  [" << tmesh.point(source(hd, tmesh)) << " -- " << tmesh.point(target(hd, tmesh)) << "]" << std::endl;

        if(walking_right_loop)
          right_border_loop.push_back(hd);
        else
          left_border_loop.push_back(hd);

        // check if we have met another non-manifold vertex
        if(non_manifold_points.find(get(cvpm, target(hd, tmesh))) != non_manifold_points.end())
        {
          std::cout << "met non-manifold vertex: " << get(cvpm, target(hd, tmesh)) << std::endl;
          ++non_manifold_vertex_counter;

          // If we meet the same non-manifold point, end the right loop. Note that this
          // doesn't guarantee that the right loop is simple.
          if(get(cvpm, target(hd, tmesh)) == get(cvpm, source(hd_with_non_manifold_source, tmesh)))
          {
            // Note that we could meet the same non-manifold vertex thrice (or even more)
            std::cout << "switching to left loop" << std::endl;
            walking_right_loop = false;
          }
        }

        vhit = (vhit == last) ? border_loop.cbegin() : ++vhit;
      }
      while(vhit != halfedge_with_non_manifold_source_it);

      std::cout << "met " << non_manifold_vertex_counter << " non manifold vertices" << std::endl;

      // If the loop was simple (only one non-manifold vertex or two different non manifold vertices)
      if(non_manifold_vertex_counter <= 2 && left_border_loop.empty()) // @tmp
      {
        std::cout << "found simple loop with border_loop.size() = " << border_loop.size() << std::endl;

        hole_fillers.resize(hole_fillers.size() + 1);
        fill_hole_with_polyline(border_loop, tmesh, hole_fillers.back(), gt);
      }
      else
      {
        std::cout << left_border_loop.size() << " on the left" << std::endl;
        std::cout << right_border_loop.size() << " on the right" << std::endl;

        // the loop was not simple enough, split it and start again
        if(left_border_loop.empty())
        {
          // if we're unlucky, all the non-manifold vertices are on the same side of the loop, so
          // do some random perturbations (@todo something smart, see what was done in repair_PS...)
          right_border_loop.push_back(right_border_loop.front());
          right_border_loop.pop_front();
          border_loops.push(right_border_loop);
        }
        else
        {
          border_loops.push(right_border_loop);
          border_loops.push(left_border_loop);
        }
      }
    }

    for(const TriangleMesh& sm : hole_fillers)
    {
      static int i = 0;
      std::stringstream oss;
      oss << "results/hole_filler_" << i++ << ".off" << std::ends;
      output_mesh(sm, oss);
    }

    // Now, plug back in all the patches and stitch
    CGAL_assertion(hole_fillers.size() >= 2);
    while(hole_fillers.size() != 1)
    {
      copy_face_graph(hole_fillers.back(), hole_fillers.front());
      hole_fillers.resize(hole_fillers.size() - 1);
    }

    typedef std::pair<face_descriptor, face_descriptor>               Face_pair;

    std::vector<Face_pair> f2f;
    copy_face_graph(hole_fillers.front(), tmesh,
                    parameters::face_to_face_output_iterator(std::back_inserter(f2f)));

    for(const Face_pair& fp : f2f)
      local_new_faces.push_back(fp.second);

    PMP::stitch_borders(tmesh);

    output_mesh(tmesh, "results/non_manifold_hole_filled.off"); // @fixme static int name (others too)
  }

  if(PMP::does_self_intersect(local_new_faces, tmesh))
  {
    output_mesh(tmesh, "results/bad_fill.off");
    std::cerr << "Issues during hole filling" << std::endl; // @tmp
    return false;
  }

  return true;
}

template <typename TriangleMesh, typename GT>
bool fill_holes(TriangleMesh& tmesh,
                const GT& gt)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor             halfedge_descriptor;

  std::cout << "Hole filling! (" << PMP::number_of_borders(tmesh) << " borders)" << std::endl;

  CGAL_warning(!PMP::does_self_intersect(tmesh));

  // During the snapping, we merge duplicate points on the boundary, now we need to split them
  // again to subdivide non manifold holes into simpler, smaller, manifold holes
  PMP::duplicate_non_manifold_vertices(tmesh);

  // @todo this is very brute forcy, it should be restricted to the halfedges in the star
  // (2 ring in the offset mesh should be enough)
  unsigned int iter = 0;
  unsigned int nb = PMP::number_of_borders(tmesh);

  bool all_went_well = true;
  for(;;)
  {
    std::cout << "Hole filling iteration: " << iter << std::endl;

    bool filled_a_hole = false;

    for(halfedge_descriptor hd : halfedges(tmesh)) // @todo local
    {
      if(is_border(hd, tmesh))
      {
        if(!fill_hole(hd, tmesh, gt))
        {
          all_went_well = false;
          output_mesh(tmesh, "results/bad_fill.off");
        }

        filled_a_hole = true;
      }
    }

    if(!filled_a_hole)
      break;

    if(iter++ >= 0/*10 * nb*/) // safeguard in case we're continuously failing to fill a hole @tmp
      break;
  }

  return all_went_well && tmesh.is_valid();
}

template <typename TriangleMesh, typename NamedParameters>
void repair_flatness(TriangleMesh& tmesh,
                     const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor                 face_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                 GT;

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh>                      Primitive;
  typedef CGAL::AABB_traits<GT, Primitive>                                            Traits;
  typedef CGAL::AABB_tree<Traits>                                                     Tree;

  typedef CGAL::dynamic_face_property_t<bool>                                         Face_boolean_tag;
  typedef typename boost::property_map<TriangleMesh, Face_boolean_tag>::type          Face_boolean_map;

  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::type       VPM;
  typedef typename CGAL::property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type CVPM;

  CVPM cvpm = get_const_property_map(CGAL::vertex_point, tmesh);

  GT gt = choose_param(get_param(np, internal_np::geom_traits), GT());

  Tree tree(faces(tmesh).first, faces(tmesh).second, tmesh);
  tree.accelerate_distance_queries();

  Face_boolean_map marks = get(Face_boolean_tag(), tmesh);
  for(face_descriptor f : faces(tmesh))
    put(marks, f, false);

#ifdef CGAL_PMP_REPAIR_FLATNESS_INITIALIZE_FROM_DIHEDRAL_ANGLES // @todo (is there any point though?)
  // Detection helper
  typedef typename boost::property_map<FaceGraph, CGAL::face_patch_id_t<int> >::type  PatchID;

  boost::property_map<TriangleMesh, CGAL::edge_is_feature_t>::type eif = get(CGAL::edge_is_feature, tmesh);
  PatchID pid = get(CGAL::face_patch_id_t<int>(), tmesh);

  double angle = 178; // loop from 179 to 175
  Polygon_mesh_processing::sharp_edges_segmentation(tmesh, angle, eif, pid);

  // Grow selection
  for(edge_descriptor e : sharp_edges)
  {
    grow_selection(halfedge(e, tmesh), marks, tree, tmesh);
    grow_selection(opposite(halfedge(e, tmesh), tmesh), marks, tree, tmesh);
  }

#else
  for(face_descriptor f : faces(tmesh))
  {
    double lfs = lfs_at_face(f, tree, tmesh, cvpm, gt);

#ifdef CGAL_PMP_REPAIR_FLATNESS_COMPUTE_LFS_APPROXIMATION
    if(lfs < 0.005) // @fixme hardcoded
#else
    if(lfs == 0.) // in that case, we just put lfs to '0' if there is an intersection with a scaled normal
#endif
      put(marks, f, true);
  }
#endif

  // Remove all selected faces
  for(face_descriptor f : faces(tmesh))
  {
    if(get(marks, f))
    {
      std::cout << "Remove: " << f << std::endl;
      CGAL::Euler::remove_face(halfedge(f, tmesh), tmesh);
    }
  }
  std::cout << "post remove: " << faces(tmesh).size() << " total: " << num_faces(tmesh) << std::endl;


  output_mesh(tmesh, "results/punctured.off");

  // Snap
  iterative_mesh_snap(tmesh, gt);
  output_mesh(tmesh, "results/snapped.off");

  // Hole Fill
  fill_holes(tmesh, gt);
  output_mesh(tmesh, "results/filled.off");
}

template <typename TriangleMesh>
void repair_flatness(TriangleMesh& tmesh)
{
  return repair_flatness(tmesh, CGAL::parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H
