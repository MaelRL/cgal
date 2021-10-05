// Copyright (c) 2015-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/Repair/helper.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/repair_manifoldness.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTION_NO_POLYHEDRAL_ENVELOPE_CHECK
#include <CGAL/Polyhedral_envelope.h>
#endif

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/utility.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// #define CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING
// #define CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING
// #define CGAL_PMP_REMOVE_SELF_INTERSECTION_NO_POLYHEDRAL_ENVELOPE_CHECK

// Self-intersection removal is done by making a big-enough hole and filling it
//
// Local self-intersection removal is more subtle and only considers self-intersections
// within a connected component. It then tries to fix those by trying successively:
// - smoothing with the sharp edges in the area being constrained
// - smoothing without the sharp edges in the area being constrained
// - hole-filling with the sharp edges in the area being constrained
// - hole-filling without the sharp edges in the area being constrained
//
// The working area grows as long as we haven't been able to fix the self-intersection,
// up to a user-defined number of times.

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
static int unsolved_self_intersections = 0;
static int self_intersections_solved_by_constrained_smoothing = 0;
static int self_intersections_solved_by_unconstrained_smoothing = 0;
static int self_intersections_solved_by_constrained_hole_filling = 0;
static int self_intersections_solved_by_unconstrained_hole_filling = 0;
#endif

template <typename FaceRange, typename EdgeConstrainMap,
          typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
void constrain_sharp_and_border_edges(const FaceRange& faces,
                                      TriangleMesh& tmesh,
                                      EdgeConstrainMap& eif,
                                      const bool constrain_sharp_edges,
                                      const double dihedral_angle,
                                      const double /*weak_DA*/,
                                      VertexPointMap vpm,
                                      const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename GeomTraits::FT                                           FT;
  typedef typename GeomTraits::Vector_3                                     Vector;

  std::map<edge_descriptor, bool> is_border_of_selection;
  for(face_descriptor f : faces)
  {
    // @fixme what about nm vertices
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      // Default initialization is guaranteed to be `false`. Thus, meet it once will switch
      // the value to `true` and meeting it twice will switch back to `false`.
      const edge_descriptor e = edge(h, tmesh);
      is_border_of_selection[e] = !(is_border_of_selection[e]);
    }
  }

#if 0 // Until detect_features++ is integrated
  CGAL::Polygon_mesh_processing::experimental::detect_sharp_edges_pp(faces, tmesh, dihedral_angle, eif,
                                                                     parameters::weak_dihedral_angle(weak_DA));

  // borders are also constrained
  for(const auto& ep : is_border_of_selection)
    if(ep.second)
      put(eif, ep.first, true);
#else
  // this is basically the code that is in detect_features (at the very bottom)
  // but we do not want a folding to be marked as a sharp feature so the dihedral angle is also
  // bounded from above
  const double bound = dihedral_angle;
  const double cos_angle = std::cos(bound * CGAL_PI / 180.);

  for(const auto& ep : is_border_of_selection)
  {
    bool flag = ep.second;
    if(constrain_sharp_edges && !flag)
    {
      const halfedge_descriptor h = halfedge(ep.first, tmesh);
      CGAL_assertion(!is_border(edge(h, tmesh), tmesh));

      const face_descriptor f1 = face(h, tmesh);
      const face_descriptor f2 = face(opposite(h, tmesh), tmesh);

      // @todo cache normals
      const Vector n1 = compute_face_normal(f1, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
      const Vector n2 = compute_face_normal(f2, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
      const FT c = gt.compute_scalar_product_3_object()(n1, n2);

      // Do not mark as sharp edges with a dihedral angle that is almost `pi` because this is likely
      // due to a foldness on the mesh rather than a sharp edge that we wish to preserve
      // (Ideally this would be pre-treated as part of the flatness treatment)
      flag = (c <= cos_angle && c >= -cos_angle);
    }

    is_border_of_selection[ep.first] = flag; // Only needed for output, really
    put(eif, ep.first, flag);
  }
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/constrained_edges.polylines.txt");
  out << std::setprecision(17);
  for(edge_descriptor e : edges(tmesh))
    if(get(eif, e))
       out << "2 " << tmesh.point(source(e, tmesh)) << " " << tmesh.point(target(e, tmesh)) << std::endl;
  out.close();
#endif
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename PolyhedralEnvelope>
bool remove_self_intersections_with_smoothing(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& face_range,
                                              TriangleMesh& tmesh,
                                              const bool constrain_sharp_edges,
                                              const double dihedral_angle,
                                              const double weak_DA,
                                              const PolyhedralEnvelope& cc_envelope,
                                              VertexPointMap vpm,
                                              const GeomTraits& gt)
{
  namespace CP = CGAL::parameters;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: repair with smoothing... (constraining sharp edges: ";
  std::cout << std::boolalpha << constrain_sharp_edges << ")" << std::endl;
#endif

  CGAL_precondition(does_self_intersect(face_range, tmesh));

  // Rather than working directly on the mesh, copy a range and work on this instead
  const CGAL::Face_filtered_graph<TriangleMesh> ffg(tmesh, face_range);
  TriangleMesh local_mesh;
  CGAL::copy_face_graph(ffg, local_mesh, CP::vertex_point_map(vpm));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  CGAL::IO::write_polygon_mesh("results/local_mesh.off", local_mesh, CGAL::parameters::stream_precision(17));
#endif

  // Constrain sharp and border edges
  typedef CGAL::dynamic_edge_property_t<bool>                                 Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type EIFMap;
  EIFMap eif = get(Edge_property_tag(), local_mesh);

  VertexPointMap local_vpm = get_property_map(vertex_point, local_mesh);

  constrain_sharp_and_border_edges(faces(local_mesh), local_mesh, eif, constrain_sharp_edges,
                                   dihedral_angle, weak_DA, local_vpm, gt);

  // @todo choice of number of iterations? Till convergence && max of 100?
  Polygon_mesh_processing::smooth_mesh(faces(local_mesh), local_mesh, CP::edge_is_constrained_map(eif)
                                                                         .number_of_iterations(100)
                                                                         .use_safety_constraints(false));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  CGAL::IO::write_polygon_mesh("results/post_smoothing_local_mesh.off", local_mesh, CGAL::parameters::stream_precision(17));
#endif

  if(does_self_intersect(local_mesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: patch still self-intersecting after smoothing\n";
#endif
    return false;
  }
  if (!cc_envelope.is_empty() && !cc_envelope(local_mesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: patch is not in the input polyhedral envelope\n";
#endif
    return false;
  }

  std::vector<std::vector<Point> > patch;
  for(const face_descriptor f : faces(local_mesh))
  {
    halfedge_descriptor h = halfedge(f, local_mesh);
    patch.emplace_back(std::initializer_list<Point>{get(local_vpm, target(h, local_mesh)),
                                                    get(local_vpm, target(next(h, local_mesh), local_mesh)),
                                                    get(local_vpm, target(prev(h, local_mesh), local_mesh))});
  }

  std::set<face_descriptor> new_faces;
  bool success = replace_faces_with_patch(face_range, patch, tmesh, vpm, std::inserter(new_faces, new_faces.end()));
  if(!success)
    return false;

  CGAL_assertion(!does_self_intersect(new_faces, tmesh, parameters::vertex_point_map(vpm)));
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(constrain_sharp_edges)
    ++self_intersections_solved_by_constrained_smoothing;
  else
    ++self_intersections_solved_by_unconstrained_smoothing;
#endif
  return true;
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

// This function is only called when the hole is NOT subdivided into smaller holes
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename PolyhedralEnvelope>
bool fill_hole(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
               TriangleMesh& tmesh,
               const PolyhedralEnvelope& cc_envelope,
               VertexPointMap vpm,
               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Attempting hole-filling (no constraints), " << cc_faces.size() << " faces\n";
#endif

  if(!order_border_halfedge_range(cc_border_hedges, tmesh))
  {
    CGAL_assertion(false); // we shouldn't fail to orient the boundary cycle of the complete hole
    return false;
  }

  std::set<vertex_descriptor> cc_interior_vertices;
  std::set<edge_descriptor> cc_interior_edges;

  std::vector<vertex_descriptor> cc_border_vertices;
  cc_border_vertices.reserve(cc_border_hedges.size());

  std::vector<std::vector<Point> > patch;
  if(!construct_tentative_hole_patch(cc_border_vertices, cc_interior_vertices, cc_interior_edges,
                                     cc_border_hedges, cc_faces, patch, tmesh, vpm, gt) ||
     !check_patch_sanity<TriangleMesh>(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to find acceptable hole patch\n";
#endif

    return false;
  }

  if (!cc_envelope.is_empty() && !cc_envelope(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Patch is not inside the input polyhedral envelope\n";
#endif
    return false;
  }

  // Could renew the range directly within the patch replacement function
  // to avoid erasing and re-adding the same face
  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  // Plug the new triangles in the mesh, reusing previous edges and faces
  replace_faces_with_patch(cc_border_vertices, cc_interior_vertices,
                           cc_border_hedges, cc_interior_edges,
                           cc_faces, patch, tmesh, vpm,
                           std::inserter(working_face_range, working_face_range.end()));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  static int filed_hole_id = 0;
  std::stringstream oss;
  oss << "results/filled_basic_" << filed_hole_id++ << ".off" << std::ends;
  std::ofstream(oss.str().c_str()) << std::setprecision(17) << tmesh;
#endif

  CGAL_postcondition(is_valid_polygon_mesh(tmesh));

  return true;
}

// Same function as above but border of the hole is not known
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename PolyhedralEnvelope>
bool fill_hole(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
               TriangleMesh& tmesh,
               const PolyhedralEnvelope& cc_envelope,
               VertexPointMap vpm,
               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  std::vector<halfedge_descriptor> cc_border_hedges;
  for(face_descriptor fd : cc_faces)
  {
    halfedge_descriptor h = halfedge(fd, tmesh);
    for(int i=0; i<3; ++i)
    {
      if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
        cc_border_hedges.push_back(h);

      h = next(h, tmesh);
    }
  }

  if(order_border_halfedge_range(cc_border_hedges, tmesh))
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh,
                     cc_envelope,vpm, gt);
  else
    return false;
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename PolyhedralEnvelope>
bool fill_hole_with_constraints(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                TriangleMesh& tmesh,
                                const double dihedral_angle,
                                const double weak_DA,
                                const PolyhedralEnvelope& cc_envelope,
                                VertexPointMap vpm,
                                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Attempting local hole-filling with constrained sharp edges..." << std::endl;
#endif

  // If we are treating self intersections locally, first try to constrain sharp edges in the hole
  typedef CGAL::dynamic_edge_property_t<bool>                                 Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type EIFMap;
  EIFMap eif = get(Edge_property_tag(), tmesh);

  constrain_sharp_and_border_edges(cc_faces, tmesh, eif, true /*constrain_sharp_edges*/, dihedral_angle, weak_DA, vpm, gt);

  // Partition the hole using these constrained edges
  std::set<face_descriptor> visited_faces;
  std::vector<std::vector<Point> > patch;

  int cc_counter = 0;
  for(face_descriptor f : cc_faces)
  {
    if(!visited_faces.insert(f).second) // already visited that face
      continue;

    // gather the faces of the sub-hole
    std::set<face_descriptor> sub_cc;
    Polygon_mesh_processing::connected_component(f, tmesh, std::inserter(sub_cc, sub_cc.end()),
                                                 CGAL::parameters::edge_is_constrained_map(eif));

    visited_faces.insert(sub_cc.begin(), sub_cc.end());
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "CC of size " << sub_cc.size() << " (total: " << cc_faces.size() << ")" << std::endl;
#endif
    ++cc_counter;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    dump_cc(sub_cc, tmesh, "results/current_cc.off");
#endif

    // The mesh is not modified, but 'patch' gets filled
    if(!construct_tentative_sub_hole_patch(patch, sub_cc, cc_faces, tmesh, vpm, gt))
    {
      // Something went wrong while finding a potential cover for the a sub-hole --> use basic hole-filling
      return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh,
                       cc_envelope,vpm, gt);
    }
  }
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << cc_counter << " independent sub holes" << std::endl;
#endif
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/hole_fillers.off");
  out.precision(17);

  out << "OFF\n";
  out << 3*patch.size() << " " << patch.size() << " 0\n";

  for(const auto& f : patch)
  {
    for(const auto& pt : f)
      out << pt << "\n";
  }

  int id = 0;
  for(std::size_t i=0; i<patch.size(); ++i)
  {
    out << "3 " << id << " " << id+1 << " " << id+2 << "\n";
    id += 3;
  }
  out.close();
#endif

  // We're assembling multiple patches so we could have the same face appearing multiple times...
  if(!check_patch_sanity<TriangleMesh>(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Unhealthy patch, use base fill_hole" << std::endl;
#endif
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh,
                     cc_envelope, vpm, gt);
  }

  // check if the patch is inside the input polyhedral envelope
  if(!cc_envelope.is_empty() && !cc_envelope(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Patch is not entirely inside the input polyhedral envelope, use base fill_hole" << std::endl;
#endif
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh,
                     cc_envelope, vpm, gt);
  }

  // Plug the hole-filling patch in the mesh
  std::set<face_descriptor> new_faces;
  bool success = replace_faces_with_patch(cc_faces, patch, tmesh, vpm, std::inserter(new_faces, new_faces.end()));
  if(!success)
    return false;

  // Otherwise it should have failed the sanity check
  CGAL_assertion(!does_self_intersect(new_faces, tmesh, parameters::vertex_point_map(vpm)));

  // Update working range with the new faces
  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  working_face_range.insert(new_faces.begin(), new_faces.end());

  return true;
}

template <class Box, class TM, class VPM, class GT, class OutputIterator>
struct Strict_intersect_edges // "strict" as in "not sharing a vertex"
{
  typedef typename boost::graph_traits<TM>::halfedge_descriptor               halfedge_descriptor;
  typedef typename GT::Segment_3                                              Segment;

  mutable OutputIterator m_iterator;
  const TM& m_tmesh;
  const VPM m_vpmap;

  typename GT::Construct_segment_3 m_construct_segment;
  typename GT::Do_intersect_3 m_do_intersect;

  Strict_intersect_edges(const TM& tmesh, VPM vpmap, const GT& gt, OutputIterator it)
    :
      m_iterator(it),
      m_tmesh(tmesh),
      m_vpmap(vpmap),
      m_construct_segment(gt.construct_segment_3_object()),
      m_do_intersect(gt.do_intersect_3_object())
  {}

  void operator()(const Box* b, const Box* c) const
  {
    const halfedge_descriptor h = b->info();
    const halfedge_descriptor g = c->info();

    if(source(h, m_tmesh) == target(g, m_tmesh) || target(h, m_tmesh) == source(g, m_tmesh))
      return;

    const Segment s1 = m_construct_segment(get(m_vpmap, source(h, m_tmesh)), get(m_vpmap, target(h, m_tmesh)));
    const Segment s2 = m_construct_segment(get(m_vpmap, source(g, m_tmesh)), get(m_vpmap, target(g, m_tmesh)));

    if(m_do_intersect(s1, s2))
      *m_iterator++ = std::make_pair(b->info(), c->info());
  }
};

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool is_simple_3(const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                 const TriangleMesh& tmesh,
                 VertexPointMap vpm,
                 const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor                       halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference                            Point_ref;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS                                         Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor, Box_policy> Box;

  std::vector<Box> boxes;
  boxes.reserve(cc_border_hedges.size());

  for(halfedge_descriptor h : cc_border_hedges)
  {
    const Point_ref p = get(vpm, source(h, tmesh));
    const Point_ref q = get(vpm, target(h, tmesh));
    CGAL_assertion(!gt.equal_3_object()(p, q));

    boxes.emplace_back(p.bbox() + q.bbox(), h);
  }

  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(boxes.size());

  for(Box& b : boxes)
    box_ptr.push_back(&b);

  typedef boost::function_output_iterator<CGAL::internal::Throw_at_output>          Throwing_output_iterator;
  typedef internal::Strict_intersect_edges<Box, TriangleMesh, VertexPointMap,
                                           GeomTraits, Throwing_output_iterator>    Throwing_filter;
  Throwing_filter throwing_filter(tmesh, vpm, gt, Throwing_output_iterator());

  try
  {
    const std::ptrdiff_t cutoff = 2000;
    CGAL::box_self_intersection_d<Parallel_if_available_tag>(box_ptr.begin(), box_ptr.end(), throwing_filter, cutoff);
  }
  catch(CGAL::internal::Throw_at_output_exception&)
  {
    return false;
  }

  return true;
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename PolyhedralEnvelope>
bool remove_self_intersections_with_hole_filling(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                                 std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                                 std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                                 TriangleMesh& tmesh,
                                                 bool local_self_intersection_removal,
                                                 const double strong_dihedral_angle,
                                                 const double weak_dihedral_angle,
                                                 const PolyhedralEnvelope& cc_envelope,
                                                 VertexPointMap vpm,
                                                 const GeomTraits& gt)
{
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/zone_border.polylines.txt");
  out << std::setprecision(17);
  for(const auto& h : cc_border_hedges)
    out << "2 " << tmesh.point(source(h, tmesh)) << " " << tmesh.point(target(h, tmesh)) << std::endl;
  out.close();
#endif

  if(!is_simple_3(cc_border_hedges, tmesh, vpm, gt))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Hole filling cannot handle non-simple border" << std::endl;
#endif
    return false;
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING
  // Do not try to impose sharp edge constraints if we are not doing local-only self intersections removal
  local_self_intersection_removal = false;
#endif

  bool success = false;
  if(local_self_intersection_removal)
  {
    success = fill_hole_with_constraints(cc_border_hedges, cc_faces, working_face_range, tmesh,
                                         strong_dihedral_angle, weak_dihedral_angle,
                                         cc_envelope, vpm, gt);
  }
  else
  {
    success = fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh,
                        cc_envelope, vpm, gt);
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(success)
  {
    if(local_self_intersection_removal)
      ++self_intersections_solved_by_constrained_hole_filling;
    else
      ++self_intersections_solved_by_unconstrained_hole_filling;
  }
#endif

  return success;
}

// the parameter `step` controls how many extra layers of faces we take around the range `faces_to_remove`
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename Visitor>
std::pair<bool, bool>
remove_self_intersections_one_step(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& faces_to_remove,
                                   std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                   TriangleMesh& tmesh,
                                   const int step,
                                   const bool preserve_genus,
                                   const bool only_treat_self_intersections_locally,
                                   const double strong_dihedral_angle,
                                   const double weak_dihedral_angle,
                                   const double containment_epsilon,
                                   VertexPointMap vpm,
                                   const GeomTraits& gt,
                                   Visitor& visitor)
{
  typedef boost::graph_traits<TriangleMesh>                               graph_traits;
  typedef typename graph_traits::vertex_descriptor                        vertex_descriptor;
  typedef typename graph_traits::halfedge_descriptor                      halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                          face_descriptor;

  std::set<face_descriptor> faces_to_remove_copy = faces_to_remove;

#if defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG) || defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT)
  static int call_id = -1;
  ++call_id;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "##### running remove_self_intersections_one_step (#" << call_id << "), step " << step
            << " with " << faces_to_remove.size() << " intersecting faces\n";
#endif

  bool something_was_done = false; // indicates if a region was successfully remeshed
  bool all_fixed = true; // indicates if all removal went well
  // indicates if a removal was not possible because the region handle has
  // some boundary cycle of halfedges
  bool topology_issue = false;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: is_valid in one_step(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";
  std::cout.flush();

  unsolved_self_intersections = 0;
#endif

  CGAL_precondition(is_valid_polygon_mesh(tmesh));
#if defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG) || defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT)
  int cc_id = -1;
#endif

  while(!faces_to_remove.empty())
  {
    if (visitor.stop())
      return std::make_pair(false, false);
    visitor.start_component_handling();
    visitor.status_update(faces_to_remove);
#if defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG) || defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT)
    ++cc_id;
#endif
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: --------------- Removal per CC: " << faces_to_remove.size() << " remaining faces to remove (CC: " << cc_id << " - step: " << step << ")\n";
#endif

    // Process a connected component of faces to remove.
    // collect all the faces from the connected component
    std::set<face_descriptor> cc_faces;
    std::vector<face_descriptor> queue(1, *faces_to_remove.begin()); // temporary queue
    cc_faces.insert(queue.back());
    while(!queue.empty())
    {
      face_descriptor top = queue.back();
      queue.pop_back();
      halfedge_descriptor h = halfedge(top, tmesh);
      for(int i=0; i<3; ++i)
      {
        face_descriptor adjacent_face = face(opposite(h, tmesh), tmesh);
        if(adjacent_face!=boost::graph_traits<TriangleMesh>::null_face())
        {
          if(faces_to_remove.count(adjacent_face) != 0 && cc_faces.insert(adjacent_face).second)
            queue.push_back(adjacent_face);
        }

        h = next(h, tmesh);
      }
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in the current CC\n";
    std::cout << "  DEBUG: first face: " << get(vpm, source(halfedge(*(cc_faces.begin()), tmesh), tmesh)) << " "
              << get(vpm, target(halfedge(*(cc_faces.begin()), tmesh), tmesh)) << " "
              << get(vpm, target(next(halfedge(*(cc_faces.begin()), tmesh), tmesh), tmesh)) << "\n";
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    std::string fname = "results/initial_r_"+std::to_string(call_id)+"_cc_" + std::to_string(cc_id)+"_s_"+std::to_string(step)+".off";
    std::cout << "  DEBUG: Writing initial CC #" << cc_id << " in " << fname << std::endl;
    dump_cc(cc_faces, tmesh, fname);

    fname="results/mesh_at_r_"+std::to_string(call_id)+"_cc_"+std::to_string(cc_id)+"_s_"+std::to_string(step)+".off";
    std::cout << "  DEBUG: Writing current mesh in " << fname << std::endl;
    std::ofstream mout(fname);    mout << std::setprecision(17) << tmesh;
    mout.close();
#endif

    // expand the region to be filled
    if(step > 0)
    {
      expand_face_selection(cc_faces, tmesh, step,
                            make_boolean_property_map(cc_faces),
                            Emptyset_iterator());
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    fname="results/expanded_r_"+std::to_string(call_id)+"_cc_"+std::to_string(cc_id)+"_s_"+std::to_string(step)+".off";
    std::cout << "  DEBUG: Writing expanded CC #" << cc_id << " in " << fname << std::endl;
    dump_cc(cc_faces, tmesh, fname);
#endif

    // try to compactify the selection region by also selecting all the faces included
    // in the bounding box of the initial selection
    std::vector<halfedge_descriptor> stack_for_expension;
    Bbox_3 bb;
    for(face_descriptor fd : cc_faces)
    {
      for(halfedge_descriptor h : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
        bb += get(vpm, target(h, tmesh)).bbox();
        face_descriptor nf = face(opposite(h, tmesh), tmesh);
        if(nf != boost::graph_traits<TriangleMesh>::null_face() && cc_faces.count(nf) == 0)
        {
          stack_for_expension.push_back(opposite(h, tmesh));
        }
      }
    }

    while(!stack_for_expension.empty())
    {
      halfedge_descriptor h = stack_for_expension.back();
      stack_for_expension.pop_back();
      if(cc_faces.count(face(h, tmesh)) == 1)
        continue;

      if(do_overlap(bb, get(vpm, target(next(h, tmesh), tmesh)).bbox()))
      {
        cc_faces.insert(face(h, tmesh));
        halfedge_descriptor candidate = opposite(next(h, tmesh), tmesh);
        if(face(candidate, tmesh) != boost::graph_traits<TriangleMesh>::null_face())
          stack_for_expension.push_back(candidate);

        candidate = opposite(prev(h, tmesh), tmesh);
        if(face(candidate, tmesh) != boost::graph_traits<TriangleMesh>::null_face())
          stack_for_expension.push_back(candidate);
      }
    }

    if(only_treat_self_intersections_locally)
    {
      if(!does_self_intersect(cc_faces, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt)))
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: No self-intersection in CC\n";
#endif

        for(const face_descriptor f : cc_faces)
          faces_to_remove.erase(f);

        continue;
      }
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in expanded and compactified CC\n";
#endif

    // remove faces from the set to process
    for(const face_descriptor f : cc_faces)
      faces_to_remove.erase(f);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    fname="results/expanded_compactified_r_"+std::to_string(call_id)+"_cc_"+std::to_string(cc_id)+"_s_"+std::to_string(step)+".off";
    std::cout << "  DEBUG: Writing expanded and compactified CC #" << cc_id << " in " << fname << std::endl;
    dump_cc(cc_faces, tmesh, fname);
#endif

    //Check for non-manifold vertices in the selection and remove them by selecting all incident faces:
    //  extract the set of halfedges that is on the boundary of the holes to be
    //  made. In addition, we make sure no hole to be created contains a vertex
    //  visited more than once along a hole border (pinched surface)
    //  We save the size of boundary_hedges to make sur halfedges added
    //  from non_filled_hole are not removed.
    bool non_manifold_vertex_remaining_in_selection = false;
    do
    {
      bool non_manifold_vertex_removed = false; //here non-manifold is for the 1D polyline
      std::vector<halfedge_descriptor> boundary_hedges;
      for(face_descriptor fh : cc_faces)
      {
        halfedge_descriptor h = halfedge(fh, tmesh);
        for(int i=0; i<3; ++i)
        {
          if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
            boundary_hedges.push_back(h);

          h = next(h, tmesh);
        }
      }

      // detect vertices visited more than once along
      // a hole border. We then remove all faces incident
      // to such a vertex to force the removal of the vertex.
      // Actually even if two holes are sharing a vertex, this
      // vertex will be removed. It is not needed but since
      // we do not yet have one halfedge per hole it is simpler
      // and does not harm
      std::set<vertex_descriptor> border_vertices;
      for(halfedge_descriptor h : boundary_hedges)
      {
        if(!border_vertices.insert(target(h, tmesh)).second)
        {
          bool any_face_added = false;
          for(halfedge_descriptor hh : halfedges_around_target(h, tmesh))
          {
            if(!is_border(hh, tmesh))
            {
              // add the face to the current selection
              any_face_added |= cc_faces.insert(face(hh, tmesh)).second;
              faces_to_remove.erase(face(hh, tmesh));
            }
          }

          if(any_face_added)
            non_manifold_vertex_removed = true;
          else
            non_manifold_vertex_remaining_in_selection = true;
        }
      }

      if(!non_manifold_vertex_removed)
        break;
    }
    while(true);

    if(preserve_genus && non_manifold_vertex_remaining_in_selection)
    {
      topology_issue = true;
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: CC not handled due to the presence at least one non-manifold vertex\n";

      ++unsolved_self_intersections;
#endif

      visitor.end_component_handling();
      continue; // cannot replace a patch containing a nm vertex by a disk
    }

    // before running this function if preserve_genus=false, we duplicated
    // all of them
    CGAL_assertion(!non_manifold_vertex_remaining_in_selection);

    // Collect halfedges on the boundary of the region to be selected
    // (pointing inside the domain to be remeshed)
    std::vector<halfedge_descriptor> cc_border_hedges;
    for(face_descriptor fd : cc_faces)
    {
      halfedge_descriptor h = halfedge(fd, tmesh);
      for(int i=0; i<3; ++i)
      {
        if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
          cc_border_hedges.push_back(h);

        h = next(h, tmesh);
      }
    }

    if(cc_faces.size() == 1) // it is a triangle nothing better can be done
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
     ++unsolved_self_intersections;
#endif
      visitor.end_component_handling();
      continue;
    }

    working_face_range.insert(cc_faces.begin(), cc_faces.end());

    // Now, we have a proper selection that we can work on.

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTION_NO_POLYHEDRAL_ENVELOPE_CHECK
    Polyhedral_envelope<GeomTraits> cc_envelope;
    if (containment_epsilon!=0)
      cc_envelope = Polyhedral_envelope<GeomTraits>(cc_faces, tmesh, containment_epsilon);
#else
    struct Return_true {
      bool is_empty() const { return true; }
      bool operator()(const std::vector<std::vector<typename GeomTraits::Point_3> >&) const { return true; }
      bool operator()(const TriangleMesh&) const { return true; }
    };

    Return_true cc_envelope;
    CGAL_USE(containment_epsilon);
#endif

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING
    // First, try to smooth if we only care about local self-intersections
    // Two different approaches:
    // - First, try to constrain edges that are in the zone to smooth and whose dihedral angle is large,
    //   but not too large (we don't want to constrain edges that are created by foldings)
    // - If that fails, try to smooth without any constraints, but make sure that the deviation from
    //   the first zone is small
    //
    // If smoothing fails, the face patch is restored to its pre-smoothing state.
    //
    // Note that there is no need to update the working range because smoothing doesn`t change
    // the number of faces (and old faces are re-used).
    bool fixed_by_smoothing = false;

    if(only_treat_self_intersections_locally)
    {
      fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, tmesh, true,
                                                                    strong_dihedral_angle,
                                                                    weak_dihedral_angle,
                                                                    cc_envelope,
                                                                    vpm, gt);

      if(!fixed_by_smoothing) // try again, but without constraining sharp edges
      {
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: Could not be solved via smoothing with constraints\n";
 #endif

        fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, tmesh, false,
                                                                      strong_dihedral_angle,
                                                                      weak_dihedral_angle,
                                                                      cc_envelope,
                                                                      vpm, gt);
      }
    }

    if(fixed_by_smoothing)
    {
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Solved with smoothing!\n";
 #endif

      something_was_done = true;
      visitor.end_component_handling();
      continue;
    }
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    else
    {
      std::cout << "  DEBUG: Could not be solved via smoothing\n";
    }
 #endif
#endif // ndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Trying hole-filling based approach\n";
#endif

    int selection_chi = euler_characteristic_of_selection(cc_faces, tmesh);

    if( selection_chi!=1 ) // not a topological disk
    {
      if (preserve_genus)
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: CC not handled, selection is not a topological disk (preserve_genus=true)\n";
        ++unsolved_self_intersections;
#endif
        topology_issue = true;
        continue;
      }

      // check if the topological issue is created by border edges
      std::vector<halfedge_descriptor> mesh_non_border_hedges;
      std::set<halfedge_descriptor> mesh_border_hedge;

      for(halfedge_descriptor h : cc_border_hedges)
      {
        if(!is_border(opposite(h, tmesh), tmesh))
          mesh_non_border_hedges.push_back(h);
        else
          mesh_border_hedge.insert(opposite(h, tmesh));
      }

      if (mesh_border_hedge.empty())
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: CC not handled, selection is not a topological disk (preserve_genus=false)\n";
        ++unsolved_self_intersections;
#endif
        topology_issue = true;
        continue;
      }

      // we look for cycles of border halfedges and update selection_chi
      while(!mesh_border_hedge.empty())
      {
        // we must count the number of cycle of boundary edges
        halfedge_descriptor h_b = *mesh_border_hedge.begin(), h=h_b;
        mesh_border_hedge.erase(mesh_border_hedge.begin());
        do
        {
          h = next(h, tmesh);
          if(h == h_b)
          {
            // found a cycle
            selection_chi += 1;
            break;
          }
          else
          {
            typename std::set<halfedge_descriptor>::iterator it = mesh_border_hedge.find(h);
            if(it == mesh_border_hedge.end())
              break; // not a cycle does not count

            mesh_border_hedge.erase(it);
          }
        }
        while(true);
      }

      if(selection_chi!=1)
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: CC not handled, selection is not a topological disk even if"
                  << " boundary cycles are removed: chi=" << selection_chi << "\n";
      ++unsolved_self_intersections;
#endif
        topology_issue = true;
        visitor.end_component_handling();
        continue;
      }
      else
      {
        cc_border_hedges.swap(mesh_non_border_hedges);

        // count the number of cycles of halfedges of the boundary
        std::map<vertex_descriptor, vertex_descriptor> bhs;
        for(halfedge_descriptor h : cc_border_hedges)
          bhs[source(h, tmesh)] = target(h, tmesh);

        int nbc=0;
        while(!bhs.empty())
        {
          ++nbc;
          std::pair<vertex_descriptor, vertex_descriptor > top=*bhs.begin();
          bhs.erase(bhs.begin());

          do
          {
            typename std::map<vertex_descriptor, vertex_descriptor>::iterator it_find = bhs.find(top.second);
            if(it_find == bhs.end())
              break;

            top = *it_find;
            bhs.erase(it_find);
          }
          while(true);
        }

        if(nbc != 1)
        {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
          std::cout << "  DEBUG: CC not handled because it is not a topological disk("
                    << nbc << " boundary cycles)\n";
          ++unsolved_self_intersections;
#endif

          all_fixed = false;
          visitor.end_component_handling();
          continue;
        }
        else
        {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
          std::cout << "  DEBUG: CC that is not a topological disk but has only one boundary cycle(preserve_genus=false)\n";
#endif
        }
      }
    }

    if(!remove_self_intersections_with_hole_filling(cc_border_hedges, cc_faces, working_face_range,
                                                    tmesh, only_treat_self_intersections_locally,
                                                    strong_dihedral_angle, weak_dihedral_angle,
                                                    cc_envelope,
                                                    vpm, gt))
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Failed to fill hole\n";
      ++unsolved_self_intersections;
#endif

      all_fixed = false;
    }
    else
    {
      something_was_done = true;
    }
    visitor.end_component_handling();
  }

  if(!something_was_done)
  {
    faces_to_remove.swap(faces_to_remove_copy);
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Nothing was changed during this step, self-intersections won`t be recomputed." << std::endl;
#endif
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::stringstream oss;
  oss << "results/after_step_" << step << ".off" << std::ends;
  std::ofstream(oss.str().c_str()) << std::setprecision(17) << tmesh;
#endif

  return std::make_pair(all_fixed, topology_issue);
}

} // namespace internal

namespace experimental {

template <class TriangleMesh>
struct Remove_self_intersection_default_visitor
{
  bool stop() const { return false; }
  template <class FaceContainer>
  void status_update(const FaceContainer&) {}
  void start_main_loop() {}
  void end_main_loop() {}
  void start_iteration() {}
  void end_iteration() {}
  void start_component_handling() {}
  void end_component_handling() {}
  void parameters_used( bool /* parameters_used(preserve_genus */,
                        bool /* only_treat_self_intersections_locally */,
                        int /* max_steps */,
                        double /* strong_dihedral_angle */,
                        double /* weak_dihedral_angle */,
                        double /* containment_epsilon */ ) {}
};

template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_self_intersections(const FaceRange& face_range,
                               TriangleMesh& tmesh,
                               const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef boost::graph_traits<TriangleMesh>                                 graph_traits;
  typedef typename graph_traits::face_descriptor                            face_descriptor;

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       GeomTraits;
  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));

  bool preserve_genus = choose_parameter(get_parameter(np, internal_np::preserve_genus), true);
  const bool only_treat_self_intersections_locally = choose_parameter(get_parameter(np, internal_np::apply_per_connected_component), false);

  // When treating intersections locally, we don't want to grow the working range too much as
  // either the solution is found fast, or it's too difficult and neither local smoothing or local
  // hole filling are going to provide nice results.
  const int default_max_step = only_treat_self_intersections_locally ? 2 : 7;
  const int max_steps = choose_parameter(get_parameter(np, internal_np::number_of_iterations), default_max_step);

  // @fixme give it its own named parameter rather than abusing 'with_dihedral_angle'?
  const double strong_dihedral_angle = choose_parameter(get_parameter(np, internal_np::with_dihedral_angle), 60.);

  // detect_feature_pp NP (unused for now)
  const double weak_dihedral_angle = 0.; // choose_parameter(get_parameter(np, internal_np::weak_dihedral_angle), 20.);

  struct Return_false {
    bool operator()(std::pair<face_descriptor, face_descriptor>) const { return false; }
  };

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::filter_t,
    NamedParameters,
    Return_false//default
  > ::type  Output_iterator_predicate;
  Output_iterator_predicate out_it_predicates
    = choose_parameter<Return_false>(get_parameter(np, internal_np::filter));

  // use containment check
  const double containment_epsilon = choose_parameter(get_parameter(np, internal_np::polyhedral_envelope_epsilon), 0.);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";
  std::cout << "\tpreserve_genus: " << preserve_genus << std::endl;
  std::cout << "\tonly_treat_self_intersections_locally: " << only_treat_self_intersections_locally << std::endl;
  std::cout << "\tmax_steps: " << max_steps << std::endl;
  std::cout << "\tstrong_dihedral_angle: " << strong_dihedral_angle << std::endl;
  std::cout << "\tweak_dihedral_angle: " << weak_dihedral_angle << std::endl;
  std::cout << "\tcontainment_epsilon: " << containment_epsilon << std::endl;
#endif

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Remove_self_intersection_default_visitor<TriangleMesh>//default
  > ::type Visitor;
  Visitor visitor = choose_parameter<Visitor>(get_parameter(np, internal_np::visitor));

  visitor.parameters_used(preserve_genus,
                          only_treat_self_intersections_locally,
                          max_steps,
                          strong_dihedral_angle,
                          weak_dihedral_angle,
                          containment_epsilon);

  if(!preserve_genus)
    duplicate_non_manifold_vertices(tmesh, np);

  // Look for self-intersections in the mesh and remove them
  int step = -1;
  bool all_fixed = true; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> faces_to_remove;
  std::set<face_descriptor> working_face_range(face_range.begin(), face_range.end());

  visitor.start_main_loop();
  while(++step < max_steps)
  {
    if (visitor.stop()) break;
    visitor.start_iteration();

    if(faces_to_remove.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;

      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(working_face_range, tmesh,
                         CGAL::filter_output_iterator(std::back_inserter(self_inter), out_it_predicates));
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: " << self_inter.size() << " intersecting pairs" << std::endl;
#endif
      for(const Face_pair& fp : self_inter)
      {
        faces_to_remove.insert(fp.first);
        faces_to_remove.insert(fp.second);
      }
    }

    if(faces_to_remove.empty() && all_fixed)
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "DEBUG: There are no more faces to remove." << std::endl;
#endif
      break;
    }

    visitor.status_update(faces_to_remove);

    std::tie(all_fixed, topology_issue) =
      internal::remove_self_intersections_one_step(
          faces_to_remove, working_face_range, tmesh,
          step, preserve_genus, only_treat_self_intersections_locally,
          strong_dihedral_angle, weak_dihedral_angle, containment_epsilon, vpm, gt, visitor);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    if(all_fixed && topology_issue)
        std::cout << "DEBUG: boundary cycles of boundary edges involved in self-intersections.\n";
#endif

    visitor.end_iteration();
  }
  visitor.end_main_loop();

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "solved by constrained smoothing: " << internal::self_intersections_solved_by_constrained_smoothing << std::endl;
  std::cout << "solved by unconstrained smoothing: " << internal::self_intersections_solved_by_unconstrained_smoothing << std::endl;
  std::cout << "solved by constrained hole-filling: " << internal::self_intersections_solved_by_constrained_hole_filling << std::endl;
  std::cout << "solved by unconstrained hole-filling: " << internal::self_intersections_solved_by_unconstrained_hole_filling << std::endl;
  std::cout << "unsolved: " << internal::unsolved_self_intersections << std::endl;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream("results/final.off") << std::setprecision(17) << tmesh;
#endif

  return step < max_steps;
}

template <typename FaceRange, typename TriangleMesh>
bool remove_self_intersections(const FaceRange& face_range, TriangleMesh& tmesh)
{
  return remove_self_intersections(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_self_intersections(TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  return remove_self_intersections(faces(tmesh), tmesh, np);
}

template <typename TriangleMesh>
bool remove_self_intersections(TriangleMesh& tmesh)
{
  return remove_self_intersections(tmesh, parameters::all_default());
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H
