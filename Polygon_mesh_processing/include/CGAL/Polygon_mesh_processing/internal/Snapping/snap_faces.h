// Copyright (c) 2018-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_FACES_H
#define CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_FACES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#ifdef CGAL_PMP_SNAP_DEBUG_PP
 #ifndef CGAL_PMP_SNAP_DEBUG
  #define CGAL_PMP_SNAP_DEBUG
 #endif
#endif

#include <CGAL/Polygon_mesh_processing/internal/Snapping/helper.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap_vertices.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap_edges.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Adapted from <CGAL/internal/AABB_tree/AABB_traversal_traits.h>
template <typename TriangleMesh,
          typename VPMS, typename VPMT, typename FacePatchMap,
          typename AABBTraits>
class Face_projection_traits
{
  typedef typename AABBTraits::FT                                          FT;
  typedef typename AABBTraits::Point_3                                     Point;
  typedef typename AABBTraits::Primitive                                   Primitive;
  typedef typename AABBTraits::Bounding_box                                Bounding_box;
  typedef typename AABBTraits::Primitive::Id                               Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id                      Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id                     Object_and_primitive_id;

  typedef CGAL::AABB_node<AABBTraits>                                      Node;

  typedef typename AABBTraits::Geom_traits                                 Geom_traits;
  typedef typename Geom_traits::Vector_3                                   Vector;
  typedef typename Geom_traits::Plane_3                                    Plane;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

public:
  explicit Face_projection_traits(const halfedge_descriptor h,
                                  const std::size_t patch_id,
                                  const FT sq_tolerance,
                                  const TriangleMesh& tm_S,
                                  const VPMS vpm_S,
                                  const TriangleMesh& tm_T,
                                  const FacePatchMap face_patch_map_T,
                                  const VPMT vpm_T,
                                  const AABBTraits& aabb_traits)
    :
      m_h(h),
      m_patch_id(patch_id),
      m_sq_tol(sq_tolerance),
      m_tm_S(tm_S), m_vpm_S(vpm_S),
      m_tm_T(tm_T), m_face_patch_map_T(face_patch_map_T), m_vpm_T(vpm_T),
      m_is_same_mesh((&tm_S == &tm_T)),
      m_continue(true),
      m_closest_point_initialized(false),
      m_traits(aabb_traits),
      m_gt(Geom_traits()) // blame AABB's traits management
  {
    Vector hv(get(vpm_S, source(h, tm_S)), get(vpm_S, target(h, tm_S)));
    m_direction = hv;
    CGAL::Polygon_mesh_processing::internal::normalize(m_direction, m_gt);
  }

  bool go_further() const { return m_continue; }

  void intersection(const Point& query, const Primitive& primitive)
  {
    face_descriptor pf = primitive.id();

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "~~~~ intersection with primitive: " << pf << std::endl;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(pf, m_tm_T)))
      std::cout << " -> V " << target(h, m_tm_T) << " (" << m_tm_T.point(target(h, m_tm_T)) << ")" << std::endl;

      std::cout << "patches: " << get(m_face_patch_map_T, pf) << " " << m_patch_id << std::endl;
#endif

    if(get(m_face_patch_map_T, pf) != m_patch_id)
      return;

    bool is_vertex_of_face = false;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(pf, m_tm_T), m_tm_T))
    {
      if(m_traits.equal_3_object()(get(m_vpm_T, target(h, m_tm_T)), query))
      {
        is_vertex_of_face = true;
        break;
      }
    }

    if(is_vertex_of_face)
    {
      // If we are NOT using the same mesh and the query is (geometrically) equal to one vertex
      // of the target face, we don't want to move the source point away from the target point
      // (because we cannot move the target point).
      if(!m_is_same_mesh)
      {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
        std::cout << "This vertex is stuck because it is equal to a vertex on the target mesh" << std::endl;
#endif
        m_closest_point_initialized = false;
        m_continue = false;
      }

      // skip the primitive if one of its endpoints is the query
      return;
    }

    // We are searching for a point on the target face pf.
    const Point new_closest_point = m_gt.construct_projected_point_3_object()(
      CGAL::internal::Primitive_helper<AABBTraits>::get_datum(primitive, m_traits), query);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "closest point to primitive: " << new_closest_point << std::endl;
#endif

    const FT sqd_to_tentative_closest_pt = CGAL::squared_distance(query, new_closest_point);

    if(!m_closest_point_initialized || sqd_to_tentative_closest_pt < m_sq_dist)
    {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << "is better!" << std::endl;
#endif

      m_closest_point_initialized = true;

      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point;
      m_sq_dist = sqd_to_tentative_closest_pt;
    }
  }

  bool do_intersect(const Point& query, const Node& node) const
  {
    // This should NOT be the anisotropic distance, as we want to test all targets within the tolerance
    return m_traits.compare_distance_object()(query, node.bbox(), m_sq_tol) == CGAL::SMALLER;
  }

  const Point& closest_point() const { return m_closest_point; }
  typename Primitive::Id closest_primitive_id() const { return m_closest_primitive; }
  bool closest_point_initialized() const { return m_closest_point_initialized; }

private:
  halfedge_descriptor m_h;
  const std::size_t m_patch_id;
  const FT m_sq_tol;
  Vector m_direction;

  const TriangleMesh& m_tm_S;
  VPMS m_vpm_S;
  const TriangleMesh& m_tm_T;
  FacePatchMap m_face_patch_map_T;
  VPMT m_vpm_T;
  const bool m_is_same_mesh;

  bool m_continue;
  bool m_closest_point_initialized;
  typename Primitive::Id m_closest_primitive;
  Point m_closest_point;
  FT m_sq_dist;

  const AABBTraits& m_traits;
  Geom_traits m_gt;
};

template <typename ConcurrencyTag>
struct Faces_to_split_map_inserter // Parallel
{
#ifdef CGAL_LINKED_WITH_TBB
  template <typename FacesToSplitMap, typename FaceDescriptor, typename Point>
  void operator()(FacesToSplitMap& Faces_to_split,
                  const FaceDescriptor closest_f,
                  const FaceDescriptor f,
                  const Point& closest_p)
  {
    typename FacesToSplitMap::accessor acc;
    Faces_to_split.insert(acc, closest_f);
    acc->second.emplace_back(f, closest_p);
  }
#endif
};

template <>
struct Faces_to_split_map_inserter<CGAL::Sequential_tag>
{
  template <typename FacesToSplitMap, typename FaceDescriptor, typename Point>
  void operator()(FacesToSplitMap& Faces_to_split,
                  const FaceDescriptor closest_f,
                  const FaceDescriptor f,
                  const Point& closest_p)
  {
    Faces_to_split[closest_f].emplace_back(f, closest_p);
  }
};

// The UniqueVertex is a pair of a container of vertex_descriptor and FT, representing
// vertices with the same geometric position and their associated snapping tolerance
// (defined as the minimum of the tolerances of the vertices of the bunch)
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename VertexWithTolerance, typename TriangleMesh,
          typename FacesToSplitMap, typename AABBTree,
          typename VertexPatchMap_S, typename FacePatchMap_T,
          typename VPMS, typename VPMT, typename GT>
void find_splittable_face(const VertexWithTolerance& vertex_with_tolerance,
                          FacesToSplitMap& faces_to_split,
                          const AABBTree* aabb_tree_ptr,
                          const TriangleMesh& tm_S,
                          VertexPatchMap_S vertex_patch_map_S,
                          VPMS vpm_S,
                          const TriangleMesh& tm_T,
                          FacePatchMap_T face_patch_map_T,
                          VPMT vpm_T,
                          const GT& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor             edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef typename GT::FT                                                         FT;
  typedef typename boost::property_traits<VPMS>::value_type                       Point;
  typedef typename boost::property_traits<VPMS>::reference                        Point_ref;

  typedef typename AABBTree::AABB_traits                                          AABB_traits;

  typedef internal::Face_projection_traits<TriangleMesh, VPMS, VPMT, FacePatchMap_T, AABB_traits> Projection_traits;

  // by construction the whole range has the same position
  const halfedge_descriptor h = vertex_with_tolerance.first;
  const vertex_descriptor v = target(h, tm_S);
  const Point_ref query = get(vpm_S, v);
  const FT sq_tolerance = CGAL::square(vertex_with_tolerance.second);
  const std::size_t patch_id = get(vertex_patch_map_S, v);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "--------------------------- Query: " << v << " (" << query << ")" << std::endl;
#endif

  Projection_traits traversal_traits(h, patch_id, sq_tolerance, tm_S, vpm_S,
                                     tm_T, face_patch_map_T, vpm_T, aabb_tree_ptr->traits());
  aabb_tree_ptr->traversal(query, traversal_traits);

  if(!traversal_traits.closest_point_initialized())
  {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "Couldn't find any single projection point" << std::endl;
#endif
    return;
  }

  const Point& closest_p = traversal_traits.closest_point();

  // The filtering in the AABB tree checks the dist query <-> node bbox, which might be smaller than
  // the actual distance between the query <-> closest point
  face_descriptor closest_f = traversal_traits.closest_primitive_id();
  halfedge_descriptor closest_h = halfedge(closest_f, tm_T);
  bool is_close_enough =
    gt.compare_squared_distance_3_object()(query,
                                           gt.construct_triangle_3_object()(get(vpm_T, source(closest_h, tm_T)),
                                                                            get(vpm_T, target(closest_h, tm_T)),
                                                                            get(vpm_T, target(next(closest_h, tm_T), tm_T))),
                                           sq_tolerance) != LARGER;

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "  Closest point: (" << closest_p << ")" << std::endl
            << "  on edge: " << traversal_traits.closest_primitive_id()
            << "  at sq distance " << gt.compute_squared_distance_3_object()(query, closest_p)
            << " with squared tolerance: " << sq_tolerance
            << " && close enough? " << is_close_enough << std::endl;
#endif

  if(!is_close_enough)
    return;

  // Using a map because we need to know if the same halfedge is split multiple times
  Faces_to_split_map_inserter<ConcurrencyTag>()(faces_to_split, closest_f, vertex_with_tolerance.first, closest_p);
}

#ifdef CGAL_LINKED_WITH_TBB
template <typename PointWithToleranceContainer,
          typename TriangleMesh, typename FacesToSplitMap, typename AABBTree,
          typename VertexPatchMap_S, typename FacePatchMap_T,
          typename VPMS, typename VPMT, typename GT>
struct find_splittable_face_for_parallel_for
{
  find_splittable_face_for_parallel_for(const PointWithToleranceContainer& points_with_tolerance,
                                        FacesToSplitMap& faces_to_split,
                                        const AABBTree* aabb_tree_ptr,
                                        const TriangleMesh& tm_S,
                                        const VertexPatchMap_S vertex_patch_map_S,
                                        const VPMS vpm_S,
                                        const TriangleMesh& tm_T,
                                        const FacePatchMap_T face_patch_map_T,
                                        const VPMT vpm_T,
                                        const GT& gt)
    :
      m_points_with_tolerance(points_with_tolerance),
      m_faces_to_split(faces_to_split), m_aabb_tree_ptr(aabb_tree_ptr),
      m_tm_S(tm_S), m_vertex_patch_map_S(vertex_patch_map_S), m_vpm_S(vpm_S),
      m_tm_T(tm_T), m_face_patch_map_T(face_patch_map_T), m_vpm_T(vpm_T),
      m_gt(gt)
  { }

  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    for(std::size_t i=r.begin(); i!=r.end(); ++i)
    {
      find_splittable_face<CGAL::Parallel_tag>(m_points_with_tolerance[i], m_faces_to_split, m_aabb_tree_ptr,
                                               m_tm_S, m_vertex_patch_map_S, m_vpm_S,
                                               m_tm_T, m_face_patch_map_T, m_vpm_T, m_gt);
    }
  }

private:
  const PointWithToleranceContainer& m_points_with_tolerance;
  FacesToSplitMap& m_faces_to_split;
  const AABBTree* m_aabb_tree_ptr;
  const TriangleMesh& m_tm_S;
  const VertexPatchMap_S m_vertex_patch_map_S;
  const VPMS m_vpm_S;
  const TriangleMesh& m_tm_T;
  const FacePatchMap_T m_face_patch_map_T;
  const VPMT m_vpm_T;
  const GT& m_gt;
};
#endif

// @fixme requires exact constructions and the fact that the projection of the vertices actually fall on the face
template <typename PointRange, typename Graph, typename VPM, typename GeomTraits>
std::size_t split_face_with_multiple_points(const PointRange& points,
                                            const typename boost::graph_traits<Graph>::face_descriptor f,
                                            const Graph& g,
                                            const VPM vpm,
                                            const GeomTraits& gt)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor                    vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor                  halfedge_descriptor;

  typedef typename boost::range_value<PointRange>::type                             Point;
  typedef typename GeomTraits::Vector_3                                             Vector;

  typedef Triangulation_2_projection_traits_3<GeomTraits>                           P_traits;
  typedef CGAL::Triangulation_vertex_base_with_info_2<vertex_descriptor, P_traits>  Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<bool, P_traits>                 Fbb;
  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fbb>                Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                              TDS;
  typedef CGAL::Exact_intersections_tag                                             Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>           CDT;
  typedef typename CDT::Vertex_handle                                               Vertex_handle;
  typedef typename CDT::Edge                                                        Edge;
  typedef typename CDT::Face_handle                                                 Face_handle;

  Vector n = Polygon_mesh_processing::compute_face_normal(f, g, CGAL::parameters::geom_traits(gt));

  P_traits pgt(n);
  CDT cdt(pgt);

  // Border of the face
  Vertex_handle svh;
  for(const halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, g), g))
  {
    if(svh == Vertex_handle())
    {
      svh = cdt.insert(get(vpm, source(h, g)));
      svh->info() = source(h, g);
    }

    Vertex_handle tvh = cdt.insert(get(vpm, target(h, g)));
    tvh->info() = target(h, g);

    cdt.insert_constraint(svh, tvh);

    svh = tvh;
  }

  // @todo add 'do_split' checks

  for(const Point& p : points)
  {
    vertex_descriptor v = add_vertex(g);
    Vertex_handle tvh = cdt.insert(p); // @fixme check validity of the result (new vertex and within the initial face)
    CGAL_assertion(tvh != Vertex_handle());

    tvh->info() = v;
  }

  // Color the CDT2 because the border might not be convex and there could be faces outside
  for(typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(), end=cdt.finite_faces_end(); fit!=end; ++fit)
    (*fit)->info() = false;

  Face_handle inf_fh = cdt.infinite_vertex()->face();
  int inf_i = inf_fh->index(cdt.infinite_vertex());

  bool is_in_interior = inf_fh->is_constrained(inf_i);

  std::stack<Edge> edge_stack;
  Face_handle mfh = inf_fh->neighbor(inf_i);
  edge_stack.emplace_back(mfh, CDT::mirror_index(inf_fh, inf_i));
  std::set<Face_handle> visited_faces; // @todo use tds_data() when available
  while(!edge_stack.empty())
  {
    Edge e = edge_stack.top();
    edge_stack.top();

    Face_handle fh = e.first;
    const int i = e.second;
    if(!visited_faces.insert(fh).second)
      continue;

    fh->info() = is_in_interior;

    if(cdt.is_constrained(CDT::cw(i)))
    {
      if(!is_in_interior)
      {
        // So far we've only been in the exterior, discard useless faces in the stack
        is_in_interior = true;
        edge_stack.clear();
        edge_stack.emplace_back(fh->neighbor(CDT::cw(i)), CDT::mirror_index(fh, CDT::cw(i)));
        continue;
      }
    }
    else
    {
      edge_stack.emplace_back(fh->neighbor(CDT::cw(i)), CDT::mirror_index(fh, CDT::cw(i)));
    }

    if(cdt.is_constrained(CDT::ccw(i)))
    {
      if(!is_in_interior)
      {
        // So far we've only been in the exterior, discard useless faces in the stack
        is_in_interior = true;
        edge_stack.clear();
        edge_stack.emplace_back(fh->neighbor(CDT::ccw(i)), CDT::mirror_index(fh, CDT::ccw(i)));
        continue;
      }
    }
    else
    {
      edge_stack.emplace_back(fh->neighbor(CDT::ccw(i)), CDT::mirror_index(fh, CDT::ccw(i)));
    }
  }

  CGAL::Euler::remove_face(f, g);

  for(typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(), end=cdt.finite_faces_end(); fit!=end; ++fit)
  {
    Face_handle fh = *fit;
    if(!fh->info())
      continue;

    // @fixme is the orientation fine?
    std::array<vertex_descriptor, 3> face_vertices = { fh->vertex(0)->info(),
                                                       fh->vertex(1)->info(),
                                                       fh->vertex(2)->info() };

    CGAL::Euler::add_face(f, g);
  }

  return points.size();
}

template <typename FacesToSplitContainer,
          typename TriangleMesh, typename GeomTraits,
          typename VPMS, typename VPMT>
std::size_t split_faces(FacesToSplitContainer& faces_to_split,
                        TriangleMesh& tm_S,
                        VPMS vpm_S,
                        TriangleMesh& tm_T,
                        VPMT vpm_T,
                        const GeomTraits& gt,
                        const bool is_source_mesh_fixed) // when snapping is B --> A and the mesh B is fixed
{
#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "split " << faces_to_split.size() << " faces" << std::endl;
#endif

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef typename boost::property_traits<VPMS>::value_type                       Point;
  typedef typename boost::property_traits<VPMT>::reference                        Point_ref;
  typedef typename GeomTraits::Vector_3                                           Vector;

  typedef std::pair<halfedge_descriptor, Point>                                   Vertex_with_new_position;
  typedef std::vector<Vertex_with_new_position>                                   Vertices_with_new_position;
  typedef std::pair<const face_descriptor, Vertices_with_new_position>            Face_and_splitters;

  std::size_t snapped_n = 0;

  CGAL::Real_timer timer;
  timer.start();

  for(Face_and_splitters& fs : faces_to_split)
  {
    const face_descriptor f_to_split = fs.first;
    Vertices_with_new_position& splitters = fs.second;

    CGAL_assertion(!splitters.empty());

    if(splitters.size() > 1)
    {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << " _______ Multiple splitting points on the same face" << std::endl;
#endif
      return split_face_with_multiple_points(splitters, f_to_split, tm_T, vpm_T, gt);
    }

    const Vertex_with_new_position& vnp = splitters[0];
    const halfedge_descriptor splitter_h = vnp.first;
    const vertex_descriptor splitter_v = target(splitter_h, tm_S);
    const Point new_position = is_source_mesh_fixed ? get(vpm_S, splitter_v) : vnp.second;

    bool do_split = true;

    // Some splits can create degenerate faces, avoid that
//      if() // @todo
//        do_split = false;

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << " -.-.-. Splitting " << f_to_split << std::endl;
    for(halfedge_descriptor hts : CGAL::halfedges_around_face(halfedge(f_to_split, tm_T)))
      std::cout << " -> V " << target(hts, tm_T) << " (" << tm_T.point(target(hts, tm_T)) << ")" << std::endl;
    std::cout << "With point: " << vnp.second << std::endl;
    std::cout << "Actually split? " << do_split << std::endl;
#endif

    // Split and update positions
    if(!do_split)
      continue;

    halfedge_descriptor res = CGAL::Euler::add_center_vertex(halfedge(f_to_split, tm_T), tm_T);
    vertex_descriptor new_v = target(res, tm_T);

    put(vpm_T, new_v, new_position); // position of the new point on the target mesh
    if(!is_source_mesh_fixed)
      put(vpm_S, splitter_v, new_position);

    ++snapped_n;
  }

  return snapped_n;
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename HalfedgeRange_S, typename HalfedgeRange_T, typename TriangleMesh,
          typename ToleranceMap_S,
          typename VertexPatchMap_S, typename FacePatchMap_T,
          typename NamedParameters_S, typename NamedParameters_T>
std::size_t snap_vertex_face_one_way(const HalfedgeRange_S& halfedge_range_S,
                                     TriangleMesh& tm_S,
                                     ToleranceMap_S tolerance_map_S,
                                     VertexPatchMap_S vertex_patch_map_S,
                                     const HalfedgeRange_T& halfedge_range_T,
                                     TriangleMesh& tm_T,
                                     FacePatchMap_T face_patch_map_T,
                                     const bool is_source_mesh_fixed,
                                     const NamedParameters_S& np_S,
                                     const NamedParameters_T& np_T)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters_T>::type           GT;
  typedef typename GT::FT                                                         FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters_S>::type       VPMS;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters_T>::type       VPMT;

  typedef typename boost::property_traits<VPMT>::value_type                       Point;

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPMT>            Primitive;
  typedef CGAL::AABB_traits<GT, Primitive>                                        AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits>                                            AABB_tree;

  typedef std::pair<halfedge_descriptor, Point>                                   Vertex_with_new_position;
  typedef std::vector<Vertex_with_new_position>                                   Vertices_with_new_position;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPMS vpm_S = choose_parameter(get_parameter(np_S, internal_np::vertex_point), get_property_map(vertex_point, tm_S));
  VPMT vpm_T = choose_parameter(get_parameter(np_T, internal_np::vertex_point), get_property_map(vertex_point, tm_T));
  const GT gt = choose_parameter<GT>(get_parameter(np_S, internal_np::geom_traits));

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Gather unique points in source range..." << std::endl;
#endif

  typedef std::pair<halfedge_descriptor, FT>                                      Vertex_with_tolerance;
  typedef std::vector<Vertex_with_tolerance>                                      Vertices_with_tolerance;

  Vertices_with_tolerance vertices_to_snap;
  vertices_to_snap.reserve(halfedge_range_S.size()); // ensures that iterators stay valid

  // Take the min tolerance for all points that have the same coordinates
  std::map<Point, FT> point_tolerance_map;

  for(halfedge_descriptor h : halfedge_range_S)
  {
//    if(get(locked_vertices_S, target(h, tm_S))) // @todo
//      continue;

    // Skip the source vertex if its two incident halfedges are geometrically identical (it means that
    // the two halfedges are already stitchable and we don't want this common vertex to be used
    // to split a halfedge somewhere else)
    if(get(vpm_S, source(h, tm_S)) == get(vpm_S, target(next(h, tm_S), tm_S)))
      continue;

    const vertex_descriptor v = target(h, tm_S);
    const FT tolerance = get(tolerance_map_S, v);

    vertices_to_snap.emplace_back(h, tolerance);

    std::pair<typename std::map<Point, FT>::iterator, bool> is_insert_successful =
      point_tolerance_map.emplace(get(vpm_S, v), tolerance);
    if(!is_insert_successful.second)
      is_insert_successful.first->second = (std::min)(is_insert_successful.first->second, tolerance);

    #ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Non-conformal query: " << v << " (" << get(vpm_S, v) << "), tolerance: " << tolerance << std::endl;
#endif
  }

  for(auto& p : vertices_to_snap)
    p.second = point_tolerance_map[get(vpm_S, target(p.first, tm_S))];

  // Since we're inserting primitives one by one, we can't pass this shared data in the constructor of the tree
  AABB_Traits aabb_traits;
  aabb_traits.set_shared_data(tm_T, vpm_T);
  AABB_tree aabb_tree(aabb_traits);

  for(halfedge_descriptor h : halfedge_range_T)
    aabb_tree.insert(Primitive(face(h, tm_T), tm_T, vpm_T));

  // Now, check which faces of the target range ought to be split by source vertices
#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "Collect faces to split with " << vertices_to_snap.size() << " vertices" << std::endl;
#endif

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << "Parallel find splittable edges!" << std::endl;
#endif

    typedef tbb::concurrent_hash_map<halfedge_descriptor,
                                     Vertices_with_new_position>                 Concurrent_face_to_split_container;
    typedef internal::Find_splittable_face_for_parallel_for<
              Vertices_with_tolerance, TriangleMesh,
              Concurrent_face_to_split_container, AABB_tree,
              VertexPatchMap_S, FacePatchMap_T, VPMS, VPMT, GT>                  Functor;

    CGAL::Real_timer timer;
    timer.start();

    Concurrent_face_to_split_container faces_to_split;
    Functor f(vertices_to_snap, faces_to_split, &aabb_tree,
              tm_S, vertex_patch_map_S, vpm_S, tm_T, face_patch_map_T, vpm_T, gt);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices_to_snap.size()), f);

    std::cout << "Time to gather faces: " << timer.time() << std::endl;

    return split_faces(faces_to_split, tm_S, vpm_S, tm_T, vpm_T, gt, is_source_mesh_fixed);
  }
  else
#endif // CGAL_LINKED_WITH_TBB
  {
#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << "Sequential find splittable faces!" << std::endl;
#endif

    std::map<face_descriptor, Vertices_with_new_position> faces_to_split;
    for(const Vertex_with_tolerance& vt : vertices_to_snap)
    {
      internal::find_splittable_face(vt, faces_to_split, &aabb_tree,
                                     tm_S, vertex_patch_map_S, vpm_S,
                                     tm_T, face_patch_map_T, vpm_T, gt);
    }

    return split_faces(faces_to_split, tm_S, vpm_S, tm_T, vpm_T, gt, is_source_mesh_fixed);
  }
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename HalfedgeRange_A, typename HalfedgeRange_B, typename TriangleMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_vertex_face_two_way(const HalfedgeRange_A& halfedge_range_A,
                                     TriangleMesh& tm_A,
                                     ToleranceMap_A tolerance_map_A,
                                     const HalfedgeRange_B& halfedge_range_B,
                                     TriangleMesh& tm_B,
                                     ToleranceMap_B tolerance_map_B,
                                     const bool is_self_snapping, // == true if range and meshes are equal
                                     const NamedParameters_A& np_A,
                                     const NamedParameters_B& np_B)
{
#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout.precision(17);
  std::cerr.precision(17);
#endif

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef CGAL::dynamic_vertex_property_t<bool>                                   Vertex_bool_tag;
  typedef CGAL::dynamic_vertex_property_t<std::size_t>                            Vertex_size_t_tag;
  typedef CGAL::dynamic_halfedge_property_t<bool>                                 Halfedge_bool_tag;

  typedef typename boost::property_map<TriangleMesh, Vertex_bool_tag>::type       Locked_vertices;
  typedef typename boost::property_map<TriangleMesh, Halfedge_bool_tag>::type     Locked_halfedges;
  typedef typename boost::property_map<TriangleMesh, Vertex_size_t_tag>::type     Vertex_patch;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::face_patch_t, NamedParameters_A,
      Constant_property_map<face_descriptor, std::size_t> /*default*/ >::type     Face_patch_map_A;
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::face_patch_t, NamedParameters_B,
      Constant_property_map<face_descriptor, std::size_t> /*default*/ >::type     Face_patch_map_B;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::is_default_parameter;
  using CGAL::parameters::get_parameter;

  const bool is_same_mesh = (&tm_A == &tm_B);
  const bool is_second_mesh_fixed = choose_parameter(get_parameter(np_B, internal_np::do_lock_mesh), false);

  // vertex-vertex and vertex-edge snapping is only considered within compatible patches
  Face_patch_map_A face_patch_map_A = choose_parameter(get_parameter(np_A, internal_np::face_patch),
                                                       Constant_property_map<face_descriptor, std::size_t>(-1));
  Face_patch_map_B face_patch_map_B = choose_parameter(get_parameter(np_B, internal_np::face_patch),
                                                       Constant_property_map<face_descriptor, std::size_t>(-1));

  Vertex_patch vertex_patch_map_A = get(Vertex_size_t_tag(), tm_A);
  for(const halfedge_descriptor h : halfedge_range_A)
  {
    halfedge_descriptor h_opp = opposite(h, tm_A);
    for(const vertex_descriptor v : vertices_around_face(h_opp, tm_A))
      put(vertex_patch_map_A, v, get(face_patch_map_A, face(h_opp, tm_A)));
  }

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Non-conformal snapping... Range sizes: "
            << std::distance(halfedge_range_A.begin(), halfedge_range_A.end()) << " and "
            << std::distance(halfedge_range_B.begin(), halfedge_range_B.end()) << std::endl;
#endif

  CGAL_expensive_precondition(is_valid_polygon_mesh(tm_A) && is_triangle_mesh(tm_A));
  CGAL_expensive_precondition(is_valid_polygon_mesh(tm_B) && is_triangle_mesh(tm_B));

  // Steps:

  // - #1 two-way vertex-edge snapping
  // - #2 two-way vertex-face snapping

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// #1 (Two-way vertex-edge snapping)
  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::size_t snapped_n = 0;

  // @fixme include all halfedges of the facerange... ?

  Vertex_patch vertex_patch_map_B = get(Vertex_size_t_tag(), tm_B);
  for(halfedge_descriptor h_B : halfedge_range_B)
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_face(h_B, tm_B))
      put(vertex_patch_map_B, target(h, tm_B), get(face_patch_map_B, face(h, tm_B)));
  }

  snapped_n += internal::snap_vertex_edge_two_way<ConcurrencyTag>(
                 halfedge_range_A, tm_A, tolerance_map_A,
                 halfedge_range_B, tm_B, tolerance_map_B,
                 is_second_mesh_fixed, np_A, np_B);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// #2 (Two one-way vertex-face snapping)
  //////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << " ///////////// Two one-way vertex-face snapping (A --> B) " << std::endl;
#endif

  // @fixme this modifies the ranges...
  snapped_n += internal::snap_vertex_face_one_way<ConcurrencyTag>(
                 halfedge_range_A, tm_A, tolerance_map_A, vertex_patch_map_A,
                 halfedge_range_B, tm_B, face_patch_map_B,
                 false /*source is never fixed*/, np_A, np_B);

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
  std::ofstream("results/vertex_edge_A.off") << std::setprecision(17) << tm_A;
  std::ofstream("results/vertex_edge_B.off") << std::setprecision(17) << tm_B;
#endif

  if(!is_self_snapping)
  {
#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << " ///////////// Two one-way vertex-face snapping (B --> A) " << std::endl;
#endif

    snapped_n += internal::snap_vertex_face_one_way<ConcurrencyTag>(
                   halfedge_range_B, tm_B, tolerance_map_B, vertex_patch_map_B,
                   halfedge_range_A, tm_A, face_patch_map_A,
                   is_second_mesh_fixed, np_B, np_A);
  }

  return snapped_n;
}

} // namespace internal

namespace experimental {

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_borders_on_faces(TriangleMesh& tm_A,
                                  ToleranceMap_A tolerance_map_A,
                                  TriangleMesh& tm_B,
                                  ToleranceMap_B tolerance_map_B,
                                  const NamedParameters_A& np_A,
                                  const NamedParameters_B& np_B)
{
  // @fixme halfedges(A) should be halfedge(vertices(A)) &&&& halfedges->halfedge(faces(B))
  return internal::snap_vertex_face_two_way<ConcurrencyTag>(halfedges(tm_A), tm_A, tolerance_map_A,
                                                            halfedges(tm_B), tm_B, tolerance_map_B,
                                                            false /*is_self_snapping*/, np_A, np_B);
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename ToleranceMap_A, typename ToleranceMap_B>
std::size_t snap_borders_on_faces(TriangleMesh& tm_A,
                                  ToleranceMap_A tolerance_map_A,
                                  TriangleMesh& tm_B,
                                  ToleranceMap_B tolerance_map_B)
{
  return snap_borders_on_faces(tm_A, tolerance_map_A, tm_B, tolerance_map_B,
                               CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

} // end namespace experimental
} // end namespace Polygon_mesh_processing
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_FACES_H
