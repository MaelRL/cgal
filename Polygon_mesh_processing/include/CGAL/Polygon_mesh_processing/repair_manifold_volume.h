  // To check:

  // Merge epsilon-close points
  // ROCK S., WOZNY M. J.: Generating topological information from a  ̈bucket of facets
  //
  // -- matching borders
  // PATEL P. S., MARCUM D. L., R EMOTIGUE M. G Stitching and filling: Creating conformal faceted geometry.
  // BAREQUET G., SHARIR M.: Filling gaps in the boundary of a polyhedron
  //
  // -- hole filling
  // ZHAO W., GAO S., LIN H.: A robust hole-filling algorithm for triangular mesh.
  // LÉVY B.: Dual domain extrapolation
  // PERNOT J. P., MORARU G., VERON P.: Filling holes in meshes using a mechanical model to simulate
  //                                    the curvature variation minimization
  // WANG J., O LIVEIRA M. M.: Filling holes on locally smooth surfaces reconstructed from point clouds
  // T EKUMALLA L. S., C OHEN E.: A hole-filling algorithm for triangular meshes.
  //
  // [Meh...] PODOLAK J., R USINKIEWICZ S.: Atomic volumes for mesh completion.
  //
  // EL-SANA J., VARSHNEY A.: Controlled simplification of genus for polygonal models
  //
  // [CDT3 gap filling] HÉTROY F., REY S., ANDUJAR C., B RUNET P., VINACUA A.: Mesh repair with user-friendly topology control
  //
  // BISCHOFF S., PAVIC D., KOBBELT L.: Automatic restoration of polygon models
  // BISCHOFF S., KOBBELT L.: Structure preserving CAD model repair.





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
#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLD_VOLUME_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLD_VOLUME_H

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/cluster_point_set.h>
#include <CGAL/Constrained_Delaunay_triangulation_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <unordered_set>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// @fixme factorize this with other versions that exist within CGAL
//
// std::vector to property map
// Not using Pointer_property_map because the underlying range is empty initially and can change
template <typename T>
struct Vector_property_map
{
  using Range = std::vector<T>;

  using key_type = std::size_t;
  using value_type = T;
  using reference = value_type&;
  using category = boost::read_write_property_map_tag;

  Vector_property_map() : m_range_ptr(std::make_shared<Range>()) { }

  inline friend void put(const Vector_property_map& map, const key_type& k, const value_type& v)
  {
    CGAL_precondition(map.m_range_ptr != nullptr);

    if(k >= map.m_range_ptr->size())
      map.m_range_ptr->resize(k+1);

    map.m_range_ptr->operator[](k) = v;
  }

  inline friend reference get(const Vector_property_map& map, const key_type& k)
  {
    CGAL_precondition(map.m_range_ptr != nullptr);
    return map.m_range_ptr->operator[](k);
  }

  Range& range() { return *m_range_ptr; }
  const Range& range() const { return *m_range_ptr; }

private:
  std::shared_ptr<Range> m_range_ptr;
};

// go through all the points and merge those are closer than epsilon
// returns the number of merged points
template <typename PointRange, typename PolygonRange, typename FT,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t merge_close_points(PointRange& points,
                               PolygonRange& polygons,
                               const FT cluster_radius,
                               const NamedParameters& np = parameters::default_values())
{
  using Point_3 = typename boost::range_value<PointRange>::type;

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Merge close points with a cluster radius of " << cluster_radius << std::endl;
#endif

  std::unordered_map<Point_3, std::size_t> cluster_ids;
  auto cluster_ids_pmap = boost::make_assoc_property_map(cluster_ids);

  std::size_t nb_clusters
    = CGAL::cluster_point_set(points, cluster_ids_pmap, np.neighbor_radius(cluster_radius));

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << nb_clusters << " clusters" << std::endl;
#endif

  if(nb_clusters == points.size())
    return 0;

  // update the point positions to the average of the clusters
  std::vector<std::size_t> cluster_sizes(nb_clusters);
  for(std::size_t i=0; i<points.size(); ++i)
    ++cluster_sizes[cluster_ids[points[i]]];

  std::vector<std::array<FT, 3> > new_points_coords(nb_clusters, {{0, 0, 0}});
  for(std::size_t i=0; i<points.size(); ++i)
  {
    const std::size_t cluster_id = cluster_ids[points[i]];
    for(std::size_t j=0; j<3; ++j)
      new_points_coords[cluster_id][j] += points[i][j];
  }

  std::vector<Point_3> new_points;
  new_points.reserve(nb_clusters);
  for(std::size_t i=0; i<nb_clusters; ++i)
  {
    new_points.emplace_back(new_points_coords[i][0] / cluster_sizes[i],
                            new_points_coords[i][1] / cluster_sizes[i],
                            new_points_coords[i][2] / cluster_sizes[i]);
  }

  // Update the point ids in the faces
  for(std::size_t i=0; i<polygons.size(); ++i)
  {
    for(std::size_t j=0; j<polygons[i].size(); ++j)
      polygons[i][j] = cluster_ids[points[polygons[i][j]]];
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << points.size() << " points before clustering" << std::endl;
  std::cout << new_points.size() << " points after clustering" << std::endl;
#endif

  // replace the old points by the new ones
  points.swap(new_points);

  return new_points.size() - points.size();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PointRange, typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t cull_interior_polygons(PointRange& points,
                                   PolygonRange& polygons,
                                   const NamedParameters& = parameters::default_values())
{
  // Do not orient

  try
  {
    // Add faces in a CDT3 and flood fill the faces from the exterior
    // everything that has only incident interior cells can be thrown away

  }
  catch(...) // Intersection_of_constraints_exception
  {
    // Identify volumes

    // Remove all faces that are entirely contained within volumes

  }

  // Remove points that are now useless
  return remove_isolated_points_in_polygon_soup(points, polygons);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t fill_hole(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                      PolygonMesh& pmesh,
                      const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
#endif
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  std::vector<vertex_descriptor> hole_vertices;
  std::vector<face_descriptor> hole_faces;

  // Use naive hole filling for now to limit failures
  PMP::triangulate_hole(pmesh, h, np.face_output_iterator(std::back_inserter(hole_faces)));
  if(hole_faces.empty())
  {
    PMP::triangulate_hole(pmesh, h, np.face_output_iterator(std::back_inserter(hole_faces))
                                      .use_delaunay_triangulation(false));
    if(hole_faces.empty())
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
      std::cerr << "Warning: unfillable hole" << std::endl;

      halfedge_descriptor done = h;
      do
      {
        std::cout << pmesh.point(source(h, pmesh)) << std::endl;
        h = next(h, pmesh);
      }
      while(h != done);
#endif
    }
  }

  return hole_faces.size();
}

// returns the number of trivial holes filled
template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t fill_trivial_holes(PolygonMesh& pmesh,
                               const std::size_t border_halfedge_count_threshold,
                               const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT border_length_threshold,
                               const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using FT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT;

  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

  PMP::merge_duplicated_vertices_in_boundary_cycles(pmesh, np);

  std::vector<halfedge_descriptor> border_cycles;
  CGAL::Polygon_mesh_processing::extract_boundary_cycles(pmesh, std::back_inserter(border_cycles));

  std::size_t filled_hole_count = 0;
  for(halfedge_descriptor h : border_cycles)
  {
    std::size_t border_halfedge_count = 0;
    FT border_length = FT(0);
    bool reject = false;

    halfedge_descriptor done = h;
    do
    {
      ++border_halfedge_count;
      if(border_halfedge_count > border_halfedge_count_threshold)
      {
        reject = true;
        break;
      }

      border_length += edge_length(h, pmesh, np);
      exact(border_length); // to avoid issues with EPECK kernel, see PMP::face_border_length()
      if(border_length > border_length_threshold)
      {
        reject = true;
        break;
      }

      h = next(h, pmesh);
    }
    while(h != done);

    if(reject)
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
      std::cout << "Ignoring border incident to h" << h<< std::endl;
#endif
      continue;
    }

    fill_hole(h, pmesh, np);
    ++filled_hole_count;
  }

  return filled_hole_count;
}

template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t fill_all_holes(PolygonMesh& pmesh,
                           const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

  PMP::merge_duplicated_vertices_in_boundary_cycles(pmesh, np);

  std::vector<halfedge_descriptor> border_cycles;
  CGAL::Polygon_mesh_processing::extract_boundary_cycles(pmesh, std::back_inserter(border_cycles));

  std::size_t filled_hole_count = 0;
  for(halfedge_descriptor h : border_cycles)
  {
    fill_hole(h, pmesh, np);
    ++filled_hole_count;
  }

  CGAL_postcondition(is_closed(pmesh));

  return filled_hole_count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// Collapse short edges on the boundary, regardless of the shape of the incident face
template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
void collapse_small_edges(TriangleMesh& tmesh,
                          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT& threshold,
                          const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using FT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT;

  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  const FT squared_threshold = CGAL::square(threshold);
  std::size_t removed = 0;

  std::unordered_set<edge_descriptor> edges_to_remove;
  std::unordered_set<face_descriptor> faces_to_remove;
  for(edge_descriptor e : edges(tmesh))
  {
    if(!is_border(e, tmesh))
      continue;

    if(PMP::squared_edge_length(e, tmesh, np) < squared_threshold)
    {
      if(CGAL::Euler::does_satisfy_link_condition(e, tmesh))
      {
        edges_to_remove.insert(e);
        ++removed;
      }
      else
      {
        halfedge_descriptor h = halfedge(e, tmesh);
        if(!is_border(h, tmesh))
          h = opposite(h, tmesh);

        if(PMP::internal::border_size(h, tmesh) == 3) // lone face
          faces_to_remove.insert(face(opposite(h, tmesh), tmesh)); // can't remove yet since we're looping edges
#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
        else
          std::cerr << "Warning: collapse of small edge " << e << " is not allowed" << std::endl;
#endif
      }
    }
  }

  for(face_descriptor f : faces_to_remove)
  {
    CGAL::Euler::remove_face(halfedge(f, tmesh), tmesh);
    removed += 3;
  }

  for(edge_descriptor e : edges_to_remove)
  {
    CGAL::Euler::collapse_edge(e, tmesh);
    ++removed;
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Removed " << removed << " edges" << std::endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t single_border_snap(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                               TriangleMesh& tmesh,
                               typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT tolerance,
                               const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using FT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Snapping within boundary starting at: " << edge(h, tmesh) << std::endl;
#endif

  // @fixme surface mesh
  using VertexFTMap = typename TriangleMesh::template Property_map<vertex_descriptor, FT>;
  VertexFTMap tol_map = tmesh.template add_property_map<vertex_descriptor, FT>("v:tol", 0).first;

  std::vector<halfedge_descriptor> border;
  halfedge_descriptor done = h;
  do
  {
    border.push_back(h);
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
    std::cout << "border edge: " << tmesh.point(source(h, tmesh)) << " " << tmesh.point(target(h, tmesh)) << std::endl;
    std::cout << "\t dist: " << PMP::edge_length(h, tmesh, np) << std::endl;
#endif
    h = next(h, tmesh);
  }
  while(h != done);

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "border of length: " << border.size() << std::endl;
#endif

  CGAL_assertion(border.size() == PMP::internal::border_size(h, tmesh));
  CGAL_assertion(border.size() >= 3);

  // Uses 'tolerance' at each vertex but bounds it by a % of the smallest incident edge length
  PMP::internal::assign_tolerance_with_local_edge_length_bound(border, tol_map, tolerance, tmesh);

  std::size_t res = PMP::internal::snap_non_conformal(border, tmesh, tol_map,
                                                      border, tmesh, tol_map,
                                                      true /*self snapping*/, np, np);

  PMP::stitch_boundary_cycle(h, tmesh, np);

  return res;
}

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t single_border_snaps(TriangleMesh& tmesh,
                                const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT tolerance,
                                const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Single-border snap" << std::endl;
#endif

  std::vector<halfedge_descriptor> boundaries;
  PMP::extract_boundary_cycles(tmesh, std::back_inserter(boundaries));

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Snapping " << boundaries.size() << " border(s)" << std::endl;
#endif

  std::size_t snapped_n = 0;
  for(halfedge_descriptor h : boundaries)
    snapped_n += single_border_snap(h, tmesh, tolerance, np);

  return snapped_n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// @fixme Horrible complexity currently as it simply groups together faces that have the same the normal
// regardless of actual CCs information in the original mesh
// it should be some kind of box_d intersection with boxes the size of the tolerance;
// for each box pair, check the normals & group
// & cache normals
template <typename TriangleMesh,
          typename FacePatchMap,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t mark_compatible_patches(FacePatchMap fpmap,
                                    const TriangleMesh& tmesh,
                                    const NamedParameters& = parameters::default_values())
{
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using GeomTraits = typename CGAL::GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  std::vector<Vector_3> distinct_patch_normals;

  for(face_descriptor f : faces(tmesh))
    put(fpmap, f, std::size_t(-1));

  std::size_t new_id = 0;
  for(face_descriptor f : faces(tmesh))
  {
    if(get(fpmap, f) != std::size_t(-1)) // already got a patch id, nothing to do
      continue;

    Vector_3 fn = CGAL::Polygon_mesh_processing::compute_face_normal(f, tmesh);
    const FT sq_fn_l = fn.squared_length();

    std::size_t id = new_id;
    for(std::size_t i=0, s=distinct_patch_normals.size(); i<s; ++i)
    {
      const Vector_3& n = distinct_patch_normals[i];
      const FT sq_n_l = n.squared_length();
      const FT sp = CGAL::scalar_product(fn, n);

      // less than 10 degrees difference ==> same patch
//      const FT sq_cos = 0.9698463103929541; // cos²(10°)

      // less than 30 degrees difference ==> same patch
      const FT sq_cos = 0.75; // cos²(30°)

      // less than 60 degrees difference ==> same patch
//      const FT sq_cos = 0.25; // cos²(60°)

      if(CGAL::square(sp) >= sq_fn_l * sq_n_l * sq_cos)
      {
        id = i;
        break;
      }
    }

    if(id == new_id)
    {
      distinct_patch_normals.push_back(fn);
      ++new_id;
    }

    put(fpmap, f, id);
  }

#if 0 //def CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  for(std::size_t i=0; i<new_id; ++i)
  {
    std::vector<face_descriptor> cc_fs;
    for(face_descriptor f : faces(tmesh))
      if(get(fpmap, f) == i)
        cc_fs.push_back(f);

    std::stringstream oss;
    oss << "results/cc_" << i << ".off" << std::ends;
    dump_cc(cc_fs, tmesh, oss.str().c_str());
  }
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << new_id << " different patch(es)" << std::endl;
#endif

  return new_id;
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename FT, typename TriangleMesh, typename NamedParameters>
std::size_t segment_and_snap(TriangleMesh& tmesh,
                             const FT snapping_tolerance,
                             const NamedParameters& np)
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Segment and snap" << std::endl;
#endif

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor          vertex_descriptor;

  typedef CGAL::dynamic_vertex_property_t<FT>                                    Vertex_tolerance_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_tolerance_tag>::type Tolerance_map;

  Tolerance_map tolerance_map = get(Vertex_tolerance_tag(), tmesh);
  for(vertex_descriptor v : vertices(tmesh))
    put(tolerance_map, v, snapping_tolerance);

  typedef CGAL::dynamic_face_property_t<std::size_t>                             Face_patch_id_tag;
  typedef typename boost::property_map<TriangleMesh, Face_patch_id_tag>::type    Face_patch_map;

  Face_patch_map fpmap = get(Face_patch_id_tag(), tmesh);

  // Color the faces of the mesh into compatible components
  mark_compatible_patches(fpmap, tmesh, np);

  std::size_t res = PMP::experimental::snap_borders<ConcurrencyTag>(tmesh, tolerance_map, np.face_patch_map(fpmap));

  PMP::stitch_borders(tmesh, np);
  PMP::merge_reversible_connected_components(tmesh, np);

  PMP::duplicate_non_manifold_vertices(tmesh); // @todo is that required

  return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t all_borders_snap(TriangleMesh& tmesh,
                             const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT& tolerance,
                             const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using FT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "Snapping ALL borders" << std::endl;
#endif

  std::vector<halfedge_descriptor> border_vertices;
  PMP::border_halfedges(tmesh, std::back_inserter(border_vertices));
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << border_vertices.size() << " total border halfedges" << std::endl;
#endif

  using VertexFTMap = typename TriangleMesh::template Property_map<vertex_descriptor, FT>;
  VertexFTMap tol_map = tmesh.template add_property_map<vertex_descriptor, FT>("v:tol", 0).first;
  PMP::internal::assign_tolerance_with_local_edge_length_bound(border_vertices, tol_map, tolerance, tmesh);

  std::size_t res = PMP::internal::snap_non_conformal(border_vertices, tmesh, tol_map,
                                                      border_vertices, tmesh, tol_map,
                                                      true/*is self snapping*/, np, np);

  PMP::remove_degenerate_faces(tmesh, np); // @todo should not be required with https://github.com/CGAL/cgal/pull/6605

  PMP::stitch_borders(tmesh, np);
  PMP::merge_reversible_connected_components(tmesh, np);

  PMP::duplicate_non_manifold_vertices(tmesh); // @todo is that required

  return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t iterative_snap(TriangleMesh& tmesh,
                           typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT tolerance,
                           const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT max_tolerance,
                           const NamedParameters& np = parameters::default_values())
{
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  int snap_id = 0;
#endif

  std::size_t snapped_n = 0;

  while(tolerance <= max_tolerance)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
    std::cout << "Snapping with tolerance " << tolerance << std::endl;
#endif

    // First, try to resolve everything within boundary cycles
    snapped_n += single_border_snaps(tmesh, tolerance, np);
    std::cout << "snapped_n [single] = " << snapped_n << std::endl;

    /// @todo Identify potential border matches (should it return min/max matched border distances?)

    /// @todo Sort matching borders (pair matching to handle things like 2.obj?)

    // Second, try to resolve everything within coplanar-ish close faces
    snapped_n += segment_and_snap(tmesh, tolerance, np);
    std::cout << "snapped_n [segment] = " << snapped_n << std::endl;

    // Third, use all borders
    snapped_n += all_borders_snap(tmesh, tolerance, np);
    std::cout << "snapped_n [all] = " << snapped_n << std::endl;

#if 0//def CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
   std::stringstream oss;
   oss << "results/snap_intermediate_" << snap_id++ << "_" << tolerance << ".off" << std::ends;
   IO::write_polygon_mesh(oss.str().c_str(), tmesh, CGAL::parameters::stream_precision(17));
#endif

    tolerance *= 10;
  }

  return snapped_n;
}

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
void snap(TriangleMesh& tmesh,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT min_snapping_threshold,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT max_snapping_threshold,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT negligible_cc_area_threshold,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT negligible_cc_volume_threshold,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT trivial_hole_halfedge_count_threshold,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT trivial_hole_border_length_threshold,
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT small_edge_threshold,
          const std::size_t negligible_cc_face_count_threshold,
          const NamedParameters& np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  auto snapping_subpipeline = [&](const bool preprocess)
  {
    bool something_happened = false;
    do
    {
      if(preprocess)
      {
        // Fill trivial holes
        internal::fill_trivial_holes(tmesh, trivial_hole_halfedge_count_threshold,
                                     trivial_hole_border_length_threshold, np);

        // Remove small CCs (that did not attach to something bigger during snap+stitch w/o preprocessing)
        PMP::remove_connected_components_of_negligible_size(tmesh, np.area_threshold(negligible_cc_area_threshold)
                                                                     .volume_threshold(negligible_cc_volume_threshold));
        PMP::keep_large_connected_components(tmesh, negligible_cc_face_count_threshold); // at least X faces in a CC

        // Remove almost degenerate elements (there should not be any real degenerate elements at this point)
        internal::collapse_small_edges(tmesh, small_edge_threshold, np);
        PMP::remove_almost_degenerate_faces(tmesh, np);

        // Stitching might have refused to act because of non-manifoldness between CCs
        PMP::stitch_borders(tmesh, np);
        PMP::merge_reversible_connected_components(tmesh, np);
      }

      // Snap matched border subsets, with increasing tolerance (zipping) within the min/max distances
      std::size_t snapped_n = internal::iterative_snap(tmesh, min_snapping_threshold,
                                                       max_snapping_threshold, np);

      something_happened = (snapped_n > 0);

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
      std::cout << "snapping pipeline (" << preprocess << "), snapped: " << snapped_n << std::endl;
#endif
    }
    while(something_happened);
  };

  // once without too many harsh preprocessing (cleaning) steps, and once with
  snapping_subpipeline(false /*preprocess*/);
  snapping_subpipeline(true /*preprocess*/);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PointRange, typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
void clean_non_manifold_incidences(PointRange&,
                                   PolygonRange&,
                                   const NamedParameters& = parameters::default_values())
{
  // @todo Purge surface CCs incident to non-manifold edges if there is a volume

  // @todo Purge small volumes incident non-manifold edges (winding number?)

}

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
void keep_outer_volumes(TriangleMesh&,
                        const NamedParameters& = parameters::default_values())
{
  // remove surfaces

  // keep outer volumes

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PointRange, typename PolygonRange,
          typename CT,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t construct_cdt3(PointRange& points,
                           PolygonRange& polygons,
                           CT& cdt,
                           const NamedParameters& = parameters::default_values())
{
  using Point_3 = typename boost::range_value<PointRange>::type;

  std::string output_filename = "cdt_dump.txt";

  int exit_code = EXIT_SUCCESS;

  auto finally = [&cdt,output_filename]() {
    {
      std::ofstream dump("dump.binary.cgal");
      CGAL::IO::save_binary_file(dump, cdt);
    }
    {
      std::ofstream dump(output_filename);
      dump.precision(17);
      cdt.write_facets(dump, cdt, std::views::filter(cdt.finite_facets(), [&](auto f) {
          return cdt.is_constrained(f);
      }));
    }
    {
      std::ofstream missing_faces("dump_missing_faces.polylines.txt");
      missing_faces.precision(17);
      cdt.recheck_constrained_Delaunay();
      if(cdt.write_missing_subfaces_file(missing_faces)) {
        std::cerr << "ERROR: Missing subfaces!\n";
      }
    }
    {
      std::ofstream missing_edges("dump_missing_segments.polylines.txt");
      missing_edges.precision(17);
      if(cdt.write_missing_segments_file(missing_edges)) {
        std::cerr << "ERROR: Missing segments!\n";
      }
    }
  };

  int poly_id = 0;
  try
  {
    for(auto polygon : polygons)
    {
      std::vector<Point_3> point_polygon;
      for(auto pi : polygon)
        point_polygon.push_back(points[pi]);

      std::cerr << "NEW POLYGON #" << poly_id << '\n';
      const auto coplanar = point_polygon.size() < 3 ||
          std::all_of(point_polygon.begin(), point_polygon.end(),
                      [p1 = point_polygon[0], p2 = point_polygon[1], p3 = point_polygon[2]](auto p) {
                        const auto coplanar =
                            CGAL::orientation(p1, p2, p3, p) == CGAL::COPLANAR;
                        if(!coplanar) {
                          std::cerr << "Non coplanar points: " << p1 << ", " << p2
                                    << ", " << p3 << ", " << p << '\n'
                                    << "  volume: " << volume(p1, p2, p3, p) << '\n';

                        }
                        return coplanar;
                      });
      if(!coplanar)
      {
        std::ofstream out(std::string("dump_noncoplanar_polygon_") + std::to_string(poly_id) + ".off");
        out.precision(17);
        out << "OFF\n" << point_polygon.size() << " 1 0\n";
        for(auto p : point_polygon)
          out << p << '\n';

        out << point_polygon.size() << ' ';
        for(std::size_t i = 0u, end = point_polygon.size(); i < end; ++i)
          out << ' ' << i;

        out << '\n';
        std::cerr << "Polygon is not coplanar\n";
      }

      try
      {
        auto id = cdt.insert_constrained_polygon(point_polygon);
        assert(id == poly_id);
        ++poly_id;
      }
      catch(int error)
      {
        exit_code = error;
      }
      // std::ofstream dump("dump.binary.cgal");
      // CGAL::Mesh_3::save_binary_file(dump, cdt);
    }

    assert(cdt.is_conforming());

    if(exit_code == EXIT_SUCCESS)
    {
      std::cout << "Restore constrained Delaunay" << std::endl;
      try {
        cdt.restore_constrained_Delaunay();
      } catch(int error) {
        exit_code = error;
      }
    }
  }
  catch(CGAL::Failure_exception&)
  {
    finally();
    std::cout << "Failure in CDT3" << std::endl;
    std::exit(1);
  }

  finally();
  assert(cdt.is_conforming());

  std::cout << "CDT: " << cdt.number_of_vertices() << " nv and " << cdt.number_of_cells() << " nc" << std::endl;

  return exit_code;
}

template <typename CT, typename InDomainPmap>
void mark_domain_in_triangulation(CT& ct,
                                  Unique_hash_map<typename CT::Cell_handle, int>& nesting_level,
                                  typename CT::Cell_handle start,
                                  int index,
                                  std::list<typename CT::Facet>& border,
                                  InDomainPmap ipm)
{
  using Cell_handle = typename CT::Cell_handle;
  using Facet = typename CT::Facet;

  if(nesting_level[start] != -1)
    return;

  std::list<Cell_handle> queue;
  queue.push_back(start);

  while(!queue.empty())
  {
    Cell_handle ch = queue.front();
    queue.pop_front();

    if(nesting_level[ch] == -1)
    {
      nesting_level[ch] = index;

      // for T2, it is index % 2 == 1, but here we want only the outside
      // (@todo we could actually just break as soon as we have walked the outer layer)
      if(index == 0)
        put(ipm, ch, false);

      for(int i=0; i<4; i++)
      {
        Facet f(ch,i);
        Cell_handle cn = ch->neighbor(i);
        if(nesting_level[cn] == -1)
        {
          if(ct.is_constrained(f))
            border.push_back(f);
          else
            queue.push_back(cn);
        }
      }
    }
  }
}

template <typename CDT, typename NestingLevel>
void dump_triangulation_faces(const std::string filename,
                              const CDT& cdt,
                              const NestingLevel& nesting_level,
                              bool only_boundary_faces = false)
{
  using Vertex_handle = typename CDT::Vertex_handle;
  using Cell_handle = typename CDT::Cell_handle;

  std::stringstream vertices_ss;
  vertices_ss.precision(17);

  std::stringstream facets_ss;
  facets_ss.precision(17);

  std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
  std::size_t nv = 0;
  std::size_t nf = 0;

  for(auto fit=cdt.finite_facets_begin(), fend=cdt.finite_facets_end(); fit!=fend; ++fit)
  {
    Cell_handle c = fit->first;
    int s = fit->second;

    Cell_handle nc = c->neighbor(s);
    if(only_boundary_faces && (nesting_level[c] == nesting_level[nc]))
      continue;

    std::array<std::size_t, 3> ids;
    for(std::size_t pos=0; pos<3; ++pos)
    {
      Vertex_handle v = c->vertex((s+pos+1)&3);
      auto insertion_res = vertex_to_id.emplace(v, nv);
      if(insertion_res.second)
      {
        vertices_ss << cdt.point(v) << "\n";
        ++nv;
      }

      ids[pos] = insertion_res.first->second;
    }

    facets_ss << "3 " << ids[0] << " " << ids[1] << " " << ids[2] << "\n";
    ++nf;
  }

  std::ofstream out(filename.c_str());
  out << "OFF\n" << nv << " " << nf << " 0\n";
  out << vertices_ss.str() << "\n" << facets_ss.str() << std::endl;
}

template <typename CT, typename InDomainPmap>
void mark_domain_in_triangulation(CT& cdt,
                                  InDomainPmap ipm)
{
  using Cell_handle = typename CT::Cell_handle;
  using Facet = typename CT::Facet;

  Unique_hash_map<Cell_handle, int> nesting_level(-1, cdt.number_of_cells());

  for(Cell_handle c : cdt.all_cell_handles())
    put(ipm, c, true);

  std::list<Facet> border;
  internal::mark_domain_in_triangulation(cdt, nesting_level, cdt.infinite_cell(), 0, border, ipm);
  while(!border.empty())
  {
    Facet f = border.front();
    border.pop_front();

    Cell_handle cn = f.first->neighbor(f.second);
    if(nesting_level[cn] == -1)
      internal::mark_domain_in_triangulation(cdt, nesting_level, cn, nesting_level[f.first]+1, border, ipm);
  }

  dump_triangulation_faces("cdt_boundary.off", cdt, nesting_level, true /*only boundary*/);
}

// same as Alpha Wrap's
template <typename CT, typename InDomainPmap>
bool is_non_manifold(typename CT::Vertex_handle v,
                     const CT& ct,
                     InDomainPmap ipm)
{
  using Cell_handle = typename CT::Cell_handle;

  CGAL_precondition(!ct.is_infinite(v));

  bool is_non_manifold = false;

  std::vector<Cell_handle> inc_cells;
  inc_cells.reserve(64);
  ct.incident_cells(v, std::back_inserter(inc_cells));

  // Flood one inside and outside CC.
  // Process both an inside and an outside CC to also detect edge pinching.
  // If there are still unprocessed afterwards, there is a non-manifoldness issue.
  //
  // Squat the conflict cell data to mark visits
  Cell_handle inside_start = Cell_handle();
  Cell_handle outside_start = Cell_handle();

  for(Cell_handle ic : inc_cells)
  {
    ic->tds_data().clear();
    if(!get(ipm, ic))
      outside_start = ic;
    else if(inside_start == Cell_handle())
      inside_start = ic;
  }

  // fully inside / outside
  if(inside_start == Cell_handle() || outside_start == Cell_handle())
    return false;

  std::stack<Cell_handle> cells_to_visit;
  cells_to_visit.push(inside_start);
  cells_to_visit.push(outside_start);
  while(!cells_to_visit.empty())
  {
    Cell_handle curr_c = cells_to_visit.top();
    cells_to_visit.pop();

    if(curr_c->tds_data().processed())
      continue;
    curr_c->tds_data().mark_processed();

    int v_pos = -1;
    CGAL_assertion_code(bool res = ) curr_c->has_vertex(v, v_pos);
    CGAL_assertion(res);

    for(int j=0; j<4; ++j)
    {
      if(j == v_pos)
        continue;

      Cell_handle neigh_c = curr_c->neighbor(j);
      CGAL_assertion(neigh_c->has_vertex(v));

      if(neigh_c->tds_data().processed() ||
         get(ipm, neigh_c) != get(ipm, curr_c)) // do not cross the boundary
        continue;

      cells_to_visit.push(neigh_c);
    }
  }

  for(Cell_handle ic : inc_cells)
  {
    if(ic->tds_data().is_clear()) // <=> unvisited cell
    {
      is_non_manifold = true;
      break;
    }
  }

  // reset the conflict flags
  for(Cell_handle ic : inc_cells)
    ic->tds_data().clear();

  return is_non_manifold;
}

template <typename K, typename CT, typename InDomainPmap>
void make_manifold_outer_surface(CT& ct,
                                 InDomainPmap ipm)
{
  using FT = typename K::FT;
  using Vertex_handle = typename CT::Vertex_handle;
  using Cell_handle = typename CT::Cell_handle;

  std::stack<Vertex_handle> non_manifold_vertices; // @todo sort somehow?
  for(Vertex_handle v : ct.finite_vertex_handles())
  {
    if(is_non_manifold(v, ct, ipm))
      non_manifold_vertices.push(v);
  }

  auto sq_longest_edge = [&](Cell_handle c) -> FT
  {
    return (std::max)({ squared_distance(ct.point(c, 0), ct.point(c, 1)),
                        squared_distance(ct.point(c, 0), ct.point(c, 2)),
                        squared_distance(ct.point(c, 0), ct.point(c, 3)),
                        squared_distance(ct.point(c, 1), ct.point(c, 2)),
                        squared_distance(ct.point(c, 1), ct.point(c, 3)),
                        squared_distance(ct.point(c, 2), ct.point(c, 3)) });
  };

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << non_manifold_vertices.size() << " initial NMV" << std::endl;
#endif

  while(!non_manifold_vertices.empty())
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
    std::cout << non_manifold_vertices.size() << " NMV in queue" << std::endl;
#endif

    Vertex_handle v = non_manifold_vertices.top();
    non_manifold_vertices.pop();

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
    std::cout << "·";
#endif

    if(!is_non_manifold(v, ct, ipm))
      continue;

    // Prioritize:
    // - cells without bbox vertices
    // - cells that already have a large number of boundary facets
    // - small cells when equal number of boundary facets
    // @todo give topmost priority to cells with > 1 non-manifold vertex?
    auto comparer = [&](Cell_handle l, Cell_handle r) -> bool
    {
      return sq_longest_edge(l) < sq_longest_edge(r);
    };

    std::vector<Cell_handle> inc_cells;
    inc_cells.reserve(64);
    ct.finite_incident_cells(v, std::back_inserter(inc_cells));

#define CGAL_PMP_REPAIR_MANIFOLD_VOLUME_USE_BRUTE_FORCE_MUTABLE_PRIORITY_QUEUE
#ifndef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_USE_BRUTE_FORCE_MUTABLE_PRIORITY_QUEUE
    std::sort(inc_cells.begin(), inc_cells.end(), comparer); // sort once
#endif

    for(auto cit=inc_cells.begin(), cend=inc_cells.end(); cit!=cend; ++cit)
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_USE_BRUTE_FORCE_MUTABLE_PRIORITY_QUEUE
      // sort at every iteration since the number of boundary facets evolves
      std::sort(cit, cend, comparer);
#endif
      Cell_handle ic = *cit;
      CGAL_assertion(!ct.is_infinite(ic));

      // This is where new material is added
      put(ipm, ic, true);

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
      static int i = 0;
      std::string step_name = "results/steps_manifold/step" + std::to_string(static_cast<int>(i++)) + ".off";
      dump_triangulation_faces(step_name, true /*only_boundary_faces*/);
#endif

      // @speed could update the manifold status while tagging
      if(!is_non_manifold(v, ct, ipm))
        break;
    }

    CGAL_assertion(!is_non_manifold(v, ct, ipm));

    std::vector<Vertex_handle> adj_vertices;
    adj_vertices.reserve(64);
    ct.finite_adjacent_vertices(v, std::back_inserter(adj_vertices));

    for(Vertex_handle nv : adj_vertices)
      if(is_non_manifold(nv, ct, ipm))
        non_manifold_vertices.push(nv);
  }

  CGAL_assertion_code(for(Vertex_handle v : ct.finite_vertex_handles()))
  CGAL_assertion(!is_non_manifold(v, ct, ipm));
}

template <typename CT, typename InDomainPmap,
          typename PointRange, typename PolygonRange>
void extract_surface(const CT& ct,
                     InDomainPmap ipm,
                     PointRange& points,
                     PolygonRange& polygons)
{
  namespace PMP = Polygon_mesh_processing;

  using Vertex_handle = typename CT::Vertex_handle;
  using Facet = typename CT::Facet;
  using Cell_handle = typename CT::Cell_handle;

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
  std::cout << "> Extract CT's outer surface()" << std::endl;
#endif

  // CT's boundary faces to polygon soup
  std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
  std::size_t nv = 0;

  for(auto fit=ct.finite_facets_begin(), fend=ct.finite_facets_end(); fit!=fend; ++fit)
  {
    Facet f = *fit;
    if(get(ipm, f.first))
      f = ct.mirror_facet(f);

    const Cell_handle c = f.first;
    const int s = f.second;
    const Cell_handle nh = c->neighbor(s);
    if(get(ipm, c) == get(ipm, nh))
      continue;

    std::array<std::size_t, 3> ids;
    for(int pos=0; pos<3; ++pos)
    {
      Vertex_handle vh = c->vertex(CT::vertex_triple_index(s, pos));
      auto insertion_res = vertex_to_id.emplace(vh, nv);
      if(insertion_res.second) // successful insertion, never-seen-before vertex
      {
        points.push_back(ct.point(vh));
        ++nv;
      }

      ids[pos] = insertion_res.first->second;
    }

    polygons.emplace_back(std::array<std::size_t, 3>{ids[0], ids[1], ids[2]});
  }
}

template <typename PointRange, typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
void make_mesh_manifold(PointRange& points,
                        PolygonRange& polygons,
                        const NamedParameters& np = parameters::default_values())
{
  using Point_3 = typename boost::range_value<PointRange>::type; // @todo np
  using K = typename CGAL::Kernel_traits<Point_3>::type;
  using Vb = CGAL::Base_with_time_stamp<CGAL::Constrained_Delaunay_triangulation_vertex_base_3<K>>;
  using Cb = CGAL::Constrained_Delaunay_triangulation_cell_base_3<K>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
  using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
  using CDT = CGAL::Constrained_Delaunay_triangulation_3<Delaunay>;
  using Cell_handle = typename CDT::Cell_handle;

  CDT cdt;

  std::size_t exit_code = internal::construct_cdt3(points, polygons, cdt, np);
  std::cout << "CDT3 construction exit code: " << exit_code << std::endl;

  std::unordered_map<Cell_handle, bool> is_in_domain;
  auto ipm = boost::make_assoc_property_map(is_in_domain);
  mark_domain_in_triangulation(cdt, ipm);

  make_manifold_outer_surface<K>(cdt, ipm);

  points.clear();
  polygons.clear();
  extract_surface(cdt, ipm, points, polygons);
}

} // namespace internal

namespace experimental {

// This function is a pipeline of PMP repair functions
//
// The input mesh is a polygon soup with all kinds of defects
//
// The output mesh is a *valid* polygon mesh in the sense that
// it is a watertight manifold surface with no self-intersections

// raw-est pipeline: fill, autoref, cdt3 outside, aw3 manifold
template <typename PointRange, typename PolygonRange,
          typename PolygonMesh,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool repair_manifold_volume(PointRange& points,
                            PolygonRange& polygons,
                            PolygonMesh& pmesh,
                            const NamedParametersIn& np_in = parameters::default_values() ,
                            const NamedParametersOut& np_out = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  using GeomTraits = typename GetGeomTraits<PolygonMesh, NamedParametersIn>::type;
  using FT = typename GeomTraits::FT;

  using Point_map = typename CGAL::GetPointMap<PointRange, NamedParametersIn>::type;
  Point_map pm = choose_parameter<Point_map>(get_parameter(np_in, internal_np::point_map));

  // ###################################################################
  // ### Parameters
  // ###################################################################

  FT bbox_diagonal = FT(0);

  Bbox_3 bb;
  for(const auto& p : points)
    bb += get(pm, p).bbox(); // @fixme conceptless stuff

  bbox_diagonal = CGAL::approximate_sqrt(CGAL::square(bb.xmax() - bb.xmin()) +
                                         CGAL::square(bb.ymax() - bb.ymin()) +
                                         CGAL::square(bb.zmax() - bb.zmin()));

  const FT epsilon = 1e-2 * bbox_diagonal;

  // this one is really just to handle numerical errors in the input
  const FT close_point_merging_threshold = 1e-5 * epsilon;

  // trivial hole filling parameters
  const std::size_t trivial_hole_halfedge_count_threshold = 4;
  const FT trivial_hole_border_length_threshold = epsilon;

  // negligible connected component parameters
  const std::size_t negligible_cc_face_count_threshold = 2;
  const FT negligible_cc_area_threshold = CGAL::square(epsilon);
  const FT negligible_cc_volume_threshold = CGAL::square(epsilon) * epsilon;

  // snapping thresholds
  const FT min_snapping_threshold = 1e-5 * epsilon;
  const FT max_snapping_threshold = epsilon;

  // Small edge threshold
  const FT small_edge_threshold = epsilon;

  // ###################################################################
  // ### Purge useless stuff
  // ###################################################################

  /// @todo initial polygon soup culling, somehow?

  internal::merge_close_points(points, polygons, close_point_merging_threshold, np_in);

  PMP::repair_polygon_soup(points, polygons, np_in);

  PMP::triangulate_polygons(points, polygons, np_in);

  // ###################################################################
  // ### Start with local stuff | @todo loop until a certain criterion?
  // ###################################################################

  // Orient (this duplicates non-manifoldness occurrences)
  PMP::orient_polygon_soup(points, polygons, np_in);

  internal::cull_interior_polygons(points, polygons, np_in);

  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(polygons));
  PMP::polygon_soup_to_polygon_mesh(points, polygons, pmesh, np_in, np_out);

  PMP::remove_degenerate_faces(pmesh, np_out);

  // One round of stitching to simplify snapping operations
  PMP::stitch_borders(pmesh, np_out);
  PMP::merge_reversible_connected_components(pmesh, np_out);

  /// @todo Remove almost-coplanar sheets

  internal::snap(pmesh,
                 min_snapping_threshold, max_snapping_threshold,
                 negligible_cc_area_threshold, negligible_cc_volume_threshold,
                 trivial_hole_halfedge_count_threshold, trivial_hole_border_length_threshold,
                 small_edge_threshold,
                 negligible_cc_face_count_threshold,
                 np_out);

  // Self-intersections (local, i.e. within CC)
  // @fixme not sure about this being here
  PMP::experimental::remove_self_intersections(pmesh, np_out.apply_per_connected_component(true)
                                                            .preserve_genus(true)
                                                            .number_of_iterations(1));

  // ###################################################################
  // ### Attempt global stuff
  // ###################################################################

  // Loop because hole filling / non-manifold repair might create new self-intersections
  while(!is_closed(pmesh) || PMP::does_self_intersect(pmesh, np_out))
  {
    /// @todo Untangle surfaces [how to detect / fix it? This also conflicts with coplanar patch removal]

    // back to a polygon soup for autorefine and non-manifoldness treatment
    points.clear();
    polygons.clear();
    PMP::polygon_mesh_to_polygon_soup(pmesh, points, polygons, np_out); // @fixme this function should use np_in

    PMP::autorefine_triangle_soup(points, polygons, np_in);

    internal::cull_interior_polygons(points, polygons, np_in);

    PMP::orient_polygon_soup(points, polygons, np_in); // @fixme is this needed?

    internal::clean_non_manifold_incidences(points, polygons, np_in);

    // back to a polygon mesh
    clear(pmesh);
    CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(polygons));
    PMP::polygon_soup_to_polygon_mesh(points, polygons, pmesh, np_in); // @fixme this function needs output named parameters

    internal::fill_all_holes(pmesh, np_out);

    /// Unpinch volumes
    internal::make_mesh_manifold(pmesh, np_out);
  }

  PMP::orient_to_bound_a_volume(pmesh, np_out);

  return is_closed(pmesh) && !PMP::does_self_intersect(pmesh);
}

template <typename PointRange, typename PolygonRange,
          typename PolygonMesh,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool repair_manifold_volume_mini(PointRange& points,
                                 PolygonRange& polygons,
                                 PolygonMesh& pmesh,
                                 const NamedParametersIn& np_in = parameters::default_values() ,
                                 const NamedParametersOut& np_out = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  // ###################################################################
  // ### Purge useless stuff
  // ###################################################################

  // @tmp
  // PMP::repair_polygon_soup(points, polygons, np_in);

  PMP::triangulate_polygons(points, polygons, np_in);

  // ###################################################################
  // ### Start with local stuff | @todo loop until a certain criterion?
  // ###################################################################

  // Orient (this duplicates non-manifoldness occurrences)
  PMP::orient_polygon_soup(points, polygons, np_in);

  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(polygons));
  PMP::polygon_soup_to_polygon_mesh(points, polygons, pmesh, np_in, np_out);
  CGAL_assertion(is_closed(pmesh));

  PMP::remove_degenerate_faces(pmesh, np_out);
  CGAL_assertion(is_closed(pmesh));

  // ###################################################################
  // ### Global stuff
  // ###################################################################

  internal::fill_all_holes(pmesh, np_out);
  CGAL_assertion(is_closed(pmesh));

  // back to a polygon soup for autorefine and non-manifoldness treatment
  points.clear();
  polygons.clear();
  PMP::polygon_mesh_to_polygon_soup(pmesh, points, polygons, np_out); // @fixme this function should use np_in

  CGAL::IO::write_OFF("before_autorefine.off", points, polygons, np_in.stream_precision(17));

  PMP::autorefine_triangle_soup(points, polygons, np_in);

  std::cout << "Soup post autoref: " << points.size() << " points and " << polygons.size() << " polygons" << std::endl;
  std::cout << "Self-intersects post autoref: " << PMP::does_triangle_soup_self_intersect(points, polygons, np_in) << std::endl;

  CGAL::IO::write_OFF("after_autorefine.off", points, polygons, np_in.stream_precision(17));

  /// Unpinch volumes
  internal::make_mesh_manifold(points, polygons, np_in);

  // back to a polygon mesh
  clear(pmesh);
  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(polygons));
  PMP::polygon_soup_to_polygon_mesh(points, polygons, pmesh, np_in); // @fixme this function needs output named parameters

  PMP::orient_to_bound_a_volume(pmesh, np_out);

  return is_closed(pmesh) && !PMP::does_self_intersect(pmesh);
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLD_VOLUME_H
