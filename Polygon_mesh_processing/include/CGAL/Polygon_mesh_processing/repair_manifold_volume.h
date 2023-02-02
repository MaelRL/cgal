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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Real_timer.h>

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

template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t fill_hole(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                      PolygonMesh& pmesh,
                      const NamedParameters& np = parameters::default_values())
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  std::vector<vertex_descriptor> hole_vertices;
  std::vector<face_descriptor> hole_faces;

  // Use naive hole filling for now to limit failures
  PMP::triangulate_hole(pmesh, h, std::back_inserter(hole_faces));
  if(hole_faces.empty())
  {
    PMP::triangulate_hole(pmesh, h, std::back_inserter(hole_faces),
                          CGAL::parameters::use_delaunay_triangulation(false));
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

  return filled_hole_count;
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
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using GeomTraits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT = typename GeomTraits::FT;

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

    if(PMP::squared_edge_length(e, tmesh) < squared_threshold)
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


} // namespace internal

namespace experimental {

// Note: pipelines already exist for: ESRI, nframes, Addup, ...

// This function is a pipeline of PMP repair functions
//
// The input mesh is a polygon soup with all kinds of defects
//
// The output mesh is a *valid* polygon mesh in the sense that
// it is a watertight manifold surface with no self-intersections
template <typename PointRange, typename PolygonRange,
          typename PolygonMesh,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool repair_manifold_volume(PointRange& points,
                            PolygonRange& polygons,
                            PolygonMesh& pmesh,
                            const NamedParametersIn& np_in = parameters::default_values() ,
                            const NamedParametersOut& np_out = parameters::default_values() )
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

  // Large CCs thresholds
  const FT threshold_value = epsilon;

  // ###################################################################
  // ### Purge useless stuff
  // ###################################################################

  internal::merge_close_points(points, polygons, close_point_merging_threshold, np_in);

  /// @todo Purge useless (unreachable from infinity) faces

  PMP::repair_polygon_soup(points, polygons, np_in);

  /// @todo Remove almost-coplanar sheets

  // ###################################################################
  // ### Start with local stuff | @todo loop until a certain criterion?
  // ###################################################################

  // Orient (this duplicates non-manifoldness occurrences)
  PMP::orient_polygon_soup(points, polygons, np_in);

  PMP::polygon_soup_to_polygon_mesh(points, polygons, pmesh, np_in, np_out);

  PMP::triangulate_faces(pmesh, np_out);

  PMP::remove_degenerate_faces(pmesh, np_out);

  // One round of stitching to simplify snapping operations
  PMP::stitch_borders(pmesh, np_out);
  PMP::merge_reversible_connected_components(pmesh, np_out);

  auto snapping_subpipeline = [&](const bool preprocess)
  {
    bool something_happened = false;
    do
    {
      if(preprocess)
      {
        // Fill trivial holes
        internal::fill_trivial_holes(pmesh, trivial_hole_halfedge_count_threshold,
                                     trivial_hole_border_length_threshold, np_out);

        // Remove small CCs (that did not attach to something bigger during snap+stitch w/o preprocessing)
        PMP::remove_connected_components_of_negligible_size(pmesh, np_out.area_threshold(negligible_cc_area_threshold)
                                                                         .volume_threshold(negligible_cc_volume_threshold));
        PMP::keep_large_connected_components(pmesh, negligible_cc_face_count_threshold); // at least X faces in a CC

        // Remove almost degenerate elements (there should not be any real degenerate elements at this point)
        internal::collapse_small_edges(pmesh, small_edge_threshold, np_out);
        PMP::remove_almost_degenerate_faces(pmesh, np_out);

        // Stitching might have refused to act because of non-manifoldness between CCs
        PMP::stitch_borders(pmesh, np_out);
        PMP::merge_reversible_connected_components(pmesh, np_out);
      }

      // Snap matched border subsets, with increasing tolerance (zipping) within the min/max distances
      std::size_t snapped_n = internal::iterative_snap(pmesh, min_snapping_threshold,
                                                       max_snapping_threshold, np_out);

      something_happened = (snapped_n > 0);

#ifdef CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG
      std::cout << "snapping pipeline (" << preprocess << "), snapped: " << snapped_n << std::endl;
#endif
    }
    while(something_happened);
  };

  // once without too many harsh preprocessing (cleaning) steps, once with
  snapping_subpipeline(false);
  snapping_subpipeline(true);

  // Self-intersections (local, i.e. within CC) | not sure about this being here
  PMP::experimental::remove_self_intersections(pmesh, np_out.apply_per_connected_component(true)
                                                            .preserve_genus(true));

  // ###################################################################
  // ### Attempt global stuff | @todo loop until valid?
  // ###################################################################

  /// @todo Untangle surfaces [how to detect / fix it? + this conflicts with coplanar patch removal]

  // Back to polygon soup for autoref
  // PMP::polygon_mesh_to_polygon_soup(pmesh, points, polygons, np_out); // @fixme this should be using in_np

  // PMP::autorefine(points, polygons, np_in);

  /// @todo Purge useless (unreachable from infinity) faces

  /// @todo Purge surface CCs incident to non-manifold edges if there is a volume
  /// @todo Purge smaller volumes incident non-manifold edges (winding number?)
  // internal::make_manifold_soup();

  // Back to a polygon mesh
  // PMP::polygon_soup_to_polygon_mesh(points, polygons, pmesh, np_in, np_out);

  internal::fill_all_holes(pmesh, np_out);

  /// @todo non-manifoldness
  // PMP::remove_geometric_non_manifold_vertices();

  PMP::orient_to_bound_a_volume(pmesh, np_out);
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLD_VOLUME_H
