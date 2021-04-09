// Copyright (c) 2019-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/internal/Repair/helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/utility.h>

#include <iterator>
#include <fstream>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Combinatorial treatment

namespace internal {

template <typename G>
struct Vertex_collector
{
  typedef typename boost::graph_traits<G>::vertex_descriptor                   vertex_descriptor;
  typedef std::map<vertex_descriptor, std::vector<vertex_descriptor> >         Collection;

  const Collection& collection() const { return _coll; }

  bool has_old_vertex(const vertex_descriptor v) const { return _coll.count(v) != 0; }
  void tag_old_vertex(const vertex_descriptor v)
  {
    CGAL_precondition(!has_old_vertex(v));
    _coll[v];
  }

  void collect_vertices(vertex_descriptor v1, vertex_descriptor v2)
  {
    std::vector<vertex_descriptor>& verts = _coll[v1];
    if(verts.empty())
      verts.push_back(v1);
    verts.push_back(v2);
  }

  template<typename OutputIterator>
  void dump(OutputIterator out)
  {
    typedef std::pair<const vertex_descriptor, std::vector<vertex_descriptor> > Pair_type;
    for(const Pair_type& p : _coll)
      *out++ = p.second;
  }

  void dump(Emptyset_iterator) { }

  Collection _coll;
};

// Replaces the current target vertex of a fan of faces, described by two halfedges, with a new vertex
template <typename PolygonMesh, typename VPM, typename ConstraintMap>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
create_new_vertex_for_sector(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_begin_h,
                             typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_last_h,
                             PolygonMesh& pmesh,
                             const VPM& vpm,
                             const ConstraintMap& cmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  vertex_descriptor old_vd = target(sector_begin_h, pmesh);
  vertex_descriptor new_vd = add_vertex(pmesh);
  put(vpm, new_vd, get(vpm, old_vd));

  put(cmap, new_vd, true);

  set_halfedge(new_vd, sector_begin_h, pmesh);
  halfedge_descriptor h = sector_begin_h;
  do
  {
    set_target(h, new_vd, pmesh);

    if(h == sector_last_h)
      break;
    else
      h = prev(opposite(h, pmesh), pmesh);
  }
  while(h != sector_begin_h); // for safety only
  CGAL_postcondition(h != sector_begin_h);

  return new_vd;
}

template <typename PolygonMesh, typename VPM, typename CMAP>
std::size_t make_umbrella_manifold(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                   PolygonMesh& pmesh,
                                   internal::Vertex_collector<PolygonMesh>& dmap,
                                   VPM vpm,
                                   CMAP cmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  std::size_t nb_new_vertices = 0;

  vertex_descriptor old_v = target(h, pmesh);
  put(cmap, old_v, true); // store the duplicates

  // count the number of borders
  int border_counter = 0;
  halfedge_descriptor ih = h, done = ih, border_h = h;
  do
  {
    if(is_border(ih, pmesh))
    {
      border_h = ih;
      ++border_counter;
    }

    ih = prev(opposite(ih, pmesh), pmesh);
  }
  while(ih != done);

  bool is_non_manifold_within_umbrella = (border_counter > 1);
  if(!is_non_manifold_within_umbrella)
  {
    const bool first_time_meeting_v = !dmap.has_old_vertex(old_v);
    if(first_time_meeting_v)
    {
      // The star is manifold, so if it is the first time we have met that vertex,
      // there is nothing to do, we just keep the same vertex.
      set_halfedge(old_v, h, pmesh); // to ensure halfedge(old_v, pmesh) stays valid
      dmap.tag_old_vertex(old_v); // so that we know we have met old_v already, next time, we'll have to duplicate
    }
    else
    {
      // This is not the canonical star associated to 'v'.
      // Create a new vertex, and move the whole star to that new vertex
      halfedge_descriptor last_h = opposite(next(h, pmesh), pmesh);
      vertex_descriptor new_v = create_new_vertex_for_sector(h, last_h, pmesh, vpm, cmap);
      dmap.collect_vertices(old_v, new_v);
      nb_new_vertices = 1;
    }
  }
  // if there is more than one sector, look at each sector and split them away from the main one
  else
  {
    // the first manifold sector, described by two halfedges
    halfedge_descriptor sector_start_h = border_h;
    CGAL_assertion(is_border(border_h, pmesh));

    bool should_stop = false;
    bool is_main_sector = true;
    do
    {
      CGAL_assertion(is_border(sector_start_h, pmesh));

      // collect the sector and split it away if it must be
      halfedge_descriptor sector_last_h = sector_start_h;
      do
      {
        halfedge_descriptor next_h = prev(opposite(sector_last_h, pmesh), pmesh);

        if(is_border(next_h, pmesh))
          break;

        sector_last_h = next_h;
      }
      while(sector_last_h != sector_start_h);
      CGAL_assertion(!is_border(sector_last_h, pmesh));
      CGAL_assertion(sector_last_h != sector_start_h);

      halfedge_descriptor next_start_h = prev(opposite(sector_last_h, pmesh), pmesh);

      // there are multiple CCs incident to this particular vertex, and we should create a new vertex
      // if it's not the first umbrella around 'old_v' or not the first sector, but only not if it's
      // both the first umbrella and first sector.
      bool must_create_new_vertex = (!is_main_sector || dmap.has_old_vertex(old_v));

      // In any case, we must set up the next pointer correctly
      set_next(sector_start_h, opposite(sector_last_h, pmesh), pmesh);

      if(must_create_new_vertex)
      {
        vertex_descriptor new_v = create_new_vertex_for_sector(sector_start_h, sector_last_h, pmesh, vpm, cmap);
        dmap.collect_vertices(old_v, new_v);
        ++nb_new_vertices;
      }
      else
      {
        // We are in the first sector and first star, ensure that halfedge(old_v, pmesh) stays valid
        set_halfedge(old_v, sector_start_h, pmesh);
      }

      is_main_sector = false;
      sector_start_h = next_start_h;
      should_stop = (sector_start_h == border_h);
    }
    while(!should_stop);
  }

  return nb_new_vertices;
}

} // end namespace internal

/// \ingroup PMP_repairing_grp
/// duplicates all the non-manifold vertices of the input mesh.
///
/// @tparam PolygonMesh a model of `HalfedgeListGraph` and `MutableHalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param pm the surface mesh to be repaired
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
///       The type of this map is model of `ReadWritePropertyMap`.
///       If this parameter is omitted, an internal property map for
///       `CGAL::vertex_point_t` should be available in `PolygonMesh`
///    \cgalParamEnd
///   \cgalParamBegin{vertex_is_constrained_map} a writable property map with `vertex_descriptor`
///     as key and `bool` as `value_type`. `put(pmap, v, true)` will be called for each duplicated
///     vertices, as well as the original non-manifold vertex in the input mesh.
///  \cgalParamEnd
///   \cgalParamBegin{output_iterator} a model of `OutputIterator` with value type
///      `std::vector<vertex_descriptor>`. The first vertex of each vector is a non-manifold vertex
///       of the input mesh, followed by the new vertices that were created to fix this precise
///       non-manifold configuration.
///  \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return the number of vertices created.
template <typename PolygonMesh, typename NamedParameters>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pmesh,
                                            const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename internal_np::Lookup_named_param_def<
                                  internal_np::vertex_is_constrained_t,
                                  NamedParameters,
                                  Static_boolean_property_map<vertex_descriptor, false> // default (no constraint pmap)
                                  >::type                            VerticesMap;
  VerticesMap cmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                      Static_boolean_property_map<vertex_descriptor, false>());

  typedef typename internal_np::Lookup_named_param_def<
                                  internal_np::output_iterator_t,
                                  NamedParameters,
                                  Emptyset_iterator>::type           Output_iterator;
  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator),
                                         Emptyset_iterator());

  std::vector<halfedge_descriptor> non_manifold_umbrellas;
  non_manifold_vertices(pmesh, std::back_inserter(non_manifold_umbrellas));

  internal::Vertex_collector<PolygonMesh> dmap;
  std::size_t nb_new_vertices = 0;
  if(!non_manifold_umbrellas.empty())
  {
    for(halfedge_descriptor h : non_manifold_umbrellas)
      nb_new_vertices += internal::make_umbrella_manifold(h, pmesh, dmap, vpm, cmap);

    dmap.dump(out);
  }

  return nb_new_vertices;
}

template <class PolygonMesh>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pmesh)
{
  return duplicate_non_manifold_vertices(pmesh, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Geometrical treatment

// main:
// @todo handle overlapping geodesic disks
//       --> in the future: identify matching zones and fill hole with screened poisson?
// @todo VPM for heat method / heat method on a FFG
//
// optimization:
// @todo can make maps of Points lighter with a vertex as key, and a custom equal comparing actual points
// @todo small sets

enum NM_TREATMENT
{
  SEPARATE = 0,
  CLIP,
  MERGE,
};

namespace internal {

template <typename PolygonMesh>
struct Work_zone
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef CGAL::dynamic_vertex_property_t<double>                             Distance_tag;
  typedef typename boost::property_map<PolygonMesh, Distance_tag>::type       Ditance_map;

  std::set<face_descriptor> faces;
  std::vector<halfedge_descriptor> border;
  Ditance_map distances;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Treatment of pinched stars

template <typename PolygonMesh>
bool is_pinched_nm_vertex(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                          const PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  int number_of_border_halfedges = 0;

  for(const halfedge_descriptor ih : halfedges_around_target(h, pmesh))
    if(is_border(ih, pmesh))
      ++number_of_border_halfedges;

  return (number_of_border_halfedges > 1);
}

// @todo something nicer?
template <typename PolygonMesh, typename VPM, typename GeomTraits>
void adjust_vertex_position(const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                            const PolygonMesh& pmesh,
                            VPM vpm,
                            const GeomTraits& gt)
{
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  const typename boost::graph_traits<PolygonMesh>::degree_size_type d = degree(v, pmesh);
  CGAL_assertion(d > 1);

  auto& pt = get(vpm, v);

  Vector_3 move{CGAL::NULL_VECTOR};
  for(auto hn : CGAL::halfedges_around_target(v, pmesh))
  {
    if(is_border(hn, pmesh))
      continue;

    const auto& opt1 = get(vpm, source(prev(hn, pmesh), pmesh));
    const auto& opt2 = get(vpm, source(hn, pmesh));
    const auto mpt = gt.construct_midpoint_3_object()(opt1, opt2);

    // - "/2" to get a point within the adjacent face
    // - (d-1) is the number of faces in the sector
    move = move + 0.5 * Vector_3(pt, mpt) / (d - 1) ;
  }

  pt += move;
}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
void treat_pinched_star(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                        NM_TREATMENT treatment,
                        PolygonMesh& pmesh,
                        VPM vpm,
                        const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Pinched star on " << h << " treatment: " << treatment << std::endl;
#endif

  if(treatment == SEPARATE)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_SEPARATE_WITH_SMOOTH
    internal::Vertex_collector<PolygonMesh> dmap;
    Static_boolean_property_map<vertex_descriptor, false> cmap; // @fixme proper cmap?
    internal::make_umbrella_manifold(h, pmesh, dmap, vpm, cmap);

    CGAL_assertion(dmap.collection().size() == 1);

    const std::vector<vertex_descriptor>& new_vs = dmap.collection().begin()->second;
    for(const vertex_descriptor nv : new_vs)
      adjust_vertex_position(nv, pmesh, vpm, gt);
#else
    std::vector<halfedge_descriptor> faces_to_delete;
    for(const halfedge_descriptor ih : halfedges_around_target(h, pmesh))
      if(!is_border(ih, pmesh))
        faces_to_delete.push_back(ih);

    for(halfedge_descriptor ih : faces_to_delete)
       Euler::remove_face(ih, pmesh);
#endif
  }
  else if(treatment == CLIP)
  {
    // @todo
    CGAL_assertion(false);
  }
  else
  {
    // treatment == MERGE

    // Add some additional faces to the star such that the center isn't a manifold vertex anymore
    std::list<halfedge_descriptor> border_halfedges;
    halfedge_descriptor ih = h, done = h;
    do
    {
      if(is_border(ih, pmesh))
        border_halfedges.push_back(ih);
      ih = opposite(next(ih, pmesh), pmesh);
    }
    while(ih != done);

    for(const halfedge_descriptor ih : border_halfedges)
    {
      const halfedge_descriptor pih = prev(ih, pmesh);
      const halfedge_descriptor nih = next(ih, pmesh);
      CGAL_assertion(source(pih, pmesh) != target(nih, pmesh));

      // Don't want to add a degenerate face
      //
      // @todo keep this? Despite this, we could theoretically create self-intersections
      if(gt.collinear_3_object()(get(vpm, target(h, pmesh)),
                                 get(vpm, source(ih, pmesh)),
                                 get(vpm, target(nih, pmesh))))
        continue;

      Euler::add_face_to_border(pih, nih, pmesh);
    }
  }
}

template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool treat_pinched_stars(UmbrellaContainer& umbrellas,
                         NM_TREATMENT treatment,
                         PolygonMesh& pmesh,
                         VPM vpm,
                         const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Treat pinched stars..." << std::endl;
#endif

  bool has_pinched_nm_vertices = false;
  std::set<halfedge_descriptor> treated_umbrellas;

  for(const halfedge_descriptor h : umbrellas)
  {
    if(!is_pinched_nm_vertex(h, pmesh))
      continue;

    has_pinched_nm_vertices = true;
    treat_pinched_star(h, treatment, pmesh, vpm, gt);
    treated_umbrellas.insert(h);
  }

  for(halfedge_descriptor h : treated_umbrellas)
    umbrellas.erase(h);

  return has_pinched_nm_vertices;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merging treatment

// Merging combinatorially changes the star, so we don't want to have two adjacent nm vertices
template <typename NMVM, typename PolygonMesh, typename VPM, typename GeomTraits>
void enforce_non_manifold_vertex_separation(NMVM nm_marks,
                                            PolygonMesh& pmesh,
                                            VPM vpm,
                                            const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor          edge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                    Point;

  std::set<edge_descriptor> edges_to_split;
  for(const edge_descriptor e : edges(pmesh))
  {
    if(get(nm_marks, source(e, pmesh)) && get(nm_marks, target(e, pmesh)))
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "edge that needs to be combinatorially separated: " << e
                << " vertices: " << source(e, pmesh) << " " << target(e, pmesh) << std::endl;
#endif
      edges_to_split.insert(e);
    }
  }

  for(const edge_descriptor e : edges_to_split)
  {
    halfedge_descriptor h = halfedge(e, pmesh);
    if(is_border(h, pmesh))
      h = opposite(h, pmesh);

    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    const Point_ref sp = get(vpm, vs);
    const Point_ref tp = get(vpm, vt);
    const Point mp = gt.construct_midpoint_3_object()(sp, tp);

    halfedge_descriptor new_h = Euler::split_edge_and_incident_faces(h, pmesh);
    put(vpm, target(new_h, pmesh), mp);
  }
}

// note that this also refines the mesh
template <typename PolygonMesh, typename VPM, typename GeomTraits>
void construct_work_zone(Work_zone<PolygonMesh>& wz,
                         const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor zh,
                         const typename GeomTraits::FT radius,
                         PolygonMesh& pmesh,
                         VPM vpm,
                         const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor          edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                    Point;

  typedef typename GeomTraits::FT                                             FT;
  typedef typename GeomTraits::Vector_3                                       Vector;

  typedef CGAL::dynamic_edge_property_t<bool>                                 Considered_tag;
  typedef typename boost::property_map<PolygonMesh, Considered_tag>::type     Considered_edge_map;

  typedef CGAL::dynamic_vertex_property_t<FT>                                 Distance_tag;

  wz.distances = get(Distance_tag(), pmesh);

  vertex_descriptor source_v = target(zh, pmesh);

  CGAL::Heat_method_3::estimate_geodesic_distances(pmesh, wz.distances, source_v); // @fixme vpm?

  // Corefine the mesh with a piecewise(face)-linear approximation of the geodesic circle
  std::set<face_descriptor> faces_to_consider;
  std::vector<halfedge_descriptor> edges_to_split;

  std::stack<halfedge_descriptor> halfedges_to_consider;
  halfedges_to_consider.push(opposite(zh, pmesh));

  Considered_edge_map considered_edges = get(Considered_tag(), pmesh);

  while(!halfedges_to_consider.empty())
  {
    const halfedge_descriptor curr_h = halfedges_to_consider.top();
    halfedges_to_consider.pop();

    const edge_descriptor curr_e = edge(curr_h, pmesh);
    put(considered_edges, curr_e, true);

    const bool is_s_in = (get(wz.distances, source(curr_e, pmesh)) <= radius);
    const bool is_t_in = (get(wz.distances, target(curr_e, pmesh)) <= radius);

    if(is_s_in != is_t_in)
    {
      // We want to keep zh pointing to the nm vertex after the edge split
      edges_to_split.push_back((curr_h == opposite(zh, pmesh)) ? zh : curr_h);
    }

    if(!is_border(curr_h, pmesh))
      faces_to_consider.insert(face(curr_h, pmesh));
    if(!is_border(opposite(curr_h, pmesh), pmesh))
      faces_to_consider.insert(face(opposite(curr_h, pmesh), pmesh));

    const vertex_descriptor curr_v = source(curr_h, pmesh);
    if(get(wz.distances, curr_v) > radius)
      continue;

    for(const halfedge_descriptor adj_h : CGAL::halfedges_around_source(curr_h, pmesh))
    {
      if(get(considered_edges, edge(adj_h, pmesh))) // already visited
        continue;

      // want the source of the new halfedges to be equal to 'source(curr_h, pmesh)'
      halfedges_to_consider.push(opposite(adj_h, pmesh));
    }
  }

  std::cout << edges_to_split.size() << " edges to split" << std::endl;

  // Actual split
  for(halfedge_descriptor h : edges_to_split)
  {
    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    const Point_ref spt = get(vpm, vs);
    const Point_ref tpt = get(vpm, vt);
    Vector tsv = gt.construct_vector_3_object()(tpt, spt);

    const FT dist_at_vs = get(wz.distances, vs);
    const FT dist_at_vt = get(wz.distances, vt);
    if(dist_at_vs == radius || dist_at_vt == radius) // nothing to do
      continue;

    Point new_p;
    if(dist_at_vs < dist_at_vt)
    {
      CGAL_assertion(dist_at_vs < radius && radius <= dist_at_vt);
      const FT lambda = (radius - dist_at_vs) / (dist_at_vt - dist_at_vs);
      new_p = gt.construct_translated_point_3_object()(
                spt, gt.construct_scaled_vector_3_object()(tsv, - lambda));
    }
    else
    {
      CGAL_assertion(dist_at_vt < radius && radius <= dist_at_vs);
      const FT lambda = (radius - dist_at_vt) / (dist_at_vs - dist_at_vt);
      new_p = gt.construct_translated_point_3_object()(
                tpt, gt.construct_scaled_vector_3_object()(tsv, lambda));
    }

    halfedge_descriptor new_h = Euler::split_edge_and_incident_faces(h, pmesh);
    put(vpm, target(new_h, pmesh), new_p);
    put(wz.distances, target(new_h, pmesh), radius);

    // @todo might be simplifiable?
    if(!is_border(h, pmesh))
    {
      faces_to_consider.insert(face(h, pmesh));
      faces_to_consider.insert(face(new_h, pmesh));
    }

    if(!is_border(opposite(h, pmesh), pmesh))
    {
      faces_to_consider.insert(face(opposite(h, pmesh), pmesh));
      faces_to_consider.insert(face(opposite(new_h, pmesh), pmesh));
    }
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << faces_to_consider.size() << " faces to consider" << std::endl;

  static int fi = 0;
  std::stringstream oss;
  oss << "results/faces_to_consider_" << fi++ << ".off" << std::ends;
  dump_cc(faces_to_consider, pmesh, oss.str().c_str());
#endif

  for(const face_descriptor f : faces_to_consider)
  {
    bool is_face_in = true;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      if(get(wz.distances, target(h, pmesh)) > radius)
      {
        is_face_in = false;
        break;
      }
    }

    if(is_face_in)
      wz.faces.insert(f);
  }
}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
Work_zone<PolygonMesh> construct_work_zone(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                           const typename GeomTraits::FT radius,
                                           PolygonMesh& pmesh,
                                           VPM vpm,
                                           const GeomTraits& gt)
{
  // 'set' complexity should be fine, it's not supposed to be a large number of faces
  Work_zone<PolygonMesh> wz;
  construct_work_zone(wz, h, radius, pmesh, vpm, gt);
  CGAL_assertion(!wz.faces.empty());

  // If the selection is not a topological disk, just don't do anything
  // because we won't know how to fill it
  if(!is_selection_a_topological_disk(wz.faces, pmesh))
  {
    std::cerr << "Warning: selection is not a topological disk" << std::endl;
    return Work_zone<PolygonMesh>();
  }

  return wz;
}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
void extract_border_of_work_zone(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                 Work_zone<PolygonMesh>& wz,
                                 const PolygonMesh& pmesh,
                                 const VPM vpm,
                                 const GeomTraits& /*gt*/) // @todo use gts properly
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;

  wz.border.reserve(wz.faces.size());

  CGAL::Face_filtered_graph<PolygonMesh> ffg(pmesh, wz.faces);

  halfedge_descriptor bh = boost::graph_traits<PolygonMesh>::null_halfedge();
  for(halfedge_descriptor h : halfedges(ffg))
  {
    if(is_border(h, ffg))
    {
      bh = h;
      break;
    }
  }

  CGAL_assertion(bh != boost::graph_traits<PolygonMesh>::null_halfedge());

  // Below walks the border of the filtered graph, with a twist: if the star is open, we don't
  // care about the part of that border incident to the vertex:

  /*
      __________________               ______________
     |                  |             |              |
     |                  |             |              |
     |     nm vertex    |  should be  |              |
     |       /  \       |  -------->  |              |
     |      /    \      |             |              |
     |____b/      \c____|             |____b    c____|
  */

  // In other words, a sequence of border halfedges that are border halfedges for pmesh is only
  // acceptable if the nm vertex is not part of it.

  // start by trying to find a halfedge that is not on the border of pmesh to get full border sequences
  const Point_ref nmp = get(vpm, target(h, pmesh));
  halfedge_descriptor nmp_ih = boost::graph_traits<PolygonMesh>::null_halfedge();

  halfedge_descriptor done = bh;
  do
  {
    if(get(vpm, target(bh, pmesh)) == nmp)
      nmp_ih = bh;

    if(!is_border(bh, pmesh))
      break;
    bh = prev(bh, ffg);
  }
  while(bh != done);

  if(is_border(bh, pmesh) && // couldn't find a non-border halfedge: 'ffg' is a full CC of 'pmesh'
     nmp_ih != boost::graph_traits<PolygonMesh>::null_halfedge()) // nm vertex is on the border
  {
    // not really sure what to do here, so let's just remove the two halfedges
    // incident to the vertex and keep the rest

    halfedge_descriptor done = bh;
    {
      if(get(vpm, source(bh, pmesh)) != nmp && get(vpm, target(bh, pmesh)) != nmp)
        wz.border.push_back(opposite(bh, pmesh));

      bh = prev(bh, ffg); // walk backwards so that opposites are in the correct order
    }
    while(bh != done);
  }
  else // standard case now
  {
    halfedge_descriptor done = bh;
    do
    {
      if(is_border(bh, pmesh))
      {
        // start of a sequence of border (for pmesh) halfedges
        // we only want to add that to the complete border if it does not contain the nm vertex
        std::vector<halfedge_descriptor> seq;
        bool add_to_sequence = true;

        do
        {
          if(!is_border(bh, pmesh)) // end of the sequence
            break;

          std::cout << "seq hd: " << get(vpm, source(bh, pmesh)) << " || " << get(vpm, target(bh, pmesh)) << std::endl;
          std::cout << "nmp: " << nmp << std::endl;
          std::cout << "get(vpm, source(bh, pmesh)) == nmp : " << (get(vpm, source(bh, pmesh)) == nmp) << std::endl;

          if(get(vpm, source(bh, pmesh)) == nmp)
            add_to_sequence = false;
          else if(add_to_sequence)
            seq.push_back(opposite(bh, pmesh));

          bh = prev(bh, ffg); // walk backwards so that opposites are in the correct order
        }
        while(bh != done);

        CGAL_assertion(!is_border(bh, pmesh) ||
                       (bh == done && is_border(prev(bh, ffg), pmesh)));

        std::cout << "add_to_sequence = " << add_to_sequence << std::endl;
        if(add_to_sequence)
          wz.border.insert(wz.border.end(), std::begin(seq), std::end(seq)); // not worth move-ing pointers
      }
      else
      {
        wz.border.push_back(opposite(bh, ffg));
        bh = prev(bh, ffg); // walk backwards so that opposites are in the correct order
      }
    }
    while(bh != done);
  }
}

// compute the faces to remove
template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
std::vector<Work_zone<PolygonMesh> >
construct_work_zones(UmbrellaContainer& umbrellas,
                     const typename GeomTraits::FT radius,
                     PolygonMesh& pmesh,
                     VPM vpm,
                     const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  std::vector<Work_zone<PolygonMesh> > wzs;
  wzs.reserve(umbrellas.size());

  for(const halfedge_descriptor h : umbrellas)
  {
    std::cout << "### treating zone incident to: " << get(vpm, source(h, pmesh)) << " " << get(vpm, target(h, pmesh)) << std::endl;

    Work_zone<PolygonMesh> wz = construct_work_zone(h, radius, pmesh, vpm, gt);
    if(!wz.faces.empty())
      wzs.push_back(wz);

    std::cout << "### extract border of zone incident to: "
              << get(vpm, source(h, pmesh)) << " " << get(vpm, target(h, pmesh)) << std::endl;
    extract_border_of_work_zone(h, wzs.back(), pmesh, vpm, gt);

    std::cout << "border of zone " << wzs.size() - 1 << std::endl;
    for(const halfedge_descriptor bh : wzs.back().border)
      std::cout << get(vpm, source(bh, pmesh)) << " " << get(vpm, target(bh, pmesh)) << std::endl;
  }

  return wzs;
}

template <typename HolePointContainer, typename ThirdPointContainer, typename FaceIndexContainer>
bool fill_polyline_hole(const HolePointContainer& hole_points,
                        const ThirdPointContainer& third_points,
                        FaceIndexContainer& patch)
{
  triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch));

  if(patch.empty())
  {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Failed to fill a hole using Delaunay search space.\n";
#endif

    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch),
                              parameters::use_delaunay_triangulation(false));
#endif // CGAL_HOLE_FILLING_DO_NOT_USE_DT3
    if(patch.empty())
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "Failed to fill a hole using the whole search space.\n";
#endif
      return false;
    }
  }

  return true;
}

// Currently performed with a hack: transform the borders into a single border
// by manually adding a patch between the closest edges
template <typename HalfedgeContainer_A, typename HalfedgeContainer_B,
          typename VPM_A, typename VPM_B,
          typename PolygonMesh,
          typename OutputIterator,
          typename GeomTraits>
bool two_borders_hole_fill(const HalfedgeContainer_A& bhv_A,
                           const PolygonMesh& pmesh_A,
                           const VPM_A vpm_A,
                           const HalfedgeContainer_B& bhv_B,
                           const PolygonMesh& pmesh_B,
                           const VPM_B vpm_B,
                           OutputIterator out,
                           const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename boost::property_traits<VPM_A>::value_type                  Point;
  typedef typename boost::property_traits<VPM_A>::reference                   Point_ref_A;
  typedef typename boost::property_traits<VPM_B>::reference                   Point_ref_B;

  typedef typename HalfedgeContainer_A::const_iterator                        HCit_A;
  typedef typename HalfedgeContainer_B::const_iterator                        HCit_B;

  typedef CGAL::Triple<int, int, int>                                         Face_indices;

  CGAL_precondition(bhv_A.size() >= 3);
  CGAL_precondition(bhv_B.size() >= 3);

  HCit_A canon_A;
  HCit_B canon_B;

  // @todo avoid the O(n^2) complexity (but the borders are likely small, so...)
  // @todo this type of matching might not be the best idea

  double best_score = std::numeric_limits<double>::max();
  for(HCit_A it_A=std::begin(bhv_A), A_end=std::end(bhv_A); it_A!=A_end; ++it_A)
  {
    const Point_ref_A sap = get(vpm_A, source(*it_A, pmesh_A));
    const Point_ref_A tap = get(vpm_A, target(*it_A, pmesh_A));

    for(HCit_B it_B=std::begin(bhv_B), B_end=std::end(bhv_B); it_B!=B_end; ++it_B)
    {
      const Point_ref_B sbp = get(vpm_B, source(*it_B, pmesh_B));
      const Point_ref_B tbp = get(vpm_B, target(*it_B, pmesh_B));

      const double score = gt.compute_squared_distance_3_object()(sap, tbp)
                         + gt.compute_squared_distance_3_object()(tap, sbp);
      if(score < best_score)
      {
        best_score = score;
        canon_A = it_A;
        canon_B = it_B;
      }
    }
  }

  // polyline
  std::vector<Point> hole_points, third_points;

  // _____   _________   ________
  //      | | <------ | |
  //      | | canon_A | |
  //      | |         | |
  //      | | canon_B | |
  //      | | ------> | |
  // -----   --------   -------

  const Point_ref_A sap = get(vpm_A, source(*canon_A, pmesh_A));
  const Point_ref_A tap = get(vpm_A, target(*canon_A, pmesh_A));
  const Point_ref_B sbp = get(vpm_B, source(*canon_B, pmesh_B));
  const Point_ref_B tbp = get(vpm_B, target(*canon_B, pmesh_B));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "best score:" << std::endl;
  std::cout << "sap: " << sap << std::endl;
  std::cout << "tap: " << tap << std::endl;
  std::cout << "sbp: " << sbp << std::endl;
  std::cout << "tbp: " << tbp << std::endl;
#endif

  std::ofstream border_out("results/hole_border.polylines.txt");
  border_out.precision(17);

  // Only a single ID is needed, the rest is 0, hole.size()-1, and quad_id + 1
  std::size_t quad_id = static_cast<std::size_t>(-1);

  // Walk A's border
  HCit_A last_A = std::prev(bhv_A.end());
  HCit_A it_A = (canon_A == last_A) ? std::begin(bhv_A) : std::next(canon_A);
  do
  {
    const halfedge_descriptor h = *it_A;

    if(!hole_points.empty())
      border_out << "2 " << hole_points.back() << " " << get(vpm_A, source(h, pmesh_A)) << std::endl;

    hole_points.push_back(get(vpm_A, source(h, pmesh_A)));
    std::cout << "point on hole: " << hole_points.back() << std::endl;

    const halfedge_descriptor oh = opposite(h, pmesh_A);
    if(is_border(oh, pmesh_A))
      third_points.push_back(construct_artificial_third_point(h, pmesh_A, vpm_A, gt));
    else
      third_points.push_back(get(vpm_A, target(next(oh, pmesh_A), pmesh_A)));
    it_A = (it_A == last_A) ? std::begin(bhv_A) : std::next(it_A);
  }
  while(it_A != canon_A);

  // vertical
  border_out << "2 " << hole_points.back() << " " << sap << std::endl;

  quad_id = hole_points.size();
  hole_points.push_back(sap);
  third_points.push_back(CGAL::midpoint(sbp, tap));

  // Walk B's border
  HCit_B last_B = std::prev(bhv_B.end());
  HCit_B it_B = (canon_B == last_B) ? std::begin(bhv_B) : std::next(canon_B);
  do
  {
    const halfedge_descriptor h = *it_B;

    border_out << "2 " << hole_points.back() << " " << get(vpm_B, source(h, pmesh_B)) << std::endl;
    hole_points.push_back(get(vpm_B, source(h, pmesh_B)));

    const halfedge_descriptor oh = opposite(h, pmesh_B);
    if(is_border(oh, pmesh_A))
      third_points.push_back(construct_artificial_third_point(h, pmesh_B, vpm_B, gt));
    else
      third_points.push_back(get(vpm_B, target(next(opposite(h, pmesh_B), pmesh_B), pmesh_B)));
    it_B = (it_B == last_B) ? std::begin(bhv_B) : std::next(it_B);
  }
  while(it_B != canon_B);

  border_out << "2 " << hole_points.back() << " " << sbp << std::endl;
  hole_points.push_back(sbp);
  third_points.push_back(CGAL::midpoint(tbp, sap));

  CGAL_assertion(hole_points.size() == third_points.size());

  std::vector<Face_indices> patch;
  fill_polyline_hole(hole_points, third_points, patch);

  // add the missing quad
  patch.emplace_back(quad_id, 0, quad_id+1);
  patch.emplace_back(quad_id+1, 0, hole_points.size() - 1);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream pout("results/patch.off");
  pout << std::setprecision(17);
  pout << 3 * patch.size() << " " << patch.size() << " 0\n";

  for(const Face_indices& f : patch)
  {
    pout << hole_points[f.first] << "\n";
    pout << hole_points[f.second] << "\n";
    pout << hole_points[f.third] << "\n";
  }

  for(std::size_t i=0, ps=patch.size(); i<ps; ++i)
    pout << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
#endif

  for(const Face_indices& face : patch)
  {
    *out++ = std::initializer_list<Point>{ hole_points[face.first],
                                           hole_points[face.second],
                                           hole_points[face.third] };
  }

  return true;
}

template <typename Patch,
          typename PolygonMesh,
          typename VPM,
          typename GeomTraits>
bool fix_patch_orientation(Patch& point_patch,
                           const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h1,
                           const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h2,
                           const PolygonMesh& pmesh,
                           const VPM vpm,
                           const GeomTraits& gt)
{
  typedef typename boost::property_traits<VPM>::reference                       Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  bool ok_orientation_1 = false, ok_orientation_2 = false;

  const Point_ref h1sp = get(vpm, source(h1, pmesh));
  const Point_ref h1tp = get(vpm, target(h1, pmesh));
  const Point_ref h2sp = get(vpm, source(h2, pmesh));
  const Point_ref h2tp = get(vpm, target(h2, pmesh));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "h1sp: " << h1sp << std::endl;
  std::cout << "h1tp: " << h1tp << std::endl;
  std::cout << "h2sp: " << h2sp << std::endl;
  std::cout << "h2tp: " << h2tp << std::endl;
#endif

  for(const auto& face : point_patch)
  {
    for(int i=0; i<3; ++i)
    {
      const Point& p1 = face[i];
      const Point& p2 = face[(i+1)%3];

      if(gt.equal_3_object()(p1, h1sp) && gt.equal_3_object()(p2, h1tp))
        ok_orientation_1 = true;
      if(gt.equal_3_object()(p1, h2sp) && gt.equal_3_object()(p2, h2tp))
        ok_orientation_2 = true;
    }
  }

  std::cout << "orientations: " << ok_orientation_1 << " " << ok_orientation_2 << std::endl;

  if(ok_orientation_1 != ok_orientation_2)
    return false;

  if(!ok_orientation_1)
  {
    CGAL_assertion(false); // @tmp
    for(auto& face : point_patch)
      std::swap(face[0], face[1]);
  }

  return true;
}

template <typename WorkZone, typename PolygonMesh, typename VPM, typename GeomTraits>
bool merge_zones(const WorkZone& wz_1,
                 const WorkZone& wz_2,
                 PolygonMesh& pmesh,
                 const VPM vpm,
                 const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename GeomTraits::Point_3                                        Point;

  const std::vector<halfedge_descriptor>& bhv_1 = wz_1.border;
  const std::vector<halfedge_descriptor>& bhv_2 = wz_2.border;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Merge holes (" << wz_1.faces.size() << " and " << wz_2.faces.size() << ")" << std::endl;
  std::cout << "holes of size: " << bhv_1.size() << " && " << bhv_2.size() << std::endl;
#endif

  // make sure that the holes are topological disks
  std::vector<std::vector<Point> > point_patch;
  point_patch.reserve(2 * bhv_1.size());

  bool success = two_borders_hole_fill(bhv_1, pmesh, vpm, bhv_2, pmesh, vpm,
                                       std::back_inserter(point_patch), gt);
  if(!success)
    return false;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << point_patch.size() << " new faces" << std::endl;
  static int merge_id = 0;
  std::stringstream oss;
  oss << "results/tentative_merge_patch_" << merge_id << ".off" << std::ends;
  dump_tentative_hole(point_patch, oss.str().c_str());
#endif

  if(!check_patch_sanity<PolygonMesh>(point_patch))
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "this patch is INSANE!" << std::endl;
#endif
    return false;
  }

  // @todo shouldn't the orientation always be correct by construction?
  // Do this with combinatorics if it turns out it can be not correct
  success = fix_patch_orientation(point_patch, bhv_1.front(), bhv_2.front(), pmesh, vpm, gt);
  if(!success)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Incompatible orientations" << std::endl;
#endif
    return false; // can't find an orientation compatible with both borders
  }

  std::set<face_descriptor> fs(std::begin(wz_1.faces), std::end(wz_1.faces));
  fs.insert(std::begin(wz_2.faces), std::end(wz_2.faces));

  return replace_faces_with_patch(fs, point_patch, pmesh, vpm);
}

/*
 The idea is to walk the polylines incident to the nm vertex that are on the border of the mesh,
 and modify the boundary so that it doesn't pass through that vertex anymore.
 for example:

*/
template <typename PolygonMesh, typename VPM, typename GeomTraits>
bool fill_zone_with_open_border(const Work_zone<PolygonMesh>& wz,
                                PolygonMesh& pmesh,
                                const VPM vpm,
                                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                    Point;
  typedef typename boost::property_traits<VPM>::reference                     Point_ref;

  typedef CGAL::Triple<int, int, int>                                         Face_indices;

  std::vector<Point> hole_points, third_points;
  for(const halfedge_descriptor bh : wz.border)
  {
    const halfedge_descriptor obh = opposite(bh, pmesh);

    hole_points.push_back(get(vpm, source(bh, pmesh)));

    if(is_border(obh, pmesh))
      third_points.push_back(construct_artificial_third_point(bh, pmesh, vpm, gt));
    else
      third_points.push_back(get(vpm, target(next(obh, pmesh), pmesh)));
  }

  // last one
  const halfedge_descriptor lh = wz.border.back();
  hole_points.push_back(get(vpm, target(lh, pmesh)));

  const halfedge_descriptor olh = opposite(lh, pmesh);
  if(is_border(olh, pmesh))
    third_points.push_back(construct_artificial_third_point(lh, pmesh, vpm, gt));
  else
    third_points.push_back(get(vpm, target(next(olh, pmesh), pmesh)));

  CGAL_assertion(hole_points.size() == third_points.size());

  std::vector<Face_indices> patch;
  fill_polyline_hole(hole_points, third_points, patch);

  std::vector<std::vector<Point> > point_patch;
  point_patch.reserve(patch.size());
  for(const Face_indices& face : patch)
  {
    point_patch.emplace_back(std::initializer_list<Point>{ hole_points[face.first],
                                                           hole_points[face.second],
                                                           hole_points[face.third] });
  }

  return replace_faces_with_patch(wz.faces, point_patch, pmesh, vpm);
}

template <typename WorkZoneContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool fill_zones(const WorkZoneContainer& wzs,
                PolygonMesh& pmesh,
                VPM vpm,
                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor        halfedge_descriptor;

  for(std::size_t i=0, n=wzs.size(); i<n; ++i)
  {
    CGAL_assertion(!wzs[i].faces.empty());
    CGAL_assertion(!wzs[i].border.empty());

    if(wzs[i].faces.size() == 1)
    {
      Euler::remove_face(halfedge(*(std::begin(wzs[i].faces)), pmesh), pmesh);
      continue;
    }

    const halfedge_descriptor f = wzs[i].border.front();
    const halfedge_descriptor b = wzs[i].border.back();

    const bool is_loop = source(f, pmesh) == target(b, pmesh);
    if(!is_loop)
    {
      fill_zone_with_open_border(wzs[i], pmesh, vpm, gt);
    }
    else if(!fill_hole(wzs[i].faces, pmesh, vpm, gt))
    {
      std::cerr << "Failed to fill hole of work zone #" << i << std::endl;
      return false;
    }
  }

  return true;
}

template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
void smooth_umbrellas(UmbrellaContainer& umbrellas,
                      PolygonMesh& pmesh,
                      VPM vpm,
                      const GeomTraits& gt)
{
  for(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h : umbrellas)
    adjust_vertex_position(target(h, pmesh), pmesh, vpm, gt);
}

template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool treat_umbrellas(UmbrellaContainer& umbrellas,
                     const NM_TREATMENT treatment,
                     const typename GeomTraits::FT radius,
                     PolygonMesh& pmesh,
                     VPM vpm,
                     const GeomTraits& gt)
{
  // can't merge if there were pinched stars
  const bool can_employ_merge_strategy = !treat_pinched_stars(umbrellas, treatment, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_SEPARATE_WITH_SMOOTH
  // If smoothing, we don't need any particular work zones
  if(treatment == SEPARATE)
  {
    smooth_umbrellas(umbrellas, pmesh, vpm, gt);
    return true;
  }
#endif

  std::vector<Work_zone<PolygonMesh> > wzs = construct_work_zones(umbrellas, radius, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << wzs.size() << " work zones" << std::endl;
  static int i = 0;
  for(const Work_zone<PolygonMesh>& wz : wzs)
  {
    std::stringstream oss;
    oss << "results/zone_" << i++ << ".off" << std::ends;
    dump_cc(wz.faces, pmesh, oss.str().c_str());
  }

  CGAL::write_polygon_mesh("results/post_construction.off", pmesh, CGAL::parameters::stream_precision(17));
#endif

  if(treatment == SEPARATE)
  {
    // Separate treatment here implies no smoothing but cut&fill
    return fill_zones(wzs, pmesh, vpm, gt);
  }

  // Merging is not (yet?) supported for > 2 umbrellas
  if(!can_employ_merge_strategy || wzs.size() != 2)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Merging treatment requested, but configuration makes it impossible" << std::endl;
#endif
    return fill_zones(wzs, pmesh, vpm, gt);
  }
  else
  {
    const bool success = merge_zones(wzs.front(), wzs.back(), pmesh, vpm, gt);
    if(!success)
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "merging failed, falling back to separating strategy" << std::endl;
#endif
      return fill_zones(wzs, pmesh, vpm, gt);
    }

    return true;
  }
}

} // namespace internal

template <typename PolygonMesh, typename NamedParameters>
void repair_non_manifold_vertices(PolygonMesh& pmesh,
                                  const NM_TREATMENT treatment,
                                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type      VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type          Geom_traits;
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  typedef typename Geom_traits::FT                                            FT;
  typedef typename boost::property_traits<VertexPointMap>::value_type         Point;

  const FT radius = 0.3; // @todo automatically computed or np

  // Collect the non-manifold vertices
  typedef CGAL::dynamic_vertex_property_t<bool>                               Mark;
  typedef typename boost::property_map<PolygonMesh, Mark>::type               Marked_vertices;
  Marked_vertices nm_marks = get(Mark(), pmesh);
  for(vertex_descriptor v : vertices(pmesh))
    put(nm_marks, v, false);

  std::unordered_map<Point, std::set<halfedge_descriptor> > nm_points;
  geometrically_non_manifold_vertices(pmesh, nm_points, nm_marks, np);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << nm_points.size() << " case(s) to treat:" << std::endl;
  for(auto& e : nm_points)
    std::cout << "\t" << e.first << std::endl;
#endif

  internal::enforce_non_manifold_vertex_separation(nm_marks, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::write_polygon_mesh("results/enforced_separation.off", pmesh, CGAL::parameters::stream_precision(17));
#endif

  for(auto& e : nm_points)
  {
    // All these umbrellas have the same center
    std::set<halfedge_descriptor>& umbrellas = e.second;
    CGAL_assertion(!umbrellas.empty());

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "NM vertex at pos (" << get(vpm, target(*(std::begin(umbrellas)), pmesh))
              << ") with " << umbrellas.size() << " incident umbrella(s):" << std::endl;
    for(const halfedge_descriptor h : umbrellas)
      std::cout << h << " (" << get(vpm, source(h, pmesh)) << ") (" << get(vpm, target(h, pmesh)) << ")" << std::endl;
#endif

    internal::treat_umbrellas(umbrellas, treatment, radius, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "done with that nm vertex" << std::endl << std::endl;
    CGAL::write_polygon_mesh("results/intermediary.off", pmesh, CGAL::parameters::stream_precision(17));
#endif
  }

  CGAL_postcondition(is_valid_polygon_mesh(pmesh, true));

  CGAL_postcondition_code(nm_points.clear();)
  CGAL_postcondition_code(geometrically_non_manifold_vertices(pmesh, nm_points, nm_marks, np);)
  CGAL_postcondition(nm_points.empty());
}

template <typename PolygonMesh>
void repair_non_manifold_vertices(PolygonMesh& pmesh,
                                  const NM_TREATMENT treatment)
{
  return repair_non_manifold_vertices(pmesh, treatment, parameters::all_default());
}

template <typename PolygonMesh>
void repair_non_manifold_vertices(PolygonMesh& pmesh)
{
  return repair_non_manifold_vertices(pmesh, SEPARATE, parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
