// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl
//                 Daniel Zint

#ifndef CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H
#define CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Dual_contouring_internal.h>

#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief creates an indexed face set that represents an isosurface using the Dual Contouring algorithm.
 *
 * @todo summary; citation; link to user manual.
 * @todo Positioning requirements are not clear, need a concept or not be documented
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Parallel_if_available_tag`, `Parallel_tag`, or `Sequential_tag`.
 * \tparam Domain_ must be a model of `IsosurfacingDomainWithGradient`.
 * \tparam PointRange must be a model of the concept `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type can be constructed from the point type of the domain.
 * \tparam PolygonRange must be a model of the concept `RandomAccessContainer` and `BackInsertionSequence`
 *                      whose value type is itself a model of the concepts `RandomAccessContainer`
 *                      and `BackInsertionSequence` whose value type is `std::size_t`.
 * \tparam Positioning is a functor containing the function `position` that takes `domain`, `isovalue`,
 *                     `cell`, and `position` as input and returns a boolean that is `true`
 *                     if the isosurface intersects the cell.
 *
 * \param domain the domain providing input data and its topology
 * \param isovalue value of the isosurface
 * \param points points of the polzgons in the created indexed face set
 * \param polygons each element in the vector describes a polygon using the indices of the points in `points`
 * \param positioning the functor dealing with vertex positioning inside a cell
 */
template <typename Concurrency_tag = Sequential_tag,
          typename Domain_,
          typename PointRange,
          typename PolygonRange,
          typename Positioning = internal::Positioning::QEM_SVD<true> >
void dual_contouring(const Domain_& domain,
                     const typename Domain_::FT isovalue,
                     PointRange& points,
                     PolygonRange& polygons,
                     const Positioning& positioning = Positioning())
{
  // create vertices in each relevant cell
  internal::Dual_contouring_vertex_positioning<Domain_, Positioning> pos_func(domain, isovalue, positioning);
  domain.template iterate_cells<Concurrency_tag>(pos_func);

  // connect vertices around an edge to form a face
  internal::Dual_contouring_face_generation<Domain_> face_generation(domain, isovalue);
  domain.template iterate_edges<Concurrency_tag>(face_generation);

  // copy vertices to point range
  points.resize(pos_func.points_counter);
  for(const auto& vtop : pos_func.map_voxel_to_point)
    points[pos_func.map_voxel_to_point_id[vtop.first]] = vtop.second;

  // copy faces to polygon range
  polygons.reserve(face_generation.faces.size());
  for(const auto& q : face_generation.faces)
  {
    std::vector<std::size_t> vertex_ids;
    for(const auto& v_id : q.second)
    {
      // ignore voxels that are outside the valid region and are not stored in the map
      if(pos_func.map_voxel_to_point_id.count(v_id) > 0)
        vertex_ids.push_back(pos_func.map_voxel_to_point_id[v_id]);
    }

    // ignore degenerated faces
    if(vertex_ids.size() > 2)
      polygons.push_back(vertex_ids);
  }
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H
