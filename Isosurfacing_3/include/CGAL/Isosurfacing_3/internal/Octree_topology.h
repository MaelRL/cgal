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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_OCTREE_TOPOLOGY_H
#define CGAL_ISOSURFACING_3_INTERNAL_OCTREE_TOPOLOGY_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>

#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename GeomTraits> // @todo should not be necessary
class Octree_topology
{
public:
  using Geom_traits = GeomTraits;
  using Octree = Octree_wrapper<Geom_traits>;
  using Vertex_descriptor = typename Octree::Vertex_handle;
  using Edge_descriptor = typename Octree::Edge_handle;
  using Cell_descriptor = typename Octree::Voxel_handle;

  static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
  static constexpr std::size_t VERTICES_PER_CELL = 8;
  static constexpr std::size_t EDGES_PER_CELL = 12;

  using Vertices_incident_to_edge = std::array<Vertex_descriptor, 2>;
  using Cells_incident_to_edge = std::array<Cell_descriptor, 4>;  // @todo: not always 4
  using Cell_vertices = std::array<Vertex_descriptor, 8>;
  using Cell_edges = std::array<Edge_descriptor, 12>;

public:
  Octree_topology(const Octree& octree)
    : octree(octree)
  { }

  Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const
  {
    return octree.edge_vertices(e);
  }

  Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const
  {
    return octree.edge_voxels(e);
  }

  Cell_vertices cell_vertices(const Cell_descriptor& c) const
  {
    return octree.voxel_vertices(c);
  }

  Cell_edges cell_edges(const Cell_descriptor& c) const
  {
    return octree.voxel_edges(c);
  }

  template <typename Functor>
  void iterate_vertices(Functor& f, Sequential_tag) const
  {
    for(const Vertex_descriptor& v : octree.leaf_vertices())
      f(v);
  }

  template <typename Functor>
  void iterate_edges(Functor& f, Sequential_tag) const
  {
    for(const Edge_descriptor& e : octree.leaf_edges())
      f(e);
  }

  template <typename Functor>
  void iterate_cells(Functor& f, Sequential_tag) const
  {
    for(const Cell_descriptor& v : octree.leaf_voxels())
      f(v);
  }

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Functor>
  void iterate_vertices(Functor& f, Parallel_tag) const
  {
    const auto& vertices = octree.leaf_vertices();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(vertices[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices.size()), iterator);
  }

  template <typename Functor>
  void iterate_edges(Functor& f, Parallel_tag) const
  {
    const auto& edges = octree.leaf_edges();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(edges[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()), iterator);
  }

  template <typename Functor>
  void iterate_cells(Functor& f, Parallel_tag) const
  {
    const auto& cells = octree.leaf_voxels();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(cells[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()), iterator);
  }
#endif // CGAL_LINKED_WITH_TBB

private:
  const Octree& octree;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_OCTREE_TOPOLOGY_H
