// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
#define CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_topology.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

namespace CGAL {
namespace Isosurfacing {

/*
 * \ingroup PkgIsosurfacing3Ref
 *
 * \cgalModels `IsosurfacingDomainWithGradient`
 *
 * \brief A domain that respesents an octree that discretizes an implicit function.
 *
 * \tparam GeomTraits must be a model of ``.
 * \tparam PointFunction the type of the implicit function. It must implement
 *                       `GeomTraits::FT operator()(const GeomTraits::Point& point) const`.
 * \tparam Gradient_ the type of the gradient functor. It must implement
 *                   `GeomTraits::Vector operator()(const GeomTraits::Point& point) const`.
 */
template <typename GeomTraits,
          typename PointFunction,
          typename Gradient_>
using Implicit_octree_domain =
  internal::Base_domain<GeomTraits,
                        internal::Octree_topology<GeomTraits>,
                        internal::Octree_geometry<GeomTraits>,
                        internal::Implicit_function_with_geometry<GeomTraits,
                                                                  internal::Octree_geometry<GeomTraits>,
                                                                  PointFunction>,
                        Gradient_>;

/*
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief creates a domain from an octree that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function will be evaluated on the octree vertices.
 *
 * \tparam GeomTraits must be a model of ``.
 * \tparam PointFunction the type of the implicit function. It must implement
 *                       `GeomTraits::FT operator()(const GeomTraits::Point& point) const`.
 * \tparam Gradient_ the type of the gradient functor. It must implement
 *                   `GeomTraits::Vector operator()(const GeomTraits::Point& point) const`.
 *
 * \param bbox a bounding box that specifies the size of the functions domain
 * \param spacing the distance between discretized points on the function
 * \param point_function the function with a point as argument
 * \param gradient a function that describes the gradient of the data
 *
 * \return a new `Implicit_cartesian_grid_domain`
 */
template <typename GeomTraits,
          typename PointFunction,
          typename Gradient_ = Zero_gradient>
Implicit_octree_domain<GeomTraits, PointFunction, Gradient_>
create_implicit_octree_domain(const std::shared_ptr<internal::Octree_wrapper<GeomTraits> > octree,
                              const PointFunction& point_function,
                              const Gradient_& gradient = Gradient_())
{
  using Domain = Implicit_octree_domain<GeomTraits, PointFunction, Gradient_>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;
  using Gradient = typename Domain::Gradient;
  using Point_function = typename Function::element_type::Point_function;
  using Octree = typename Topology::element_type::Octree;

  const Octree oct = octree;
  const Topology topo = std::make_shared<Topology::element_type>(oct);
  const Geometry geom = std::make_shared<Geometry::element_type>(oct);
  const Point_function point_func = std::make_shared<Point_function::element_type>(point_function);
  const Function func = std::make_shared<Function::element_type>(geom, point_func);
  const Gradient grad = std::make_shared<Gradient::element_type>(gradient);

  return Domain(topo, geom, func, grad);
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
