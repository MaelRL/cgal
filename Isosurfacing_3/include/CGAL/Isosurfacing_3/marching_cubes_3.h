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

#ifndef CGAL_ISOSURFACING_3_MARCHING_CUBES_3_H
#define CGAL_ISOSURFACING_3_MARCHING_CUBES_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/marching_cubes_functors.h>
#include <CGAL/Isosurfacing_3/internal/topologically_correct_marching_cubes_functors.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Methods_grp
 *
 * \brief creates a triangle soup that represents an isosurface using the Marching Cubes algorithm.
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_if_available_tag`, or `Parallel_tag`.
 * \tparam Domain must be a model of `IsosurfacingDomain_3`.
 * \tparam PointRange must be a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type can be constructed from the point type of the domain.
 * \tparam PolygonRange must be a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                      whose value type is itself a model of the concepts `RandomAccessContainer`
 *                      and `BackInsertionSequence` whose value type is `std::size_t`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param domain the domain providing input data and its topology
 * \param isovalue value of the isosurface
 * \param points points of the triangles in the created indexed face set
 * \param triangles each element in the vector describes a triangle using the indices of the points in `points`
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_topologically_correct_marching_cubes}
 *     \cgalParamDescription{whether the topologically correct variant of Marching Cubes should be used}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename Domain,
          typename PointRange,
          typename TriangleRange,
          typename NamedParameters = parameters::Default_named_parameters>
void marching_cubes(const Domain& domain,
                    const typename Domain::Geom_traits::FT isovalue,
                    PointRange& points,
                    TriangleRange& triangles,
                    const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // @todo test 'false'
  const bool use_tmc = choose_parameter(get_parameter(np, internal_np::use_topologically_correct_marching_cubes), true);

  if(use_tmc)
  {
    // run TMC and directly write the result to points and triangles
    internal::TMC_functor<Domain, PointRange, TriangleRange> functor(domain, isovalue, points, triangles);
    domain.template iterate_cells<ConcurrencyTag>(functor);
  }
  else
  {
    // run MC
    internal::Marching_cubes_3<Domain> functor(domain, isovalue);
    domain.template iterate_cells<ConcurrencyTag>(functor);

    // copy the result to points and triangles
    internal::triangles_to_polygon_soup(functor.triangles(), points, triangles);
  }
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_MARCHING_CUBES_3_H
