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

#ifndef CGAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Tmc_internal.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Creates a triangular indexed face set that represents an isosurface using the marching cubes algorithm.
 *
 * \details
 *
 * \tparam ConcurrencyTag determines if the algorithm is executed sequentially or in parallel. Default is sequential.
 *
 * \tparam Domain_ must be a model of `IsosurfacingDomain`.
 *
 * \tparam PointRange is a model of the concept `RandomAccessContainer` and `BackInsertionSequence` whose value type can
 * be constructed from the point type of the domain.
 *
 * \tparam PolygonRange a model of the concept
 * `RandomAccessContainer` and `BackInsertionSequence` whose value type is itself a model of the concepts
 * `RandomAccessContainer` and `BackInsertionSequence` whose value type is `std::size_t`.
 *
 * \param domain the domain providing input data and its topology
 * \param iso_value value of the isosurface
 * \param points points of the triangles in the created indexed face set
 * \param triangles each element in the vector describes a triangle using the indices of the points in `points`
 * \param topologically_correct decides whether the topologically correct variant of Marching Cubes should be used
 */
template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class TriangleRange>
void marching_cubes(const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                    TriangleRange& triangles, bool topologically_correct = true) {

    if (topologically_correct) {
        // run TMC and directly write the result to points and triangles
        internal::TMC_functor<Domain_, PointRange, TriangleRange> functor(domain, iso_value, points, triangles);
        domain.template iterate_cells<Concurrency_tag>(functor);
    } else {
        // run MC
        internal::Marching_cubes_functor<Domain_> functor(domain, iso_value);
        domain.template iterate_cells<Concurrency_tag>(functor);
        // copy the result to points and triangles
        internal::to_indexed_face_set(functor.get_triangles(), points, triangles);
    }
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_CUBES_3_H
