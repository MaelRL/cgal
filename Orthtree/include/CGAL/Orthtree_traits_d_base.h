// Copyright (c) 2023  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro

#ifndef ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_D_BASE_H
#define ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_D_BASE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_d_base` can be subclassed for easier implementation of a dd OrthtreeTraits concept.

  \tparam GeomTraits model of `Kernel`.

  \cgalModels `OrthtreeTraits`
  \sa `CGAL::Quadtree`
  \sa `CGAL::Orthtree_traits_point_d`
*/
template <typename K, typename DimensionTag>
struct Orthtree_traits_d_base {
public:

  /// \name Types
  /// @{

  using Dimension = DimensionTag;
  using FT = typename K::FT;
  using Point_d = typename K::Point_d;
  using Sphere_d = typename K::Sphere_d;
  using Cartesian_const_iterator_d = typename K::Cartesian_const_iterator_d;
  using Array = std::array<FT, Dimension::value>;

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Bbox_d; ///< Bounding box type.
#else

  class Bbox_d {
    Point_d m_min, m_max;
  public:

    Bbox_d(const Point_d& pmin, const Point_d& pmax)
      : m_min(pmin), m_max(pmax) {}

    const Point_d& min BOOST_PREVENT_MACRO_SUBSTITUTION() { return m_min; }

    const Point_d& max BOOST_PREVENT_MACRO_SUBSTITUTION() { return m_max; }
  };

#endif

  /*!
    Adjacency type.

    \note This type is used to identify adjacency directions with
    easily understandable keywords (left, right, up, etc.) and is thus
    mainly useful for `Orthtree_traits_2` and `Orthtree_traits_3`. In
    higher dimensions, such keywords do not exist and this type is
    simply an integer. Conversions from this integer to bitsets still
    work but do not provide any easier API for adjacency selection.
  */
  using Adjacency = int;

  /// @}

};

}

#endif //ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_2_BASE_H
