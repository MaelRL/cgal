// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Maxime Gimeno
#ifndef CGAL_TRANSFORM_H
#define CGAL_TRANSFORM_H
#include <CGAL/license/Polygon_mesh_processing/transform.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL{
namespace Polygon_mesh_processing{
/**
 * \ingroup PkgPolygonMeshProcessing
 * applies a `Transformation` to every vertex of a `Mesh`.
 * 
 * @tparam Transformation inherits from `CGAL::Aff_transformation_3`
 * @tparam Mesh a model of `VertexListGraph`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 * 
 * @param transformation the `Aff_transformation_3` to apply to `mesh`.
 * @param mesh the `Mesh` to transform.
 * @param np optional sequence of \ref pmp_namedparameters for `mesh`, among the ones listed below
 * 
 * * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `Mesh`\cgalParamEnd
 * \cgalNamedParamsEnd
 * 
 */
template<class Transformation, class Mesh,class NamedParameters>
void transform(const Transformation& transformation, 
               Mesh& mesh,
               const NamedParameters& np)
{
  typedef typename GetVertexPointMap<Mesh, NamedParameters>::type VPMap;
  VPMap vpm = choose_param(get_param(np, internal_np::vertex_point),
                           get_property_map(vertex_point, mesh));
  
  BOOST_FOREACH(typename boost::graph_traits<Mesh>::vertex_descriptor vd, vertices(mesh))
  {
    put(vpm, vd, transformation.transform(get(vpm, vd)));
  }
}
}
}

#endif // CGAL_TRANSFORM_H
