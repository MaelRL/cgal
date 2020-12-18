// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_WITH_INFO_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_WITH_INFO_H

#include <CGAL/Polyline_tracing/Motorcycle.h>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraphTraits, typename Info>
class Motorcycle_with_info
  : public Motorcycle<MotorcycleGraphTraits>
{
  typedef Motorcycle_with_info<MotorcycleGraphTraits, Info>         Self;
  typedef Motorcycle<MotorcycleGraphTraits>                         Base;

public:
  typedef MotorcycleGraphTraits                                     Geom_traits;
  typedef typename MotorcycleGraphTraits::Triangle_mesh             Triangle_mesh;

  typedef Motorcycle_graph_node_dictionary<Geom_traits>             Nodes;
  typedef typename Nodes::Node_ptr                                  Node_ptr;

  typedef typename Base::Point_or_location                          Point_or_location;

  typedef typename Base::result_type                                result_type;

public:
  // Access
  Info&       info()       { return _info; }
  const Info& info() const { return _info; }

  // Constructor
  template<typename Tracer,
           typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
  Motorcycle_with_info(const Point_or_location& origin, const Tracer& tracer,
                       const NamedParameters& np = CGAL::parameters::all_default())
    : Base(origin, tracer, np),
      _info(parameters::choose_parameter(parameters::get_parameter(np, internal_np::info), Info()))
  { }

  // See explanation in Motorcycle.h

  // disable copy operators
  Motorcycle_with_info& operator=(const Motorcycle_with_info& other) = delete;
  Motorcycle_with_info(const Motorcycle_with_info& other) = delete;

  // disable move operators
  Motorcycle_with_info (Motorcycle_with_info&& other) = delete;
  Motorcycle_with_info& operator=(Motorcycle_with_info&& other) = delete;

public:

private:
  Info _info;
};

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_WITH_INFO_H
