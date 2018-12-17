// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_PLY_WRITER_H
#define CGAL_IO_PLY_WRITER_H

#include <CGAL/IO/write_ply_points.h>

namespace CGAL{

  template <class Point_3, class Polygon_3>
  bool
  write_PLY(std::ostream& out,
            const std::vector< Point_3 >& points,
            const std::vector< Polygon_3 >& polygons,
            bool /* verbose */ = false)
  {

    if(!out)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

    // Write header
    out << "ply" << std::endl
        << ((get_mode(out) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
        << "comment Generated by the CGAL library" << std::endl
        << "element vertex " << points.size() << std::endl;
    
    internal::PLY::output_property_header (out,
                                           make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));

    out << "element face " << polygons.size() << std::endl;
  
    internal::PLY::output_property_header (out,
                                           std::make_pair (CGAL::Identity_property_map<Polygon_3>(),
                                                           PLY_property<std::vector<int> >("vertex_indices")));
    
    out << "end_header" << std::endl;
  
    for (std::size_t i = 0; i < points.size(); ++ i)
      internal::PLY::output_properties (out, points.begin() + i,
                                        make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));

    for (std::size_t i = 0; i < polygons.size(); ++ i)
      internal::PLY::output_properties (out, polygons.begin() + i,
                                        std::make_pair (CGAL::Identity_property_map<Polygon_3>(),
                                                        PLY_property<std::vector<int> >("vertex_indices")));

    return out.good();
  }

  template <class FaceListGraph>
  bool
  write_PLY(std::ostream& out,
            const FaceListGraph& mesh,
            bool /* verbose */ = false)
  {
    typedef typename boost::graph_traits<FaceListGraph>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<FaceListGraph>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<FaceListGraph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<FaceListGraph, boost::vertex_point_t>::type::value_type Point_3;
    typedef typename FaceListGraph::template Property_map<halfedge_descriptor,std::pair<float, float> > UV_map;
    UV_map h_uv;
    bool has_texture;
    boost::tie(h_uv, has_texture) = mesh.template property_map<halfedge_descriptor,std::pair<float, float> >("h:uv");
    if(!out)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

    // Write header
    out << "ply" << std::endl
        << ((get_mode(out) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
        << "comment Generated by the CGAL library" << std::endl
        << "element vertex " << num_vertices(mesh) << std::endl;
    
    internal::PLY::output_property_header (out,
                                           make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));

    out << "element face " << num_faces(mesh) << std::endl;
  
    internal::PLY::output_property_header (out,
                                           std::make_pair (CGAL::Identity_property_map<std::vector<std::size_t> >(),
                                                           PLY_property<std::vector<int> >("vertex_indices")));
    
    if(has_texture)
    {
      out << "element halfedge " << num_halfedges(mesh) << std::endl;
      
      internal::PLY::output_property_header (out,
                                             std::make_pair (CGAL::Identity_property_map<std::size_t >(),
                                                             PLY_property<unsigned int >("source")));
      
      internal::PLY::output_property_header (out,
                                             std::make_pair (CGAL::Identity_property_map<std::size_t >(),
                                                             PLY_property<unsigned int >("target")));
      internal::PLY::output_property_header (out,
                                             std::make_tuple (h_uv,
                                                              PLY_property<float>("u"),
                                                              PLY_property<float>("v")));
    }
    out << "end_header" << std::endl;

    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh))
    {
      Point_3 p = get(get(CGAL::vertex_point, mesh), vd);
      internal::PLY::output_properties (out, &p,
                                        make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));
    }
    

    std::vector<std::size_t> polygon;
    BOOST_FOREACH(face_descriptor fd, faces(mesh))
    {
      polygon.clear();
      
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, mesh), mesh))
        polygon.push_back (get(get(boost::vertex_index, mesh), source(hd,mesh)));

      internal::PLY::output_properties (out, &polygon,
                                        std::make_pair (CGAL::Identity_property_map<std::vector<std::size_t> >(),
                                                        PLY_property<std::vector<int> >("vertex_indices")));
    }
    
    if(has_texture)
    {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh))
      {
        typedef std::tuple<unsigned int, unsigned int, float, float> Super_tuple;
         Super_tuple t = 
            std::make_tuple(source(hd, mesh),target(hd, mesh),
                            h_uv[hd].first,
                            h_uv[hd].second);
        
        internal::PLY::output_properties (out, &t,
                                          std::make_pair (Nth_of_tuple_property_map<0,Super_tuple>(),
                                                          PLY_property<unsigned int >("source")),
                                          std::make_pair (Nth_of_tuple_property_map<1,Super_tuple>(),
                                                          PLY_property<unsigned int >("target")),
                                          std::make_pair (Nth_of_tuple_property_map<2,Super_tuple>(),
                                                          PLY_property<float>("u")),
                                          std::make_pair (Nth_of_tuple_property_map<3,Super_tuple>(),
                                                          PLY_property<float>("v")));
      }
    }
    return out.good();
  }

} // namespace CGAL

#endif // CGAL_IO_PLY_WRITER_H
