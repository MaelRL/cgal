#ifndef CGAL_ANISOTROPIC_MESH_3_COLORED_POLYHEDRON_OUTPUT_H
#define CGAL_ANISOTROPIC_MESH_3_COLORED_POLYHEDRON_OUTPUT_H

#include <CGAL/basic.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/IO/binary_file_io.h>
#include <iostream>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<class Writer>
void write_color( std::ostream& out, Writer& writer, double c) {
  if ( writer.header().binary())
  {
    I_Binary_write_big_endian_float32( out, float(c));
    I_Binary_write_big_endian_float32( out, float(c));
    I_Binary_write_big_endian_float32( out, float(c));
    I_Binary_write_big_endian_float32( out, float(1.0));
  }
  else
    out << " " << c << " " << c << " " << c << " " << 1.0; //rgba
}

template <class Colored_polyhedron, class Writer>
void
output_colored_polyhedron(std::ostream& out,
                          const Colored_polyhedron& P,
                          Writer& writer,
                          double strongest_color)
{
  typedef typename Colored_polyhedron::Vertex_const_iterator                  VCI;
  typedef typename Colored_polyhedron::Facet_const_iterator                   FCI;
  typedef typename Colored_polyhedron::Halfedge_around_facet_const_circulator HFCC;

  // Print header.
  writer.write_header(out,
                      P.size_of_vertices(),
                      P.size_of_halfedges(),
                      P.size_of_facets());

  for( VCI vi = P.vertices_begin(); vi != P.vertices_end(); ++vi)
  {
    writer.write_vertex(::CGAL::to_double( vi->point().x()),
                        ::CGAL::to_double( vi->point().y()),
                        ::CGAL::to_double( vi->point().z()));
  }

  typedef Inverse_index< VCI> Index;
  Index index( P.vertices_begin(), P.vertices_end());
  writer.write_facet_header();

  for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
  {
    HFCC hc = fi->facet_begin();
    HFCC hc_end = hc;
    std::size_t n = circulator_size( hc);
    CGAL_assertion( n >= 3);
    writer.write_facet_begin( n );
    do
    {
      writer.write_facet_vertex_index( index[ VCI(hc->vertex())]);
      ++hc;
    } while( hc != hc_end);

    double color = fi->color();
    color = color/strongest_color;
    write_color<Writer>(out, writer, color);

    writer.write_facet_end();
  }
  writer.write_footer();
}

}
}
#endif

