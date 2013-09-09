#ifndef C3T3_POLYHEDRON_BUILDER_H
#define C3T3_POLYHEDRON_BUILDER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace CGAL {

template <class C3T3, class Polyhedron_>
class Complex_3_in_triangulation_3_polyhedron_builder
  : public CGAL::Modifier_base<typename Polyhedron_::HalfedgeDS>
{
public:
  typedef C3T3                                      C3t3;
  typedef Polyhedron_                               Polyhedron;
  typedef typename Polyhedron::HalfedgeDS           HDS;

private:
  typedef typename C3t3::Triangulation              Tr;
  typedef typename Tr::Vertex_handle                Vertex_handle;
  typedef typename Tr::Point                        Point;
  typedef typename Tr::Facet                        Facet;
  typedef typename Tr::Finite_vertices_iterator     Finite_vertices_iterator;
  typedef typename C3t3::Facets_in_complex_iterator Facets_in_complex_iterator;
  typedef typename C3t3::Subdomain_index            Subdomain_index;

  const C3t3& c3t3;
  typedef CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> Base;

public:
  Complex_3_in_triangulation_3_polyhedron_builder(const C3t3& c3t3_)
    : Base(), c3t3(c3t3_)
  {
  }

  void operator()(HDS& hds)
  {
    const Tr& tr = c3t3.triangulation();
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
    const typename Tr::size_type number_of_facets = c3t3.number_of_facets_in_complex();
    builder.begin_surface(tr.number_of_vertices(), number_of_facets);

    // Finite vertices coordinates.
    std::map<Vertex_handle, int> V;
    int inum = 0;
    for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
      V[vit] = inum++;
      Point p = static_cast<Point>(vit->point());
      builder.add_vertex(p);
    }

    for( Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin();
         fit != c3t3.facets_in_complex_end(); ++fit)
    {
      Facet facet = *fit;
      Subdomain_index def_sub_index = Subdomain_index();
      if(facet.first->subdomain_index() == def_sub_index)
        facet = tr.mirror_facet(facet);

      builder.begin_facet();
      builder.add_vertex_to_facet(V[facet.first->vertex(tr.vertex_triple_index(facet.second, 0))]);
      builder.add_vertex_to_facet(V[facet.first->vertex(tr.vertex_triple_index(facet.second, 1))]);
      builder.add_vertex_to_facet(V[facet.first->vertex(tr.vertex_triple_index(facet.second, 2))]);
      builder.end_facet();
    }
    builder.end_surface();
  }

}; // end Complex_3_in_triangulation_3_polyhedron_builder

} // end namespace CGAL

#endif //C3T3_POLYHEDRON_BUILDER_H
