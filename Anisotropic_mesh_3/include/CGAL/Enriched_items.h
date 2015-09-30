#ifndef CGAL_ANISOTROPIC_MESH_3_ENRICHED_ITEMS_H
#define CGAL_ANISOTROPIC_MESH_3_ENRICHED_ITEMS_H

#include <CGAL/Origin.h>
#include <CGAL/Polyhedron_items_3.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

// a refined facet with a tag
template <class Refs, class T>
class Enriched_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
  // tag
  std::size_t m_tag;

public:
  // life cycle
  // no constructors to repeat, since only
  // default constructor mandatory
  Enriched_facet()
    : m_tag(0)  {}
  // tag
  const std::size_t& tag() const { return m_tag; }
  std::size_t& tag() { return m_tag; }
};

// a refined halfedge with a tag
template <class Refs, class Tprev, class Tvertex, class Tface>
class Enriched_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:
  // general purpose tag
  std::size_t m_tag;

public:
  // life cycle
  Enriched_halfedge()
    : m_tag(0)  {}
  // tag
  const std::size_t& tag() const { return m_tag;  }
  std::size_t& tag()             { return m_tag;  }
};

// a refined vertex with a tag
template <class Refs, class T, class P, class V>
class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
  typedef V Vector;

  // tag
  std::size_t m_tag;
  Vector m_normal;

public:
  // life cycle
  Enriched_vertex()  {}
  // repeat mandatory constructors
  Enriched_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0), m_normal(CGAL::NULL_VECTOR) {}

  Enriched_vertex(const P& pt, const Vector& n)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0), m_normal(n) {}

  // tag
  std::size_t& tag() {  return m_tag; }
  const std::size_t& tag() const {  return m_tag; }
  // normal
  Vector& normal() {  return m_normal; }
  const Vector& normal() const {  return m_normal; }
};

struct Enriched_items : public CGAL::Polyhedron_items_3
{
  // wrap vertex
  template<class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef Enriched_vertex<Refs, CGAL::Tag_true, Point, Vector> Vertex;
  };

  // wrap face
  template<class Refs, class Traits>
  struct Face_wrapper
  {
    typedef Enriched_facet<Refs, CGAL::Tag_true> Face;
  };

  // wrap halfedge
  template<class Refs, class Traits>
  struct Halfedge_wrapper
  {
    typedef Enriched_halfedge<Refs,
                              CGAL::Tag_true,
                              CGAL::Tag_true,
                              CGAL::Tag_true> Halfedge;
  };
};


}
}

#endif
