#ifndef CGAL_ANISOTROPIC_MESH_3_ANISO_POLYHEDRON_ITEMS_H
#define CGAL_ANISOTROPIC_MESH_3_ANISO_POLYHEDRON_ITEMS_H

#include <CGAL/Polyhedron_3.h>
#include <Eigen/Dense>

#include <cstdio>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template <class Refs, class T, class P>
class Metric_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
  // tag
  std::size_t m_tag;
  Eigen::Matrix3d m_metric; // THIS IS Mp, NOT Fp
  int m_contributors;
  int m_metric_origin; //1=pt on starset, 2=pt on a triangle, 3=colored in spread
  int m_colored_rank;

public:
  Metric_vertex()  {}
  Metric_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0),
      m_metric(Eigen::Matrix3d::Zero()),
      m_contributors(0),
      m_metric_origin(0),
      m_colored_rank(-1)
  {}

  // tag
  std::size_t& tag() {  return m_tag; }
  const std::size_t& tag() const {  return m_tag; }
  Eigen::Matrix3d& metric() {  return m_metric; }
  const Eigen::Matrix3d& metric() const {  return m_metric; }
  const int& contributors() const { return m_contributors; }
  int& contributors() { return m_contributors; }
  const int& metric_origin() const { return m_metric_origin; }
  int& metric_origin() { return m_metric_origin; }
  const int& colored_rank() const { return m_colored_rank; }
  int& colored_rank() { return m_colored_rank; }

  bool is_colored() const {return m_metric != Eigen::Matrix3d::Zero();}
};

template <class Refs, class T>
class Colored_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
  std::size_t m_tag; // numbering
  double m_color;
  int m_contributors;

public:
  Colored_facet()
    : m_tag(0), m_color(0.), m_contributors(0)  {}
  const std::size_t& tag() const { return m_tag; }
  std::size_t& tag()             { return m_tag; }
  const double& color() const { return m_color; }
  double& color() { return m_color; }
  const int& contributors() const { return m_contributors; }
  int& contributors() { return m_contributors; }
};

struct Aniso_items : public CGAL::Polyhedron_items_3
{
  template<class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef Metric_vertex<Refs, CGAL::Tag_true, Point> Vertex;
  };

  template<class Refs, class Traits>
  struct Face_wrapper
  {
    typedef Colored_facet<Refs, CGAL::Tag_true> Face;
  };
};

}
}

#endif
