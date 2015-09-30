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
  std::size_t m_tag; //numbering
  Eigen::Matrix3d m_metric; // !! THIS IS Mp, NOT Fp !!
  int m_contributors; //number of points of the starset in coherent zones with that vertex closest
  int m_metric_origin; //1=colored from a pt on starset, 2=from a pt within a triangle, 3=colored during spread()
  int m_colored_rank; //order of coloration during spread. -1=originally colored
  int m_red, m_green; //number of times the point was in coherent (green) or incoherent (red) zone during the probing
  bool m_colored_this_pass; // == was it already present in a coherent zone this pass

public:
  Metric_vertex()  {}
  Metric_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0),
      m_metric(Eigen::Matrix3d::Zero()),
      m_contributors(0),
      m_metric_origin(0),
      m_colored_rank(-1),
      m_red(0),
      m_green(0),
      m_colored_this_pass(false)
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
  const int& red() const { return m_red; }
  int& red() { return m_red; }
  const int& green() const { return m_green; }
  int& green() { return m_green; }
  const bool& colored_this_pass() const { return m_colored_this_pass; }
  bool& colored_this_pass() { return m_colored_this_pass; }

  bool is_colored() const {return m_metric != Eigen::Matrix3d::Zero();}
  double ratio() const { return ((double)m_green)/((double) m_green+m_red);}
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

  /*
  template<class Refs, class Traits>
  struct Face_wrapper
  {
    typedef Colored_facet<Refs, CGAL::Tag_true> Face;
  };
  */
};

}
}

#endif
