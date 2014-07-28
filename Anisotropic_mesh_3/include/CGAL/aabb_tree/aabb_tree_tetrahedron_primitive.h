// model of the concept AABBPrimitive

#ifndef AABB_TREE_TETRAHEDRON_PRIMITIVE_H
#define AABB_TREE_TETRAHEDRON_PRIMITIVE_H

#include <CGAL/helpers/combinatorics_helper.h>

#include <Eigen/Dense>

#include <vector>

namespace CGAL
{

template<typename K>
struct Cell_ijkl_with_coords_access
{
  typedef typename K::Point_3         Point;
  typedef typename K::Tetrahedron_3   Tet;

  Cell_ijkl m_cell;
  std::vector<Point>* m_coords;

  Tet tetrahedron() const
  {
    return Tet((*m_coords)[m_cell.vertex(0)], (*m_coords)[m_cell.vertex(1)],
               (*m_coords)[m_cell.vertex(2)], (*m_coords)[m_cell.vertex(3)]);
  }

  Point center_point() const
  {
    return CGAL::barycenter((*m_coords)[m_cell.vertex(0)], 0.25,
                            (*m_coords)[m_cell.vertex(1)], 0.25,
                            (*m_coords)[m_cell.vertex(2)], 0.25,
                            (*m_coords)[m_cell.vertex(3)], 0.25);
  }

  Cell_ijkl_with_coords_access(const Cell_ijkl& cell_,
                               std::vector<Point>* coordinates_)
    :m_cell(cell_), m_coords(coordinates_)
  { }

  Cell_ijkl_with_coords_access(int p1, int p2, int p3, int p4,
                               std::vector<Point>* coordinates_)
    :m_cell(p1,p2,p3,p4), m_coords(coordinates_)
  { }
};

template<typename K>
class AABB_tetrahedron_primitive
{
public:
  // AABBPrimitive types
  typedef typename K::Point_3               Point;
  typedef typename K::Tetrahedron_3         Tetrahedron;
  typedef Tetrahedron                       Datum;
  typedef Cell_ijkl_with_coords_access<K>*  Id;

  // Self
  typedef AABB_tetrahedron_primitive<K>     Self;

private:
  Id m_tet;

public:
  // Constructors
  AABB_tetrahedron_primitive() {}

  AABB_tetrahedron_primitive(const AABB_tetrahedron_primitive& primitive)
    : m_tet(primitive.id()) { }

  AABB_tetrahedron_primitive(const Id& handle)
    : m_tet(handle) { }

  Self& operator=(const Self& s)
  {
    this->id() = s.id();
    return *this;
  }

  /// Returns the identifier
  const Id& id() const { return m_tet; }
  Id& id() { return m_tet; }

  /// Returns the geometric datum wrapped by the primitive
  Datum datum() const
  {
    return m_tet->tetrahedron();
  }

  /// Returns a point on the primitive
  Point reference_point() const
  {
    return m_tet->center_point();
  }

}; // end class AABB_tetrahedron_primitive
} // end namespace CGAL

#endif // AABB_TREE_TETRAHEDRON_PRIMITIVE_H
