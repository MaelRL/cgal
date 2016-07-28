// model of the concept AABBPrimitive

#ifndef AABB_TREE_BBOX_PRIMITIVE_H
#define AABB_TREE_BBOX_PRIMITIVE_H

#include <CGAL/aabb_tree/aniso_bbox_3.h>

namespace CGAL
{
  template<typename Star>
  class AABB_bbox_primitive
  {
  public:
    // AABBPrimitive types
    typedef typename Star::Point_3              Point;
    typedef typename Star::A_Bbox_3             Datum;
    typedef Star*                               Id; // Star_handle

    // Self
    typedef AABB_bbox_primitive<Star> Self;

  private:
    /// The id, here a star handle (in the star set)
    Id m_star_handle;

  public:
    // Constructors
    AABB_bbox_primitive() {}

    AABB_bbox_primitive(const AABB_bbox_primitive& primitive)
      :
      m_star_handle(primitive.id())
    { }

    AABB_bbox_primitive(const Id& handle)
      :
      m_star_handle(handle)
    { }

    Self& operator=(const Self& s)
    {
      this->id() = s.id();
      return *this;
    }
    // Default destructor, copy constructor and assignment operator are ok

    // Returns the identifier
    const Id& id() const { return m_star_handle; }
    Id& id() { return m_star_handle; }
    
    /// Returns the geometric datum wrapped by the primitive
    Datum datum() const
    {
      return m_star_handle->AABB_bbox();
    }

    /// Returns a point on the primitive
    Point reference_point() const
    {
      return m_star_handle->AABB_center_point();
    }

  }; // end class AABB_bbox_primitive
} // end namespace CGAL

#endif // BBOX_AABB_TREE_PRIMITIVE_H
