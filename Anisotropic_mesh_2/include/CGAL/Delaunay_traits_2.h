#ifndef CGAL_ANISOTROPIC_MESH_2_DELAUNAY_TRAITS_2
#define CGAL_ANISOTROPIC_MESH_2_DELAUNAY_TRAITS_2

#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/Random.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K, typename KExact = K>
class Delaunay_traits_2 : public K
{
public:
  typedef K                           Kernel;
  typedef KExact                      Exact_kernel;
  typedef typename K::FT              FT;
  typedef typename K::Point_2         Point_2;
  typedef typename K::Segment_2       Segment_2;
  typedef typename K::Triangle_2      Triangle_2;
  typedef typename K::Line_2          Line_2;
  typedef typename K::Ray_2           Ray_2;
  typedef typename K::Object_2        Object_2;

public:
  class Construct_circumcenter_2
  {
  public:
    Construct_circumcenter_2() { }

    template<typename Point>
    Point operator()(const Point& p, const Point& q, const Point& r) const
    {
      typedef typename Kernel_traits<Point>::Kernel KLocal;
      typename KLocal::Construct_circumcenter_2 o;
      return o(p,q,r);
    }

    template<typename Point>
    Point operator()(const Point& p, const Point& q) const
    {
      typedef typename Kernel_traits<Point>::Kernel KLocal;
      typename KLocal::Construct_circumcenter_2 o;
      return o(p,q);
    }
  };

public:
  class Compute_random_point_2
  {
  public:
    Compute_random_point_2() {}

    // returns a point in the sphere(center, radius)
    Point_2 operator()(const Point_2 &center, const FT &radius) const
    {
      static CGAL::Random r;
      FT px, py;
      do
      {
        px = r.get_double(-1.0, 1.0), py = r.get_double(-1.0, 1.0);
      } while (px * px + py * py > 1.0);

      return Point_2(center.x() + px * radius,
                     center.y() + py * radius);
    }

    // returns a point in the intersection of the ball(center, radius) and
    // the plane of the facet (i.e. in the "picking region" P_v(facet))
    // (center-normalp is the normal direction to the facet)
    Point_2 operator()(const Point_2 &center, const Point_2 &normalp, const FT &radius) const
    {
      static CGAL::Random r;
      FT px, py;
      do
      {
        px = r.get_double(-1.0, 1.0);
        py = r.get_double(-1.0, 1.0);
      } while (px * px + py * py > 1.0);
      FT dx = px * radius, dy = py * radius;

      FT nx = normalp.x() - center.x(),
         ny = normalp.y() - center.y();
      FT nl = sqrt(nx * nx + ny * ny);
      nx /= nl, ny /= nl;

      FT l = dx * nx + dy * ny;

      return Point_2(center.x() + dx - nx * l,
                     center.y() + dy - ny * l);
    }
  };


private:
  Construct_circumcenter_2 construct_circumcenter_2_object_cache;
  Compute_random_point_2 compute_random_point_2_object_cache;

private:
  typedef CGAL::Cartesian_converter<K, KExact> To_exact;
  typedef CGAL::Cartesian_converter<KExact, K> Back_from_exact;
  To_exact to_exact;
  Back_from_exact back_from_exact;

public:
  Delaunay_traits_2() :
    construct_circumcenter_2_object_cache(),
    compute_random_point_2_object_cache()
  {}
  ~Delaunay_traits_2(){}

public:
  const Construct_circumcenter_2 construct_circumcenter_2_object() const
  {
    return construct_circumcenter_2_object_cache;
  }
  const Construct_circumcenter_2 construct_weighted_circumcenter_2_object() const
  {
    return construct_circumcenter_2_object_cache;
  }
  const Compute_random_point_2 &compute_random_point_2_object() const
  {
    return compute_random_point_2_object_cache;
  }

}; // class Delaunay_traits_2
} // Anisotropic_mesh_2
} //namespace CGAL

#endif //CGAL_DELAUNAY_TRAITS_2
