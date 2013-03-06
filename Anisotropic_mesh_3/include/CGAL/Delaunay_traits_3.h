#ifndef CGAL_DELAUNAY_TRAITS_3
#define CGAL_DELAUNAY_TRAITS_3

#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {
    template<typename K, typename KExact = K>
    class Delaunay_traits_3 : public K
    {
    public:
      typedef typename K::FT              FT;
      typedef typename K::Point_3         Point_3;
      typedef typename K::Segment_3       Segment_3;
      typedef typename K::Triangle_3      Triangle_3;
      typedef typename K::Tetrahedron_3   Tetrahedron_3;
      typedef typename K::Line_3          Line_3;
      typedef typename K::Ray_3           Ray_3;
      typedef typename K::Object_3        Object_3;
      typedef KExact Exact_kernel;

    public:
      class Construct_circumcenter_3 
      {
      public:
        Construct_circumcenter_3() { }
        template<typename Point>
        Point operator()(const Point &p, const Point &q, const Point &r, const Point &s) const 
        {
          typedef typename Kernel_traits<Point>::Kernel KLocal;
          typename CGAL::Robust_circumcenter_traits_3<KLocal>::Construct_circumcenter_3 o;
          return o(p,q,r,s);
        }
        template<typename Point>
        Point operator()(const Point &p, const Point &q, const Point &r) const 
        {
          typedef typename Kernel_traits<Point>::Kernel KLocal;
          typename CGAL::Robust_circumcenter_traits_3<KLocal>::Construct_circumcenter_3 o;
          return o(p,q,r);
        }
        template<typename Point>
        Point operator()(const Point &p, const Point &q) const 
        {
          typedef typename Kernel_traits<Point>::Kernel KLocal;
          typename CGAL::Robust_circumcenter_traits_3<KLocal>::Construct_circumcenter_3 o;
          return o(p,q);
        }
      };

    public:
      class Compute_volume 
      {
      public:
        Compute_volume(){ }
        FT operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s) const
        {
          FT x[3], y[3], z[3];
          Point_3 tp[4];
          tp[0] = p; tp[1] = q; tp[2] = r; tp[3] = s;
          for(int i = 1; i <= 3; i++)
          {
            x[i - 1] = tp[i].x() - tp[0].x();
            y[i - 1] = tp[i].y() - tp[0].y();
            z[i - 1] = tp[i].z() - tp[0].z();
          }
          FT determinant = x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1] 
          - x[0] * y[2] * z[1] - x[2] * y[1] * z[0] - x[1] * y[0] * z[2];
          return fabs(determinant / 6.0);
        }
        FT operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const 
        {
          Point_3 tp[3];
          tp[0] = p; tp[1] = q; tp[2] = r;
          FT l[3];
          for (int i = 0; i < 3; i++) 
          {
            int j = (i + 1) % 3;
            l[i] = sqrt((tp[i].x() - tp[j].x()) * (tp[i].x() - tp[j].x()) +
              (tp[i].y() - tp[j].y()) * (tp[i].y() - tp[j].y()) +
              (tp[i].z() - tp[j].z()) * (tp[i].z() - tp[j].z()));
          }
          FT s = (l[0] + l[1] + l[2]) / 2.0;
          return sqrt(s * (s - l[0]) * (s - l[1]) * (s - l[2]));
        }
      };

    public:
      class Compute_random_point_3 
      {
      public:
        Compute_random_point_3() {}

        // returns a point in the sphere(center, radius)
        Point_3 operator()(const Point_3 &center, const FT &radius) const 
        {
          static CGAL::Random r;
          FT px, py, pz;
          do
          { 
            px = r.get_double(-1.0, 1.0), py = r.get_double(-1.0, 1.0), pz = r.get_double(-1.0, 1.0);
          } while (px * px + py * py + pz * pz > 1.0);

          return Point_3(center.x() + px * radius, 
                         center.y() + py * radius, 
                         center.z() + pz * radius);
        }

        // returns a point in the intersection of the ball(center, radius) and
        // the plane of the facet (i.e. in the "picking region" P_v(facet))
        // (center-normalp is the normal direction to the facet)
        Point_3 operator()(const Point_3 &center, const Point_3 &normalp, const FT &radius) const 
        {
          static CGAL::Random r;
          FT px, py, pz;
          do
          {
            px = r.get_double(-1.0, 1.0);
            py = r.get_double(-1.0, 1.0);
            pz = r.get_double(-1.0, 1.0);
          } while (px * px + py * py + pz * pz > 1.0);
          FT dx = px * radius, dy = py * radius, dz = pz * radius;

          FT nx = normalp.x() - center.x(), 
             ny = normalp.y() - center.y(), 
             nz = normalp.z() - center.z();
          FT nl = sqrt(nx * nx + ny * ny + nz * nz);
          nx /= nl, ny /= nl, nz /= nl;
                    
          FT l = dx * nx + dy * ny + dz * nz;

          return Point_3(center.x() + dx - nx * l, 
                         center.y() + dy - ny * l, 
                         center.z() + dz - nz * l);
        }
      };


    private:
      Construct_circumcenter_3 construct_circumcenter_3_object_cache;
      Compute_volume compute_volume_object_cache;
      Compute_random_point_3 compute_random_point_3_object_cache;

    private:
      typedef CGAL::Cartesian_converter<K, KExact> To_exact;
      typedef CGAL::Cartesian_converter<KExact, K> Back_from_exact;
      To_exact to_exact;
      Back_from_exact back_from_exact;
      
    public:
      Delaunay_traits_3() :
          construct_circumcenter_3_object_cache(),
          compute_volume_object_cache(),
          compute_random_point_3_object_cache()
          {}
      ~Delaunay_traits_3(){}

    public:
      const Construct_circumcenter_3 construct_circumcenter_3_object() const 
      { 
        return construct_circumcenter_3_object_cache; 
      }
      const Compute_volume &compute_volume_object() const 
      { 
        return compute_volume_object_cache; 
      }      
      const Compute_random_point_3 &compute_random_point_3_object() const 
      { 
        return compute_random_point_3_object_cache; 
      }


    }; // class Delaunay_traits_3
  } // Anisotropic_mesh_3
} //namespace CGAL

#endif //CGAL_DELAUNAY_TRAITS_3
