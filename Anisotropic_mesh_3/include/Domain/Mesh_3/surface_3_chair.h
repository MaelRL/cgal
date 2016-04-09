#ifndef CGAL_ANISOTROPIC_MESH_3_SURFACE_3_CHAIR_H
#define CGAL_ANISOTROPIC_MESH_3_SURFACE_3_CHAIR_H

namespace CGAL
{

template<typename K>
typename K::FT chair_function(const typename K::Point_3& p)
{
  typedef typename K::FT FT;

  FT x = p.x();
  FT y = p.y();
  FT z = p.z();
  const FT a = 0.8;
  const FT b = 0.4;
  const FT k = 1.0;

  return (x*x+y*y+z*z-a*k*k) * (x*x+y*y+z*z-a*k*k)
          - b * (((z-k)*(z-k)-2*x*x) * ((z+k)*(z+k)-2*y*y));
}

template<typename Domain>
Domain const * chair_domain()
{
  typedef typename Domain::R                     Traits;
  typedef typename Traits::Point_3               Point_3;
  typedef typename Traits::Sphere_3              Sphere_3;

  typename Traits::FT radius = 10.;

  return new Domain(chair_function<Traits>,
                    Sphere_3(Point_3(0., 0., 0.), radius),
                    1e-6);
}

} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_SURFACE_3_chair_H
