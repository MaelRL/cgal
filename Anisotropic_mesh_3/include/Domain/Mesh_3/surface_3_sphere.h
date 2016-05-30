#ifndef CGAL_ANISOTROPIC_MESH_3_SURFACE_3_SPHERE_H
#define CGAL_ANISOTROPIC_MESH_3_SURFACE_3_SPHERE_H

namespace CGAL
{

template<typename K>
typename K::FT sphere_function(const typename K::Point_3& p)
{
  double x2 = p.x()*p.x(), y2 = p.y()*p.y(), z2 = p.z()*p.z();
  return x2 + y2 + z2 - 1.;
}

template<typename Domain>
Domain const * sphere_domain()
{
  typedef typename Domain::R                     Traits;
  typedef typename Traits::Point_3               Point_3;
  typedef typename Traits::Sphere_3              Sphere_3;

  typename Traits::FT radius = 2. * 2.;
  return new Domain(sphere_function<Traits>,
                    Sphere_3(Point_3(1., 1., 1.), radius),
                    1e-6);
}

} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_SURFACE_3_SPHERE_H
