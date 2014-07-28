#ifndef IMPLICIT_SURFACE_TYPE_FWD_H
#define IMPLICIT_SURFACE_TYPE_FWD_H

#include <vector>

namespace CGAL
{
class Epick;

template<typename K>
class Robust_circumcenter_traits_3;

template<class R>
class Point_3;

namespace Anisotropic_mesh_3
{

template<typename K,
         typename Pt_container /*= std::vector<typename K::Point_3>*/ >
class Constrain_surface_3_implicit;

} // end namespace Anisotropic_mesh_3
} // end namespace CGAL

typedef CGAL::Robust_circumcenter_traits_3<CGAL::Epick> Kernel;
typedef CGAL::Point_3<CGAL::Epick> Point_3;
typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_implicit<Kernel,
                                                               std::vector<Point_3> > Implicit_surface;

#endif // IMPLICIT_SURFACE_TYPE_FWD_H
