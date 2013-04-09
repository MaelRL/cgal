#ifndef IMPLICIT_SURFACE_TYPE_FWD_H
#define IMPLICIT_SURFACE_TYPE_FWD_H

namespace CGAL {

  class Epick;

  namespace Anisotropic_mesh_3 {

    template<typename K>
    class Constrain_surface_3_implicit;

  } // end namespace Anisotropic_mesh_3

} // end namespace CGAL

typedef CGAL::Epick Kernel;
typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_implicit<Kernel> Implicit_surface;

#endif // IMPLICIT_SURFACE_TYPE_FWD_H
