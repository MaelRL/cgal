#ifndef POLYHEDRAL_SURFACE_TYPE_FWD_H
#define POLYHEDRAL_SURFACE_TYPE_FWD_H

#if 1//def USE_FORWARD_DECL

#include <CGAL/Filtered_kernel_fwd.h>
#include <memory>

#include "Polyhedron_type_fwd.h"

namespace CGAL
{

namespace Anisotropic_mesh_3
{

template<typename K,
         typename Polyhedron /* = CGAL::Polyhedron_3<K, Enriched_items> */ >
class Constrain_surface_3_polyhedral;

} // end namespace Anisotropic_mesh_3
} // end namespace CGAL

// kernel
typedef CGAL::Robust_circumcenter_traits_3<CGAL::Epick> Kernel;

// surface mesh
typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_polyhedral<Kernel,
                                                                 Polyhedron> Polyhedral_surface;
#else
#include "Polyhedral_surface_type.h"
#endif

#endif // POLYHEDRAL_SURFACE_TYPE_FWD_H
