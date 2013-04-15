#ifndef IMPLICIT_SURFACE_TYPE_H
#define IMPLICIT_SURFACE_TYPE_H

// CGAL
// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef CGAL::Robust_circumcenter_traits_3<Epick> Kernel;


// surface mesh
#include <CGAL/Constrain_surface_3_implicit.h>

// simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;
typedef Kernel::Plane_3 Plane_3;

// implicit surface 
typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_implicit<Kernel> Implicit_surface;

#include "Implicit_surface_type_fwd.h"

#endif // IMPLICIT_SURFACE_TYPE_H
