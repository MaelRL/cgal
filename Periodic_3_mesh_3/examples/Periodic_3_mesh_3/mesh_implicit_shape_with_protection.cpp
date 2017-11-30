#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
//#define CGAL_MESH_3_DEBUG_FACET_CRITERIA
//#define CGAL_MESH_3_DEBUG_CELL_CRITERIA
#define CGAL_MESH_3_PROFILING
//#define CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
#define CGAL_MESH_3_PROTECTION_DEBUG 1111

#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Medit_IO.h>
#include <CGAL/Implicit_periodic_3_mesh_domain_3.h>
#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

namespace P3M3_IO = CGAL::Periodic_3_mesh_3::IO;

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;

// Domain
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point;
typedef K::Iso_cuboid_3                                             Iso_cuboid;

typedef FT (Function)(const Point&);
typedef CGAL::Mesh_domain_with_polyline_features_3<
          CGAL::Implicit_periodic_3_mesh_domain_3<Function,K> >     Mesh_domain;

// Polyline
typedef std::vector<Point>                                          Polyline_3;
typedef std::list<Polyline_3>                                       Polylines;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Mesh_domain>::type    Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
          Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>  C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                                   Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Periodic domain
static const int domain_size = 1;
static const Iso_cuboid periodic_domain = Iso_cuboid(0, 0, 0, domain_size, domain_size, domain_size);

// Used to artifically create a periodic function.
static const Tr dummy_tr(periodic_domain);
Point canonicalize_point(const Point& p)
{
  // Compute the representative of 'p' in the fundamental domain
  return dummy_tr.canonicalize_point(p);
}

// Implicit function
static const FT cx = 0.5, cy = 0.5, cz = 0.5;
static const FT scale = 0.9;
FT cone_function(const Point& p)
{
  const Point& cp = canonicalize_point(p);
  const FT x = cp.x(), y = cp.y(), z = cp.z();

  if (((x-cx)*(x-cx) + (y-cy)*(y-cy)) > scale * (z-cz)*(z-cz))
    return 1; // outside
  else
    return -1; // inside
}

// To obtain a good looking mesh at the base of the cone, we protect the base circle
void cone_polylines(Polylines& polylines)
{
  const FT z = 0.;
  const FT radius_at_z = CGAL::sqrt(scale * cz * cz);

  Polyline_3 polyline;
  for(int i = 0; i < 360; ++i)
  {
    polyline.push_back(Point(0.5 + radius_at_z * std::sin(i*CGAL_PI/180),
                             0.5 + radius_at_z * std::cos(i*CGAL_PI/180),
                             z));
  }
  polyline.push_back(polyline.front()); // close the line

  polylines.push_back(polyline);
}

int main()
{
  // Domain
  Mesh_domain domain(cone_function,
                     CGAL::Iso_cuboid_3<K>(0, 0, 0, domain_size, domain_size, domain_size));

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.02 * domain_size,
                         facet_angle = 0.05 * domain_size,
                         facet_size = 0.02 * domain_size,
                         cell_radius_edge_ratio = 2, cell_size = 0.5);

  // Create the features that we want to preserve
  Polylines polylines;
  cone_polylines(polylines);

  // Insert the features in the domain
  domain.add_features(polylines.begin(), polylines.end());

  // Insert a corner to make sure the apex of the cone is present in the mesh
  domain.add_corner(Point(0.5, 0.5, 0.5));

  // Mesh generation with feature preservation (and no optimizers)
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria, manifold(),
                                                     no_exude(), no_perturb());

  // Output
  std::ofstream medit_file_bis("output_protected_implicit_shape.mesh");
  P3M3_IO::write_complex_to_medit(medit_file_bis, c3t3_bis);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
