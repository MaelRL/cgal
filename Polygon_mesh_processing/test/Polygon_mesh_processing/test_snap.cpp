#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>

#include <iostream>
#include <fstream>
#include <set>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_exact_constructions_kernel             EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel           EPICK;

typedef CGAL::Polyhedron_3<EPECK, CGAL::Polyhedron_items_with_id_3>   Exact_polyhedron;
typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_with_id_3>   Polyhedron;
typedef CGAL::Surface_mesh<EPICK::Point_3>                            Surface_mesh;

template <typename Kernel, typename Mesh>
void test()
{
  typedef typename Kernel::FT                                         FT;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;

  Mesh fg_source, fg_target;

  // empty meshes
  std::size_t res = PMP::internal::snap_border(fg_source, fg_target);
  assert(res == 0);

  std::ifstream source_input("data_snapping/border_snapping_source.off");
  if(!source_input || !(source_input >> fg_source))
  {
    std::cerr << "Error: cannot open source mesh\n";
    return;
  }

  // one empty mesh
  res = PMP::internal::snap_border(fg_source, fg_target);
  assert(res == 0);

  res = PMP::internal::snap_border(fg_target, fg_source);
  assert(res == 0);

  std::ifstream target_input("data_snapping/border_snapping_target.off");
  if(!target_input || !(target_input >> fg_target))
  {
    std::cerr << "Error: cannot open target mesh\n";
    return;
  }

  Mesh fg_source_cpy = fg_source;

  // this epsilon value is too small, nothing happens
  std::cout << "*********************** EPS = 0.000000001 *************** " << std::endl;
  res = PMP::internal::snap_border(fg_source_cpy, fg_target, 0.000000001);
  std::cout << "Moved: " << res << " vertices" << std::endl;
  assert(res == 0);

  // this epsilon value is too big, and the snapping is rejected
  std::cout << "*********************** EPS = 0.1 *************** " << std::endl;
  fg_source_cpy = fg_source;
  res = PMP::internal::snap_border(fg_source_cpy, fg_target, 0.1,
                                   params::geom_traits(Kernel()), params::all_default());
  std::cout << "Moved: " << res << " vertices" << std::endl;
  assert(res == 0);

  // this is a good value of 'epsilon', but not all expected vertices are projected
  // because the sampling of the border of the source mesh is not uniform
  std::cout << "*********************** EPS = 0.001 *************** " << std::endl;
  fg_source_cpy = fg_source;
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_map(0.001);
  res = PMP::internal::snap_border(fg_source_cpy, fg_target, tol_map,
                                   params::all_default(), params::all_default());
  std::cout << "Moved: " << res << " vertices" << std::endl;
  assert(res == 76);

  std::ofstream partial_snap_out("partially_snapped_mesh.off");
  partial_snap_out << std::setprecision(17) << fg_source_cpy;

  // this one automatically computes an epsilon bound at each vertex
  std::cout << "*********************** EPS = LOCALLY COMPUTED *************** " << std::endl;
  fg_source_cpy = fg_source;
  res = PMP::internal::snap_border(fg_source_cpy, fg_target);
  std::cout << "Moved: " << res << " vertices" << std::endl;
  assert(res == 77);

  std::ofstream full_snap_out("snapped_mesh.off");
  full_snap_out << std::setprecision(17) << fg_source_cpy;
}

int main(int, char**)
{
  std::cout << "TEST EPECK POLYHEDRON" << std::endl;
  test<EPECK, Exact_polyhedron>();

  std::cout << std::endl << "TEST EPICK POLYHEDRON" << std::endl;
  test<EPICK, Polyhedron>();

  std::cout << std::endl << "TEST EPICK SURFACE MESH" << std::endl;
  test<EPICK, Surface_mesh>();

  return EXIT_SUCCESS;
}

