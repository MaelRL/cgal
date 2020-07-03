//#define CGAL_PMP_SNAP_DEBUG_PP
//#define CGAL_PMP_SNAP_DEBUG_OUTPUT

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/STL_reader.h>

#include <array>
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_exact_constructions_kernel             EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel           EPICK;

typedef CGAL::Polyhedron_3<EPECK, CGAL::Polyhedron_items_with_id_3>   Exact_polyhedron;
typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_with_id_3>   Polyhedron;
typedef CGAL::Surface_mesh<EPICK::Point_3>                            Surface_mesh;

template <typename Kernel, typename Mesh>
void read_mesh(const char* filename,
               Mesh& sm)
{
  typedef typename Kernel::Point_3                                    Point;

  std::ifstream in(filename, std::ios::binary);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::array<int, 3> > faces;
    CGAL::read_STL(in, points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else // off reading
  {
    if(!in || !(in >> sm))
    {
      std::cerr << "Error: cannot open mesh\n";
      return;
    }
  }
}

template <typename Kernel, typename Mesh>
void test(const char* filename,
          const double large_tolerance,
          const double good_tolerance,
          const double small_tolerance)
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;
  typedef typename Kernel::FT                                         FT;

  Mesh sm;
  read_mesh<Kernel>(filename, sm);

  std::cout << "------------------" << std::endl;
  std::cout << "num v/f: " << num_vertices(sm) << " " << num_faces(sm) << std::endl;

  std::size_t res;
  Mesh sm_cpy = sm;

  // zero tolerance, just to test the API
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_zero(0);
  PMP::experimental::snap_borders(sm_cpy, sm_cpy);
  PMP::experimental::snap_borders(sm_cpy, tol_pmap_zero, sm_cpy, tol_pmap_zero);
  PMP::experimental::snap_borders(sm_cpy, tol_pmap_zero, sm_cpy, tol_pmap_zero,
                                  CGAL::parameters::geom_traits(Kernel()),
                                  CGAL::parameters::geom_traits(Kernel()));

  // too big, creates wrong snaps
  sm_cpy = sm;
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_large(large_tolerance);
  res = PMP::experimental::snap_borders(sm_cpy, tol_pmap_large);
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out1("too_large.off");
  out1.precision(17);
  out1 << sm_cpy;
  out1.close();

  // too small
  sm_cpy = sm;
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_small(small_tolerance);
  res = PMP::experimental::snap_borders(sm_cpy, tol_pmap_small,
                                        CGAL::parameters::geom_traits(Kernel()));
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out2("too_small.off");
  out2.precision(17);
  out2 << sm_cpy;
  out2.close();

  sm_cpy = sm;
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_good(good_tolerance);
  res = PMP::experimental::snap_borders(sm_cpy, tol_pmap_good);
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out3("good.off");
  out3.precision(17);
  out3 << sm_cpy;
  out3.close();

  // automatically computed, custom tolerance at each vertex
  sm_cpy = sm;
  res = PMP::experimental::snap_borders(sm_cpy);
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out4("custom.off");
  out4.precision(17);
  out4 << sm_cpy;
  out4.close();
}

template <typename Kernel, typename Mesh>
void test_no_border()
{
  typedef typename Kernel::FT                                         FT;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor     halfedge_descriptor;

  Mesh fg_source, fg_target;

  std::ifstream source_input("/home/mrouxell/DATA/Customers/Hitachi/Smesh_1.off");
  if(!source_input || !(source_input >> fg_source))
  {
    std::cerr << "Error: cannot open source mesh\n";
    return;
  }

  std::ifstream target_input("/home/mrouxell/DATA/Customers/Hitachi/Smesh_2.off");
  if(!target_input || !(target_input >> fg_target))
  {
    std::cerr << "Error: cannot open target mesh\n";
    return;
  }

  std::vector<halfedge_descriptor> srange;
  PMP::internal::vertices_as_halfedges(vertices(fg_source), fg_source, std::back_inserter(srange));
  std::vector<halfedge_descriptor> trange;
  PMP::internal::edges_as_halfedges(edges(fg_target), fg_target, std::back_inserter(trange));

  CGAL::Constant_property_map<vertex_descriptor, FT> tol_map(0.00001);

  std::size_t res = PMP::internal::snap_vertex_edge_two_way(srange, fg_source, tol_map,
                                                            trange, fg_target, tol_map,
                                                            false /*self snapping*/,
                                                            CGAL::parameters::all_default(),
                                                            CGAL::parameters::all_default());
  std::cout << "res: " << res << " vertices" << std::endl;

  std::ofstream("out_A.off") << fg_source;
  std::ofstream("out_B.off") << fg_target;
}

void test(const char* filename,
          const double large_tolerance,
          const double good_tolerance,
          const double small_tolerance)
{
  std::cout << "######################## TEST FILE: " << filename << " ################## " << std::endl;

  std::cout << "~~~~~~~~~~~ TEST EPECK POLYHEDRON ~~~~~~~~~~~" << std::endl;
  test<EPECK, Exact_polyhedron>(filename, large_tolerance, good_tolerance, small_tolerance);

  std::cout << std::endl << "~~~~~~~~~~~ TEST EPICK POLYHEDRON ~~~~~~~~~~~" << std::endl;
  test<EPICK, Polyhedron>(filename, large_tolerance, good_tolerance, small_tolerance);

  std::cout << std::endl << "~~~~~~~~~~~ TEST EPICK SURFACE MESH ~~~~~~~~~~~" << std::endl;
  test<EPICK, Surface_mesh>(filename, large_tolerance, good_tolerance, small_tolerance);
}

int main(int, char**)
{
  test("data_snapping/non_conform_snapping.off", 0.2, 0.01, 0.0001);
  test("data_snapping/non-conform_snapping-hole.off", 0.2, 0.01, 0.0001);
  test("data_snapping/non_conform_snapping-multiple_ccs.off", 0.02, 0.01, 0.0001);
  test("data_snapping/non-conform_snapping-overlap.off", 0.2, 0.01, 0.0001);

  test("data_snapping/real_data.off", 1., 0.05, 0.0008);
  test("data_snapping/real_data_2.off", 2, 0.05, 0.000001);

  test("data_snapping/pig.stl", 20, 0.3, 0.001);

  return EXIT_SUCCESS;
}

