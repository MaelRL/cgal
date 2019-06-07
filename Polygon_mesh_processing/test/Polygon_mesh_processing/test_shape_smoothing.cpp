#define CGAL_PMP_SMOOTHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

#include <CGAL/utility.h>

#include <fstream>
#include <iostream>
#include <set>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef Kernel::Point_3                                                     Point;
typedef CGAL::Surface_mesh<Point>                                           SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>        Mesh_with_id;

namespace PMP = CGAL::Polygon_mesh_processing;

bool equal_doubles(double d1, double d2, double e)
{
  return (d1 > d2 - e) && (d1 < d2 + e);
}

template <typename Mesh>
void test_implicit_constrained_devil(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_implicit_constrained_devil --" << std::endl;
#endif

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, mesh);

  // z max is 20 in the devil
  std::set<vertex_descriptor> selected_vertices;
  for(vertex_descriptor v : vertices(mesh))
  {
    const double z = get(vpmap, v).z();
    if(is_border(v, mesh) || z > 19.0)
      selected_vertices.insert(v);
  }

  CGAL::Boolean_property_map<std::set<vertex_descriptor> > vcmap(selected_vertices);

  std::vector<Point> fixed_points(selected_vertices.size());
  std::size_t i = 0;
  for(vertex_descriptor v : selected_vertices)
    fixed_points[i++] = get(vpmap, v);

  const double time_step = 1.0;
  PMP::smooth_along_curvature_flow(mesh, time_step, CGAL::parameters::vertex_is_constrained_map(vcmap)
                                                                     .number_of_iterations(2));

  i = 0;
  for(vertex_descriptor v : selected_vertices)
  {
    const Point p = get(vpmap, v);
    CGAL_assertion(equal_doubles(p.x(), fixed_points[i].x(), 1e-10));
    CGAL_assertion(equal_doubles(p.y(), fixed_points[i].y(), 1e-10));
    CGAL_assertion(equal_doubles(p.z(), fixed_points[i].z(), 1e-10));
    ++i;
  }

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("output_implicit_constrained_devil.off");
  out << mesh;
  out.close();
#endif
}

template <typename Mesh>
void test_implicit_constrained_elephant(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_implicit_constrained_elephant --" << std::endl;
#endif

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, mesh);

  std::set<vertex_descriptor> selected_vertices;
  for(vertex_descriptor v : vertices(mesh))
  {
    const double y = get(vpmap, v).z();
    if(y > 0.14)
      selected_vertices.insert(v);
  }

  CGAL::Boolean_property_map<std::set<vertex_descriptor> > vcmap(selected_vertices);

  Point fixed_point = get(vpmap, *selected_vertices.begin());

  const double time_step = 1.0;
  PMP::smooth_along_curvature_flow(mesh, time_step,
                                   CGAL::parameters::vertex_is_constrained_map(vcmap)
                                                    .number_of_iterations(1));

  CGAL_assertion(equal_doubles(get(vpmap, *selected_vertices.begin()).x(), fixed_point.x(), 1e-14));
  CGAL_assertion(equal_doubles(get(vpmap, *selected_vertices.begin()).y(), fixed_point.y(), 1e-14));
  CGAL_assertion(equal_doubles(get(vpmap, *selected_vertices.begin()).z(), fixed_point.z(), 1e-14));

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("output_implicit_constrained_elephant.off");
  out << mesh;
  out.close();
#endif
}

template <typename Mesh>
void test_curvature_flow_time_step(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_curvature_flow_time_step --" << std::endl;
#endif

  const double time_step = 1e-15;
  PMP::smooth_along_curvature_flow(mesh, time_step);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("output_devil_time_step.off");
  out << mesh;
  out.close();
#endif
}

template <typename Mesh>
void test_curvature_flow(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_curvature_flow --" << std::endl;
#endif

  const double time_step = 1.0;
  PMP::smooth_along_curvature_flow(mesh, time_step);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("output_precision_elephant.off");
  out << mesh;
  out.close();
#endif
}

int main(int, char**)
{
  const char* filename_devil = "data/mannequin-devil.off";
  const char* filename_elephant = "data/elephant.off";

  std::ifstream input1(filename_devil);
  SurfaceMesh mesh_devil;
  if(!input1 || !(input1 >> mesh_devil))
  {
    std::cerr << "Error: can not read file.";
    return 1;
  }
  input1.close();

  std::ifstream input2(filename_elephant);
  SurfaceMesh mesh_elephant;
  if(!input2 || !(input2 >> mesh_elephant))
  {
    std::cerr << "Error: can not read file.";
    return 1;
  }
  input2.close();

  test_curvature_flow_time_step<SurfaceMesh>(mesh_devil);
  test_curvature_flow<SurfaceMesh>(mesh_elephant);
  test_implicit_constrained_elephant<SurfaceMesh>(mesh_elephant);
  test_implicit_constrained_devil<SurfaceMesh>(mesh_devil);

  input1.open(filename_devil);
  Mesh_with_id pl_mesh_devil;
  if(!input1 || !(input1 >> pl_mesh_devil)){
    std::cerr << "Error: can not read file.";
    return EXIT_FAILURE;
  }
  input1.close();

  // Polyhedron

  input2.open(filename_elephant);
  Mesh_with_id pl_mesh_elephant;
  if(!input2 || !(input2 >> pl_mesh_elephant))
  {
    std::cerr << "Error: can not read file.";
    return EXIT_FAILURE;
  }
  input2.close();

  set_halfedgeds_items_id(pl_mesh_devil);
  set_halfedgeds_items_id(pl_mesh_elephant);

  test_curvature_flow_time_step<Mesh_with_id>(pl_mesh_devil);
  test_curvature_flow<Mesh_with_id>(pl_mesh_elephant);
  test_implicit_constrained_elephant<Mesh_with_id>(pl_mesh_elephant);
  test_implicit_constrained_devil<Mesh_with_id>(pl_mesh_devil);

  return EXIT_SUCCESS;
}
