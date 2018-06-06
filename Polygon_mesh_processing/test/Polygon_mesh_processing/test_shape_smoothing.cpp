#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/utility.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <boost/graph/graph_traits.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Mesh_with_id;

template<typename Mesh>
struct Constraints_pmap
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

  typedef vertex_descriptor                   key_type;
  typedef bool                                value_type;
  typedef value_type&                         reference;
  typedef boost::read_write_property_map_tag  category;

  std::set<vertex_descriptor>* set_ptr_;

public:
  Constraints_pmap(std::set<vertex_descriptor>* set_ptr)
    : set_ptr_(set_ptr)
  {}
  Constraints_pmap()
    : set_ptr_(NULL)
  {}

  friend value_type get(const Constraints_pmap& map, const key_type& e)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    return !map.set_ptr_->empty()
         && map.set_ptr_->count(e);
  }
  friend void put(Constraints_pmap& map
                , const key_type& e, const value_type is)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    if (is)                map.set_ptr_->insert(e);
    else if(get(map, e))   map.set_ptr_->erase(e);
  }
};

template <typename Mesh>
void test_implicit_constrained_devil(Mesh mesh)
{
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_implicit_constrained_devil --" << std::endl;
  #endif

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);
  // z max is 20 in the devil;
  std::set<vertex_descriptor> selected_vertices;
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    const double z = get(vpmap, v).z();
    if(z  > 19.0)
      selected_vertices.insert(v);
  }
  Constraints_pmap<Mesh> vcmap(&selected_vertices);

  const double time_step = 1.0;
  CGAL::Polygon_mesh_processing::smooth_along_curvature_flow(mesh, time_step,
                            CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap).
                                                                       number_of_iterations(2));
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_implicit_constrained_devil.off");
  out << mesh;
  out.close();
  #endif
}

template <typename Mesh>
void test_implicit_constrained_pyramid(Mesh mesh)
{
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_implicit_constrained_pyramid --" << std::endl;
  #endif

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  // z max is 20 in the devil;
  std::set<vertex_descriptor> selected_vertices;

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    const double z = get(vpmap, v).z();
    if(z  > 0.8)
      selected_vertices.insert(v);
  }

  Constraints_pmap<Mesh> vcmap(&selected_vertices);

  const double time_step = 1.0;
  CGAL::Polygon_mesh_processing::smooth_along_curvature_flow(mesh, time_step,
                            CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap).
                                                                       number_of_iterations(5));
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_implicit_constrained_pyramid.off");
  out << mesh;
  out.close();
  #endif
}

template <typename Mesh>
void test_explicit_scheme(Mesh mesh)
{
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_explicit_scheme --" << std::endl;
  #endif

  const double time_step = 1;
  const unsigned int iterations = 5;
  CGAL::Polygon_mesh_processing::smooth_along_curvature_flow(mesh, time_step,
                            CGAL::Polygon_mesh_processing::parameters::use_explicit_scheme(true).
                                                                       number_of_iterations(iterations));

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_explicit.off");
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
  CGAL::Polygon_mesh_processing::smooth_along_curvature_flow(mesh, time_step);

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_devil_time_step.off");
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
  CGAL::Polygon_mesh_processing::smooth_along_curvature_flow(mesh, time_step);

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_precision_pyramid.off");
  out << mesh;
  out.close();
  #endif
}

template <typename Mesh>
void test_demo_helpers(Mesh mesh)
{
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout << "-- test_demo_helpers --" << std::endl;
  #endif

  const double time_step = 1e-2;
  std::vector<CGAL::Triple<int, int, double> > stiffness;

  bool compute_stiffness = true;
  CGAL::Polygon_mesh_processing::internal::solve_mcf(faces(mesh), mesh, time_step,
                                           stiffness, compute_stiffness,
                                           CGAL::Polygon_mesh_processing::parameters::all_default());
  compute_stiffness = false;
  CGAL::Polygon_mesh_processing::internal::solve_mcf(faces(mesh), mesh, time_step,
                                           stiffness, compute_stiffness,
                                           CGAL::Polygon_mesh_processing::parameters::all_default());

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_devil_demo_helpers.off");
  out << mesh;
  out.close();
  #endif
}

int main(int argc, char* argv[])
{
  const char* filename_devil = "data/mannequin-devil.off";
  const char* filename_pyramid = "data/simple_pyramid.off";

  std::ifstream input1(filename_devil);
  SurfaceMesh mesh_devil;
  if (!input1 || !(input1 >> mesh_devil)){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  input1.close();

  std::ifstream input2(filename_pyramid);
  SurfaceMesh mesh_pyramid;
  if (!input2 || !(input2 >> mesh_pyramid)){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  input2.close();

  test_demo_helpers<SurfaceMesh>(mesh_devil);
  test_curvature_flow_time_step<SurfaceMesh>(mesh_devil);
  test_curvature_flow<SurfaceMesh>(mesh_pyramid);
  test_implicit_constrained_pyramid<SurfaceMesh>(mesh_pyramid);
  test_implicit_constrained_devil<SurfaceMesh>(mesh_devil);
  test_explicit_scheme<SurfaceMesh>(mesh_devil);

  input1.open(filename_devil);
  Mesh_with_id pl_mesh_devil;
  if (!input1 || !(input1 >> pl_mesh_devil)){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  input1.close();

  input2.open(filename_pyramid);
  Mesh_with_id pl_mesh_pyramid;
  if (!input2 || !(input2 >> pl_mesh_pyramid)){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  input2.close();

  set_halfedgeds_items_id(pl_mesh_devil);
  set_halfedgeds_items_id(pl_mesh_pyramid);

  test_demo_helpers<Mesh_with_id>(pl_mesh_devil);
  test_curvature_flow_time_step<Mesh_with_id>(pl_mesh_devil);
  test_curvature_flow<Mesh_with_id>(pl_mesh_pyramid);
  test_implicit_constrained_pyramid<Mesh_with_id>(pl_mesh_pyramid);
  test_implicit_constrained_devil<Mesh_with_id>(pl_mesh_devil);
  test_explicit_scheme<Mesh_with_id>(pl_mesh_devil);

  return 0;
}
