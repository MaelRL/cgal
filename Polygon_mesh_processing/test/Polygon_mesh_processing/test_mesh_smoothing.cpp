#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <boost/graph/graph_traits.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef Kernel::Point_3 Point;

struct Constraints_pmap
{
  std::set<vertex_descriptor>* set_ptr_;

  typedef vertex_descriptor                   key_type;
  typedef bool                                value_type;
  typedef value_type&                         reference;
  typedef boost::read_write_property_map_tag  category;

  public:
  Constraints_pmap(std::set<vertex_descriptor>* set_ptr)
    : set_ptr_(set_ptr)
  {}
  Constraints_pmap()
    : set_ptr_(NULL)
  {}

  friend value_type get(const Constraints_pmap& map, const key_type& v)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    return !map.set_ptr_->empty()
         && map.set_ptr_->count(v);
  }

  friend void put(Constraints_pmap& map
                , const key_type& v, const value_type is)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    if (is)                map.set_ptr_->insert(v);
    else if(get(map, v))   map.set_ptr_->erase(v);
  }
 };

void test_angle_smoothing(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh);

  for (vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.7203429230262004);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 0);
      break;
    }
  }
}

void test_area_smoothing(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_areas(mesh);

  for (vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.6691415930575334);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 0);
      break;
    }
  }
}

void test_angle_smoothing_without_projection(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh, CGAL::Polygon_mesh_processing::parameters::do_project(false));

  for (vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.59571844622769954);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 1.0652302640732678);
      break;
    }
  }
}

void test_area_smoothing_without_projection(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_areas(mesh, CGAL::Polygon_mesh_processing::parameters::do_project(false));

  for (vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.42183982448892759);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 0.87816017551107273);
      break;
    }
  }
}

void test_constrained_vertices(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  double x_init, y_init, z_init;
  std::set<vertex_descriptor> selected_vertices;
  for(vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      selected_vertices.insert(v);
      x_init = get(vpmap, v).x();
      y_init = get(vpmap, v).y();
      z_init = get(vpmap, v).z();
    }
  }

  Constraints_pmap vcmap(&selected_vertices);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh,
        CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap));
  CGAL::Polygon_mesh_processing::smooth_areas(mesh,
        CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap));

  for(vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      assert(x_init == get(vpmap, v).x());
      assert(y_init == get(vpmap, v).y());
      assert(z_init == get(vpmap, v).z());
    }
  }
}

int main(){

  const char* filename_polygon = "data/simple_polygon.off";
  test_angle_smoothing(filename_polygon);
  test_area_smoothing(filename_polygon);
  test_constrained_vertices(filename_polygon);
  const char* filename_pyramid = "data/simple_pyramid.off";
  test_angle_smoothing_without_projection(filename_pyramid);
  test_area_smoothing_without_projection(filename_pyramid);

  return 0;
}
