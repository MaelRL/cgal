#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_distance_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/partition.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <tbb/task_group.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                                      Kernel;

typedef Kernel::Point_3                                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                                         Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::halfedge_descriptor             halfedge_descriptor;
typedef boost::graph_traits<Triangle_mesh>::face_descriptor                 face_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef typename SMS::GarlandHeckbert_policies<Triangle_mesh, Kernel>       GH_policies;
typedef typename GH_policies::Get_cost                                      GH_cost;
typedef typename GH_policies::Get_placement                                 GH_placement;

typedef typename Kernel::Triangle_3                                         Triangle;
typedef std::vector<Triangle>                                               Triangle_container;
typedef typename Triangle_container::iterator                               TC_iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, TC_iterator>                  Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive>                                Traits;
typedef CGAL::AABB_tree<Traits>                                             AABB_tree;

struct Border_is_constrained_edge_map
{
  typedef Triangle_mesh::Edge_index edge_descriptor;
  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  Border_is_constrained_edge_map(const Triangle_mesh& sm)
    :
      sm(sm)
  { }

  friend bool get(Border_is_constrained_edge_map m, const edge_descriptor e) {
    // @todo use flags instead
    return (CGAL::is_border(source(e, m.sm), m.sm) || CGAL::is_border(target(e, m.sm), m.sm)); // 1 ring
    // return CGAL::is_border(e, m.sm); // only border
  }

private:
  const Triangle_mesh& sm;
};

template <typename FacePatchMap>
struct Patch_boundary_is_constrained_edge_map
{
  typedef Triangle_mesh::Edge_index edge_descriptor;

  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  Patch_boundary_is_constrained_edge_map(const Triangle_mesh& sm,
                                         const FacePatchMap fpatch_map)
    : sm(sm), fpatch_map(fpatch_map)
  { }

  friend bool get(Patch_boundary_is_constrained_edge_map m, const edge_descriptor e)
  {
    const Triangle_mesh& sm = m.sm;
    if(is_border(e, sm))
      return false;

    return get(m.fpatch_map, face(halfedge(e, sm), sm)) ==
           get(m.fpatch_map, face(opposite(halfedge(e, sm), sm), sm));
  }

  const Triangle_mesh& sm;
  FacePatchMap fpatch_map;
};

struct Simplify
{
  Simplify(Triangle_mesh& tm, const double edge_length, const double bounding_distance, const AABB_tree& tree)
    : tm_ptr(&tm), edge_length(edge_length), bounding_distance(bounding_distance), tree(tree)
  {}

  void operator()() const
  {
    // There is no need to constrain the placement here because we have (combinatorially)
    // constrained a 1-ring around the border, hence there can be no geometrical shift
    typedef SMS::Bounded_distance_placement<GH_placement, AABB_tree>   Bounded_GH_placement;
    typedef SMS::Bounded_normal_change_placement<Bounded_GH_placement> Placement;

    GH_policies gh_policies(*tm_ptr);
    const GH_cost& gh_cost = gh_policies.get_cost();
    const GH_placement& gh_placement = gh_policies.get_placement();

    Border_is_constrained_edge_map ecm(*tm_ptr);
    Placement placement(Bounded_GH_placement(bounding_distance, tree, gh_placement));

    std::cout << "bounding_distance: " << bounding_distance << std::endl;
    std::cout << "edge_length: " << edge_length << std::endl;

    SMS::Edge_length_stop_predicate<double> stop(edge_length);

    SMS::edge_collapse(*tm_ptr, stop,
                       CGAL::parameters::edge_is_constrained_map(ecm)
                                        .get_placement(placement)
                                        .get_cost(gh_cost));
  }

private:
  Triangle_mesh* tm_ptr;
  double edge_length;
  double bounding_distance;
  const AABB_tree& tree;
};

int main(int argc, char** argv)
{
  Triangle_mesh tm;

  std::ifstream in((argc>1) ? argv[1] : "data/elephant.off");
  if(!in || !(in >> tm))
  {
    std::cerr << "Failed to read input mesh: " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: #V = " << num_vertices(tm)
                  << " #E = " << num_edges(tm)
                  << " #F = " << num_faces(tm) << " sec." << std::endl;

  // The idea is to let the edge_length be large, and (ab)use the bounded_distance placement wrapper
  // such that it becomes a kind of stopping criterion. Since this is a local value, it is
  // compatible with parallel approaches (on the contrary to other stop criteria such as
  // % of initial edges)

  // can't use edge_length which is also local because the queue is sorted by GH cost,
  // not by edge length and the inconsistency means bad results
  const double edge_length = (argc>2) ? boost::lexical_cast<double>(argv[2]) : 0.01;
  const double bounding_distance = (argc>3) ? boost::lexical_cast<double>(argv[3]) : 0.001;
  const int number_of_parts = (argc>4) ? boost::lexical_cast<int>(argv[4]) : 8;

  CGAL::Real_timer t;
  t.start();

  // Build the tree that will be used in distance bound checks
  std::vector<Triangle> input_triangles;

  for(face_descriptor f : faces(tm))
  {
    halfedge_descriptor h = halfedge(f, tm);
    CGAL_assertion(!is_border(h, tm));

    input_triangles.push_back(Triangle(tm.point(source(h, tm)),
                                       tm.point(target(h, tm)),
                                       tm.point(target(next(h, tm), tm))));
  }

  AABB_tree tree(input_triangles.begin(), input_triangles.end());

  Triangle_mesh final_mesh;

  if(number_of_parts > 1)
  {
    // Partition ID map
    Triangle_mesh::Property_map<face_descriptor, std::size_t> fpmap
        = tm.add_property_map<face_descriptor, std::size_t>("f:partition").first;

    // Set some custom options for METIS
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1; // need to have contiguous subdomains
    options[METIS_OPTION_UFACTOR] = 1;

    CGAL::METIS::partition_dual_graph(tm, number_of_parts,
                                      CGAL::parameters::METIS_options(&options)
                                      .face_partition_id_map(fpmap));

    std::cerr << "Built partition in " << t.time() << " sec."<< std::endl;
    t.reset();

    // Simplify each subregion in parallel
    std::vector<Triangle_mesh> meshes(number_of_parts);
    CGAL::Face_filtered_graph<Triangle_mesh> fg(tm, 0, fpmap);
    for(int i=0; i<number_of_parts; ++i)
    {
      if(i != 0)
        fg.set_selected_faces(i, fpmap);

      CGAL::copy_face_graph(fg, meshes[i]);
      std::cout << "Patch of size: " << num_vertices(meshes[i]) << "nv and "
                << num_faces(meshes[i]) << " nf" << std::endl;
    }
    std::cerr << "Meshes copied in " << t.time() << "s."<< std::endl;

    // Parallel block
    tbb::task_group tasks;
    for(int id=0; id<number_of_parts; ++id)
      tasks.run(Simplify(meshes[id], edge_length, bounding_distance, tree));

    tasks.wait();

    std::cerr << "  Total time with parallel simplification: " << t.time() << "s."<< std::endl;

    for(int i=0; i<number_of_parts; ++i)
      std::cout << "Patch of size: " << vertices(meshes[i]).size() << "nv and "
                                     << faces(meshes[i]).size() << " nf" << std::endl;

    std::size_t nv=0, nf=0, ne=0;
    for(int i=0; i< number_of_parts; ++i)
    {
      meshes[i].collect_garbage();
      nv += num_vertices(meshes[i]);
      nf += num_faces(meshes[i]);
      ne += num_edges(meshes[i]);
    }

    final_mesh.reserve(nv, ne, nf);

    Triangle_mesh::Property_map<face_descriptor, std::size_t> fpatch_map
        = final_mesh.add_property_map<face_descriptor, std::size_t>("f:patch_id").first;

    for(int i=0; i<number_of_parts; ++i)
      CGAL::copy_face_graph(meshes[i], final_mesh,
                            CGAL::parameters::face_to_face_output_iterator(
                              boost::make_function_output_iterator(
                                [&fpatch_map, i]
                                (const std::pair<Triangle_mesh::Face_index, Triangle_mesh::Face_index>& p) {
                                  fpatch_map[p.second] = i;
                                }
                              )));

    std::ofstream("reassembled.off") << final_mesh;

    PMP::stitch_borders(final_mesh);
    std::cerr << "Total time with mesh reassembly:" << t.time() << "s."<< std::endl;
  }
  else
  {
    final_mesh = tm; // @fixme ugly
  }

  CGAL_assertion(CGAL::is_valid_polygon_mesh(final_mesh));

  std::cout << "POST FIRST PASS: #V = " << num_vertices(final_mesh)
                            << " #E = " << num_edges(final_mesh)
                            << " #F = " << num_faces(final_mesh) << " sec." << std::endl;

  // final pass
  GH_policies gh_policies(final_mesh);
  const GH_cost& gh_cost = gh_policies.get_cost();
  const GH_placement& gh_placement = gh_policies.get_placement();

  SMS::Edge_length_stop_predicate<double> stop(edge_length);

  typedef SMS::Bounded_distance_placement<GH_placement, AABB_tree>                     Bounded_GH_placement;

  // Could - on paper - constrain the interior of patches to improve speed

//  typedef Patch_boundary_is_constrained_edge_map<decltype(fpatch_map)>                 Patch_boundary_ECM;
//  typedef SMS::Constrained_placement<Bounded_GH_placement,
//                                     Patch_boundary_ECM>                               Constrained_bounded_GH_placement;
//  typedef SMS::Bounded_distance_placement<Constrained_bounded_GH_placement, AABB_tree> Placement;
//  Patch_boundary_ECM ecm(final_mesh, fpatch_map);
//  Placement placement(bounding_distance, tree, Constrained_bounded_GH_placement(ecm, Bounded_GH_placement(gh_placement)));
  // and re-add CGAL::parameters::edge_is_constrained_map(ecm)

  typedef SMS::Bounded_normal_change_placement<Bounded_GH_placement> Placement;
  Placement placement(Bounded_GH_placement(bounding_distance, tree, gh_placement));

  SMS::edge_collapse(final_mesh, stop,
                     CGAL::parameters::get_placement(placement)
                                      .get_cost(gh_cost));

  final_mesh.collect_garbage();

  std::cerr << "  Total time: " << t.time() << " sec.\n"
            << "  Result: #V = " << num_vertices(final_mesh)
            << " #E = " << num_edges(final_mesh)
            << " #F = " << num_faces(final_mesh) << std::endl;

  std::ofstream("simplified_v1_mesh.off") << final_mesh;
}
