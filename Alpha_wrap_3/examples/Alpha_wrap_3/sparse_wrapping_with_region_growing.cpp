#define CGAL_AW3_DEBUG_QUEUE
#define CGAL_AW3_DEBUG_DUMP_EVERY_STEP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <string>

namespace AW3 = CGAL::Alpha_wraps_3;
namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using FT = K::FT;

using Mesh = CGAL::Surface_mesh<Point_3>;

using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

template <typename RegionMap>
struct Filtering_visitor
  : public AW3::internal::Wrapping_default_visitor
{
  using Region_ID = typename boost::property_traits<RegionMap>::value_type;

  using AABB_face_graph_primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using AABB_face_graph_traits = CGAL::AABB_traits<K, AABB_face_graph_primitive>;
  using Tree = CGAL::AABB_tree<AABB_face_graph_traits>;

  const Mesh& mesh;
  RegionMap rfm;
  Tree tree; // @todo There's already a tree in the oracle

public:
  Filtering_visitor(const Mesh& mesh,
                    RegionMap rfm)
    : mesh(mesh), rfm(rfm)
  {
    PMP::build_AABB_tree(mesh, tree);
  }

public:
  template <typename Wrapper, typename Facet>
  bool consider_facet(const Wrapper& wrapper,
                      const Facet& f) const
  {
    if(wrapper.triangulation().is_infinite(f))
      return true;

    std::array<std::set<Region_ID>, 3> region_ids; // one per vertex of the face

    auto fill_set = [&](const int i) -> bool
    {
      const Point_3& p = wrapper.triangulation().point(f.first, (f.second + 1 + i)%4);

      PMP::Face_location<Mesh, FT> loc = PMP::locate_with_AABB_tree(p, tree, mesh);

      // avoid numerical errors creating negative barycentric coordinates
      PMP::internal::snap_location_to_border(loc, mesh, 1e-10);

      // only check faces whose points are on the offset
      // @todo this could simply be a check on the vertex type, aka "if f->vertex(i)->type is DEFAULT"
      const Point_3 cp = PMP::construct_point(loc, mesh);
      if(CGAL::squared_distance(p, cp) > 1.5 * CGAL::square(wrapper.offset()))
        return false;

      // If the closest point on the input is on a vertex or an edge,
      // we need to get all the incident regions
      std::set<face_descriptor> ifs;
      PMP::internal::incident_faces(loc, mesh, std::inserter(ifs, ifs.begin()));

      for(const face_descriptor fd : ifs)
        region_ids[i].insert(get(rfm, fd));

      return true;
    };

    if(!fill_set(0))
      return true;
    if(!fill_set(1))
      return true;
    if(!fill_set(2))
      return true;

    // check if there is a common region within the three sets
    std::set<Region_ID> common01;
    std::set_intersection(region_ids[0].begin(), region_ids[0].end(),
                          region_ids[1].begin(), region_ids[1].end(),
                          std::inserter(common01, common01.begin()));
    std::set<Region_ID> common012;
    std::set_intersection(common01.begin(), common01.end(),
                          region_ids[2].begin(), region_ids[2].end(),
                          std::inserter(common012, common012.begin()));

    std::cout << " ***" << std::endl;
    std::cout << "Considering\n"
              << "\t" << wrapper.triangulation().point(f.first, (f.second + 1)%4) << "\n"
              << "\t" << wrapper.triangulation().point(f.first, (f.second + 2)%4) << "\n"
              << "\t" << wrapper.triangulation().point(f.first, (f.second + 3)%4) << "\n"
              << "S0, " << region_ids[0].size() << ": ";
    for(const Region_ID& id : region_ids[0])
      std::cout << id << " ";
    std::cout << "\nS1, " << region_ids[1].size() << ": ";
    for(const Region_ID& id : region_ids[1])
      std::cout << id << " ";
    std::cout << "\nS2, " << region_ids[2].size() << ": ";
    for(const Region_ID& id : region_ids[2])
      std::cout << id << " ";
    std::cout << "\ncommon, " << common012.size() << ": ";
    for(const Region_ID& id : common012)
      std::cout << id << " ";
    std::cout << std::endl;

    if(!common012.empty())
      std::cout << "Skipping!" << std::endl;

    // Do not consider the facet if there is a common region ID at all vertices
    // @fixme, this can have false positives since the region is not convex
    return common012.empty();
  }
};

int main(int argc, char** argv)
{
  std::cout.precision(17);

  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh) || !is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;
  const double angle = (argc > 4) ? std::stod(argv[4]) : 30.;

  const CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;

  // Region Growing to partition the mesh
  const double cos_threshold = std::cos(angle * CGAL_PI / 180);

  using Patch_id_pmap = boost::property_map<Mesh, CGAL::face_patch_id_t<int> >::type;
  Patch_id_pmap rfm = get(CGAL::face_patch_id_t<int>(), mesh);
  boost::vector_property_map<CGAL::Epick::Vector_3> normal_map;

  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(
      mesh, rfm,
      CGAL::parameters::cosine_of_maximum_angle(cos_threshold)
                       .region_primitive_map(normal_map)
                       .maximum_distance(1e10) // tmp
                       .postprocess_regions(true));

  std::cout << "#Regions = " << nb_regions << std::endl;

  // Dump all regions
  for(std::size_t i=0; i<nb_regions; ++i)
  {
    CGAL::Face_filtered_graph<Mesh> ffg(mesh, i, rfm);
    CGAL::IO::write_polygon_mesh("results/regions/region_" + std::to_string(i) + ".off", ffg, CGAL::parameters::stream_precision(17));
  }

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  Filtering_visitor visitor(mesh, rfm);
  Mesh wrap;

  CGAL::alpha_wrap_3(mesh, alpha, offset,
                     -1. /*tolerance*/, -1. /*max length*/, false /*parsimonious*/, false /*dump*/,
                     wrap,
                     CGAL::parameters::visitor(visitor));

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  std::string input_name = std::string(filename);
  input_name = input_name.substr(input_name.find_last_of("/") + 1, input_name.length() - 1);
  input_name = input_name.substr(0, input_name.find_last_of("."));
  std::string output_name = input_name
                            + "_" + std::to_string(static_cast<int>(relative_alpha))
                            + "_" + std::to_string(static_cast<int>(relative_offset)) + ".off";
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
