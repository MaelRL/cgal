// #define CGAL_AW3_DEBUG_DUMP_EVERY_STEP
#define CGAL_AW3_DEBUG
#define CGAL_AW3_DEBUG_QUEUE
#define CGAL_AW3_DEBUG_STEINER_COMPUTATION

// {
     // #define CGAL_AW3_ORIGINAL_STEINER_CONSTRUCTION
// }
// OR
// {
    // main strategy for sharpening
    #define CGAL_AW3_SHARPEN_WITH_SAMPLING
    // #define CGAL_AW3_SHARPEN_WITH_ITERATIVE_SAMPLING
    // #define CGAL_AW3_SHARPEN_WITH_ITERATIVE_BALL_QUERIES

    // pick one of these
    // #define CGAL_AW3_SAMPLE_WITH_RANDOM_WALKS
    #define CGAL_AW3_SAMPLE_WITH_FOUNTAINS
    // #define CGAL_AW3_SAMPLE_WITH_BARYCENTRIC_SEEDS

    // #define CGAL_AW3_RANDOM_FOUNTAIN_STREAMS
    #define CGAL_AW3_SAMPLE_WITH_FOUNTAINS_RECURSIVE

    // pick one of these (meaningful only when using random walks)
    #define CGAL_AW3_RANDOM_WALK_USE_WALK_ON_SPHERES
    // #define CGAL_AW3_RANDOM_WALK_USE_FIXED_STEP // not yet implemented

    #define CGAL_AW3_ENHANCE_SAMPLES_WITH_QEM_POINTS

    // pick one of these
    // #define CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_QEM_OPTIMAL_POINT
    // #define CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_QEM_MINIMIZER
    #define CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_CANDIDATE_CLOSEST_TO_CC
// }

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <string>
#include <vector>

namespace AW3 = CGAL::Alpha_wraps_3;
namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

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
  const double relative_tolerance = (argc > 4) ? std::stod(argv[4]) : 1000.;
  const bool parsimonious = (argc > 5) ? (std::string(argv[5]) == "true") : false;
  const bool dump = (argc > 6) ? (std::string(argv[6]) == "true") : false;

  const CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  const double tolerance = alpha / relative_tolerance;

  std::cout << "Diagonal, rel_alpha, rel_offset, rel_tol: " << diag_length << ", " << relative_alpha << ", " << relative_offset << ", " << relative_tolerance << std::endl;
  std::cout <<  "actual alpha, actual offset: " << alpha << ", " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  Mesh wrap;
  CGAL::alpha_wrap_3(mesh, alpha, offset, tolerance, parsimonious, dump, wrap);

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << ", " << num_faces(wrap) << std::endl;
  std::cout << "Duration: " << t.time() << std::endl;

  // Save the result
  std::string input_name = std::string(filename);
  input_name = input_name.substr(input_name.find_last_of("/") + 1, input_name.length() - 1);
  input_name = input_name.substr(0, input_name.find_last_of("."));
  std::string output_name = input_name
                            + "_" + std::to_string(static_cast<int>(relative_alpha))
                            + "_" + std::to_string(static_cast<int>(relative_offset))
                            + "_" + std::to_string(static_cast<int>(relative_tolerance)) + ".off";
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
