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
    // #define CGAL_AW3_SHARPEN_WITH_DYNAMIC_BALL_QUERIES

    // suboptions for sampling strategies
    // #define CGAL_AW3_SAMPLE_WITH_RANDOM_WALKS
    // #define CGAL_AW3_SAMPLE_WITH_FOUNTAINS
    // #define CGAL_AW3_SAMPLE_WITH_BARYCENTRIC_SEEDS
    #define CGAL_AW3_SAMPLE_WITH_RAY_INTERSECTIONS

    // suboptions for random walks (pick one)
    #define CGAL_AW3_RANDOM_WALK_USE_WALK_ON_SPHERES
    // #define CGAL_AW3_RANDOM_WALK_USE_FIXED_STEP // not yet implemented

    // sub options for fountain sampling
    // #define CGAL_AW3_RANDOM_FOUNTAIN_STREAMS
    // #define CGAL_AW3_SAMPLE_WITH_FOUNTAINS_RECURSIVE

    // #define CGAL_AW3_ENHANCE_SAMPLES_WITH_QEM_POINTS

    // pick one of these
    #define CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_QEM_OPTIMAL_POINT
    // #define CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_QEM_MINIMIZER
    // #define CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_CANDIDATE_CLOSEST_TO_CC

    #define CGAL_AW3_USE_WEIGHTED_INSERTION_FOR_STEINER_POINTS_ON_FEATURES
// }

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh) || !is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input:" << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

  const CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  std::cout << "alpha: " << alpha << ", offset: " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  Mesh wrap;
  CGAL::alpha_wrap_3(mesh, alpha, offset, wrap);

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  const std::string output_name = generate_output_name(filename, relative_alpha, relative_offset);
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
