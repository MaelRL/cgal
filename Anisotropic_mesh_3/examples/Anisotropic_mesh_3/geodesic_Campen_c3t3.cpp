// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// the canvas used is a triangulation

#define ANISO_GEO_FILTER_SEEDS_OUTSIDE_CANVAS_BBOX

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/Campen_c3t3_canvas.h>
#include <CGAL/Canvas/Campen_c3t3_point.h>
#include <CGAL/Canvas/mesher.h>
#include <CGAL/Canvas/c3t3_canvas_generator.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>
#include <CGAL/Canvas/optimizer_cvt.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <iostream>
#include <map>
#include <vector>
#include <set>

using namespace CGAL::Anisotropic_mesh_3;

int main(int, char**)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

  typedef Metric_field<K>                                      MF;

  typedef Campen_canvas_point<K, MF>                           Campen_canvas_point;
  typedef Canvas<K, Campen_canvas_point, MF>                   Base_canvas;
  typedef Campen_canvas<K, MF>                                 Canvas;
  typedef Canvas_mesher<Base_canvas>                           Canvas_mesher;
  typedef CVT_optimizer<Canvas>                                CVT_optimizer;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  double duration;
  std::clock_t start = std::clock();
  std::srand(0);

  //----------- pick a metric field! ----
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();
//  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>(2.,2.,2.);

//  CGAL::generate_canvas<K, MF>(metric_field);
//  exit(0);

  // select the canvas
  const std::string canvas_str = "input_base_mesh";

  // select the input seeds
  std::size_t max_seeds_n = 2005;
  const std::string seeds_str = "ref_2005_primal.mesh";

  Canvas canvas(canvas_str, seeds_str, max_seeds_n, metric_field);
  canvas.initialize();
  canvas.paint();

  // Refinement
  std::size_t n_refine = 0;
  Canvas_mesher mesher(canvas, n_refine);
  mesher.refine();

  canvas.output_canvas_data_and_primal(canvas_str + "_tr");

  // Optimization
  std::size_t max_opti_iter = 0;
  CVT_optimizer optimizer(canvas, max_opti_iter);
  optimizer.optimize_seeds(canvas_str);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
