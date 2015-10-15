// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// the canvas used is a triangulation

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/Campen_triangulation_canvas.h>
#include <CGAL/Canvas/Campen_triangulation_point.h>
#include <CGAL/Canvas/mesher.h>
#include <CGAL/Canvas/canvas_mesh_3.h>
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

//  typedef Euclidean_metric_field<K>*                           MF;
  typedef Custom_metric_field<K>*                              MF;

  typedef Campen_canvas_point<K, MF>                           Campen_canvas_point;
  typedef Canvas<K, Campen_canvas_point, MF>                   Base_canvas;
  typedef Campen_canvas<K, MF>                                 Canvas;
  typedef Canvas_mesher<Base_canvas>                           Canvas_mesher;
  typedef CVT_optimizer<Canvas>                                CVT_optimizer;

  std::cout.precision(17);
  std::freopen("log.txt", "w", stdout);

  double duration;
  std::clock_t start = std::clock();
  std::srand(0);

//  MF mf = new Euclidean_metric_field<K>(1., 1., 3.);
  MF mf = new Custom_metric_field<K>();

  // select the canvas
  const std::string canvas_str = "dense_base_mesh";
  // select the input seeds
  std::size_t max_seeds_n = 150;
  const std::string seeds_str = "dense_base_mesh_tr_primal.mesh";

//  CGAL::generate_canvas<K>();
//  exit(0);

  Canvas canvas(canvas_str, seeds_str, max_seeds_n, mf);
  canvas.initialize();
  canvas.paint();
  canvas.debug();

  std::size_t n_refine = 0;
  Canvas_mesher mesher(canvas, n_refine);
  mesher.refine();

  canvas.output_canvas_data_and_primal(canvas_str + "_tr");
  canvas.output_canvas(canvas_str + "_tr");

  // Optimization
  std::size_t max_opti_iter = 1000;
  CVT_optimizer optimizer(canvas, max_opti_iter);
  optimizer.optimize_seeds(canvas_str);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
