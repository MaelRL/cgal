// This is based on the paper :
// Konukoglu et al., A recursive anisotropic fast marching approach to reaction diffusion equation

// The canvas used is an orthogonal grid

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/Konukoglu_canvas.h>
#include <CGAL/Canvas/canvas_mesher.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <cstddef>
#include <ctime>
#include <iostream>
#include <string>

using namespace CGAL::Anisotropic_mesh_3;

int main(int, char**)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
  typedef CGAL::Exact_predicates_exact_constructions_kernel        KExact;

  typedef typename K::FT                                           FT;
  typedef typename K::Point_3                                      Point_3;

  typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
  //typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

  typedef Konukoglu_canvas<K, KExact, MF>                          Canvas;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  std::clock_t start = std::clock();
  double duration;

  MF mf = new Euclidean_metric_field<K>(1., 1., 3.);
//  MF mf = new Custom_metric_field<K>();

  const std::string canvas_str = "grid";
  const std::string seeds_str = "base_mesh.mesh";
  std::size_t max_seeds_n = 1;
  std::size_t n_refine = 0;

  // canvas geometry
  Point_3 center(1., 1., 1.);
  const FT canvas_side = 2.;
  FT points_per_side = 50.; // number of points per side of the canvas

  // how big of an update the new value needs to be
  FT recursive_tolerance = 1e-5;

  Canvas canvas(canvas_str, seeds_str,
                center, canvas_side, points_per_side,
                max_seeds_n,
                mf,
                recursive_tolerance);

  canvas.initialize();
  canvas.paint();
  canvas.output_canvas_data_and_primal(canvas_str + "_tr");
  canvas.cosphericity_sweeper();

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "Total: " << duration << std::endl;
}
