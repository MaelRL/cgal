// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// the canvas used is a c2t3

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/Campen_c2t3_canvas.h>
#include <CGAL/Canvas/Campen_c2t3_point.h>
#include <CGAL/Canvas/mesher.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>
#include <CGAL/Canvas/optimizer_cvt.h>
#include <CGAL/Canvas/c2t3_canvas_generator.h>

#include <Domain/Mesh_3/surface_3_cube.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

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

  typedef typename K::FT                                       FT;
  typedef typename K::Point_3                                  Point_3;

  typedef int                          Vertex_Info; // index of the canvas point
  typedef int                          Cell_Info; // index of the subdomain

  // Define the domain type
  typedef FT (Function)(const Point_3&);
  typedef CGAL::Implicit_mesh_domain_3<Function, K>                Mesh_domain;
  typedef CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain>  Mesh_domain_WF;

  // Define a C3T3 whose underlying regular triangulation has an info member on
  // both the vertices and the cells
  typedef typename CGAL::details::Mesh_geom_traits_generator<K>::type       Gt;

  typedef CGAL::Triangulation_vertex_base_with_info_3<Vertex_Info, Gt>      Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<Cell_Info, Gt>          Cb;
  typedef CGAL::Regular_triangulation_cell_base_3<Gt, Cb>                   Rcb;
  typedef CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<Gt, Rcb>
                                                                            Rcbwc;
  typedef CGAL::Mesh_vertex_base_3<Gt, Mesh_domain_WF, Vb>                  Mvb;
  typedef CGAL::Mesh_cell_base_3<Gt, Mesh_domain_WF, Rcbwc>                 Mcb;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain_WF,
                                              K,
                                              CGAL::Sequential_tag,
                                              Mvb, Mcb>::type               Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                       C3t3;

  typedef Metric_field<K>                                  MF;
  typedef Campen_canvas_point<K, MF, C3t3>                 Campen_canvas_point;
  typedef Canvas<K, Campen_canvas_point, MF>               Base_canvas;
  typedef Campen_canvas<K, MF, C3t3>                       Canvas;
  typedef Canvas_mesher<Base_canvas>                       Canvas_mesher;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  double duration;
  std::clock_t start = std::clock();
  std::srand(0);

  //----------- pick a metric field! ------------
//  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();
  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>(1.,1.,1.);

  //----------- Generate the canvas! ------------
  const std::string canvas_str = "input_c2t3";
  C3t3 c3t3;
  CGAL::generate_canvas<C3t3, MF, Mesh_domain_WF>(c3t3, metric_field);

  // select the input seeds
  std::size_t max_seeds_n = 1;
  const std::string seeds_str = "test_input.mesh";

  Canvas canvas(c3t3, canvas_str, seeds_str, max_seeds_n, metric_field);
  canvas.initialize();
  canvas.paint();

  // Refinement
  std::size_t n_refine = 0;
  Canvas_mesher mesher(canvas, n_refine);
  mesher.refine();

  canvas.output_canvas_data_and_primal(canvas_str + "_tr");

  // Optimization
  // todo

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
