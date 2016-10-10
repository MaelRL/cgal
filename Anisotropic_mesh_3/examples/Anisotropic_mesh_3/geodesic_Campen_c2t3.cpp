// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// the canvas used is a c2t3

// #define IMPLICIT_DOMAIN

#include <CGAL/Mesh_3/global_parameters.h>

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/Campen_c2t3_canvas.h>
#include <CGAL/Canvas/Campen_c2t3_point.h>
#include <CGAL/Canvas/mesher.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>
#include <CGAL/Canvas/optimizer_cvt.h>
#include <CGAL/Canvas/c2t3_canvas_generator.h>
#include <CGAL/Canvas/geodesic_drawing_helper.h>
#include <CGAL/Canvas/optimizer_cvt_surf.h>

#include <Domain/Constrain_surface_3_ellipse.h>
#include <Domain/Constrain_surface_3_chair.h>
#include <Domain/Constrain_surface_3_cube.h>
#include <Domain/Constrain_surface_3_sphere.h>
#include <CGAL/Constrain_surface_3_polyhedral.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>

#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
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
  typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;


  typedef int                          Vertex_Info; // index of the canvas point
  typedef int                          Cell_Info; // index of the subdomain

  // Define the domain type
#ifdef IMPLICIT_DOMAIN
    // implicit
 typedef typename K::FT                               FT;
 typedef typename K::Point_3                          Point_3;
 typedef FT (Function)(const Point_3&);
 typedef CGAL::Implicit_mesh_domain_3<Function, K>    Implicit_Mesh_domain;
 typedef CGAL::Mesh_domain_with_polyline_features_3<Implicit_Mesh_domain> 
                                                      Mesh_domain;
#else
  // polyhedral
  typedef CGAL::Polyhedron_3<K>                           Polyhedron;
//  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>   Mesh_domain;
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
#endif

  // Define a C3T3 whose underlying regular triangulation has an info member on
  // both the vertices and the cells
  typedef typename CGAL::details::Mesh_geom_traits_generator<K>::type       Gt;

  typedef CGAL::Triangulation_vertex_base_with_info_3<Vertex_Info, Gt>      Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<Cell_Info, Gt>          Cb;
  typedef CGAL::Regular_triangulation_cell_base_3<Gt, Cb>                   Rcb;
  typedef CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<Gt, Rcb>
                                                                            Rcbwc;
  typedef CGAL::Mesh_vertex_base_3<Gt, Mesh_domain, Vb>                     Mvb;
  typedef CGAL::Mesh_cell_base_3<Gt, Mesh_domain, Rcbwc>                    Mcb;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain,
                                              K,
                                              CGAL::Sequential_tag,
                                              Mvb, Mcb>::type  Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>      C3t3;

  typedef Metric_field<K>                                  MF;
  typedef Campen_canvas_point<K, MF, C3t3>                 Campen_canvas_point;
  typedef Canvas<K, Campen_canvas_point, MF>               Base_canvas;
  typedef Campen_canvas<K, MF, C3t3>                       Canvas;
  typedef Canvas_mesher<Base_canvas>                       Canvas_mesher;
  typedef CVT_surf_optimizer<Canvas, C3t3>                 CVT_surf_optimizer;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  std::streambuf * old = std::cout.rdbuf();
//  std::cout.rdbuf(0);

  double duration;
  std::clock_t start = std::clock();
  std::srand(0);

  const std::string domain_str = "fandisk";
  const std::string domain_addr = "/home/mrouxell/Data/OFF/" + domain_str + ".off";

  //----------- pick a metric field! ------------
//  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();
//  Euclidean_metric_field<K>* metric_field =
//                            new Euclidean_metric_field<K>(1.0, 1.0, 1.0);
  Hyperbolic_shock_metric_field<K>* metric_field =
                            new Hyperbolic_shock_metric_field<K>(0.6);

//  Constrain_surface_3_ellipse<K>* pdomain =
//      new Constrain_surface_3_ellipse<K>(1., 1., 1.);
//  Constrain_surface_3_chair<K>* pdomain =
//                          new Constrain_surface_3_chair<K>(0.8, 0.4, 1.0);
//  Implicit_curvature_metric_field<K>* metric_field =
//                          new Implicit_curvature_metric_field<K>(*pdomain, 0.1);

//  Constrain_surface_3_polyhedral<K>* pdomain =
//                    new Constrain_surface_3_polyhedral<K>(domain_addr.c_str());
//  Polyhedral_curvature_metric_field<K>* metric_field =
//                    new Polyhedral_curvature_metric_field<K>(*pdomain);

  //----------- Generate the canvas! ------------

  const std::string canvas_str = "c2t3_" + domain_str;
  C3t3 c3t3;

  CGAL::generate_canvas<C3t3, MF, Mesh_domain>(c3t3, metric_field);

  // select the input seeds
  std::size_t max_seeds_n = 5;
  const std::string seeds_str = "input_c2t3.mesh";
//  const std::string seeds_str = "c2t3_fertility_tr_primal.mesh";

  Canvas canvas(c3t3, canvas_str, seeds_str, max_seeds_n, metric_field);
  canvas.initialize();
  canvas.paint();

  // Refinement (meshing criteria are in mesher.h atm)
  std::size_t n_refine = 200;
  Canvas_mesher mesher(canvas, n_refine);
  mesher.refine();

  canvas.output_canvas_data_and_primal(canvas_str + "_tr");

  // Optimization
  std::size_t max_opti_iter = 0;
  CVT_surf_optimizer optimizer(canvas, max_opti_iter);
  optimizer.optimize_seeds(canvas_str);

  // output geodesics
//  compute_geodesics(canvas, canvas_str);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
