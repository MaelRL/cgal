// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// the canvas used is a starset
#define ANISO_NO_CONSISTENCY
#define ANISO_OUTPUT_WIP

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Starset.h>
#include <CGAL/IO/Star_set_IO.h>
#include <CGAL/helpers/starset_merger.h>

#include <CGAL/Anisotropic_mesher_3.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Anisotropic_tet_mesher_3.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/External_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>

#include <Domain/Constrain_surface_3_ellipse.h>
#include <Domain/Constrain_surface_3_cube.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/Campen_starset_canvas.h>
#include <CGAL/Canvas/Campen_starset_point.h>
#include <CGAL/Canvas/mesher.h>

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

  typedef CGAL::Anisotropic_mesh_3::Metric_field<K>            Metric_field;
  typedef Constrain_surface_3<K>                               Constrain_surface;
  typedef Criteria_base<K>                                     Criteria;

  typedef Starset_with_info<K, Constrain_surface, Metric_field, Criteria>
                                                               Star_set;
  typedef Campen_starset_point<K, Constrain_surface, Metric_field, Criteria>
                                                               Campen_canvas_point;
  typedef Canvas<K, Campen_canvas_point, Metric_field>         Base_canvas;
  typedef Campen_starset_canvas<K, Constrain_surface, Metric_field, Criteria>
                                                               Canvas;
  typedef Canvas_mesher<Base_canvas>                           Canvas_mesher;

  // ------------------------- end of typedefs --------------------------
  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  double duration;
  std::clock_t start = std::clock();
  std::srand(0);

  //geometry
    FT a = 2.;
    FT b = 2.;
    FT c = 2.;

  //metric field
    FT epsilon = 1e-6;
    FT en_factor = 0.999;

  //facet criteria
    FT approx = 1;
    FT f_r0 = 0.1;
    FT f_rho0 = 3.0;

  //cell criteria
    FT c_r0 = 0.1;
    FT c_rho0 = 3.0;
    FT sliverity = 0.3;

  //misc
    FT gamma = 1.5;
    FT beta = 2.5;
    FT delta = 0.3;
    int max_times_to_try = 60;
    int nb = 20;

    Criteria* criteria = new Criteria(approx, f_r0, f_rho0,
                                      c_r0, c_rho0, sliverity,
                                      gamma, beta, delta, nb,
                                      max_times_to_try);

  //----------- pick a domain! ----------
//  Constrain_surface_3_ellipse<K>* pdomain =
//      new Constrain_surface_3_ellipse<K>(a, b, c);
  Constrain_surface_3_free_cube<K>* pdomain =
      new Constrain_surface_3_free_cube<K>(0., 0., 0., 3., 3., 3.);
//  Constrain_surface_3_polyhedral<K>* pdomain = new Constrain_surface_3_polyhedral<K>("../../data/Anisotropy_CMP/3DSurface/Fandisk.off");

  //----------- pick a metric field! ----
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();
//  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>(1.,1.,1.);
//  External_metric_field<K>* metric_field = new External_metric_field<K>(*pdomain, "../../data/Anisotropy_CMP/Ours_Results/Fandisk_Ours/global_smooth/modified_fandisk_metric.txt");
//  Implicit_curvature_metric_field<K>* metric_field = new Implicit_curvature_metric_field<K>(*pdomain);
//  Hyperbolic_shock_metric_field<K>* metric_field = new Hyperbolic_shock_metric_field<K>(0.6);
//  Polyhedral_curvature_metric_field<K>* metric_field = new Polyhedral_curvature_metric_field<K>(*pdomain);

  Star_set starset(pdomain, metric_field, criteria);

  //----------- pick a starset mesher! -----------





#ifdef ANISO_GEO_USE_SURFACE_ONLY
  Anisotropic_surface_mesher_3<K> ss_mesher(starset, pdomain, criteria, metric_field);
#else
  Anisotropic_mesher_3<K> ss_mesher(starset, pdomain, criteria, metric_field);
//  Anisotropic_tet_mesher_3<K> ss_mesher(starset, pdomain, criteria, metric_field);
#endif


  //------------ or simply read a dump!-------------
  // read_dump(starset, "dump.txt");
//  merge_8_small_cubes(starset);

//  ss_mesher.refine_mesh();
//  std::ofstream out_dump("dump.txt");
//  dump(starset, out_dump);

//  std::ofstream out("starset.mesh");
//  output_medit(starset, out, false/*with inconsistencies*/);
//  exit(0);

  // ---------------- aniso geo starts here -------------------------------

  // select the canvas
  const std::string canvas_str = "1D_starset";

  // select the input seeds
  std::size_t max_seeds_n = 5;
  const std::string seeds_str = "input.mesh";

  Canvas canvas(canvas_str, seeds_str, max_seeds_n, &starset);
  canvas.initialize();
  canvas.paint();
  canvas.fix_distance_and_seed_id_at_corners();

  // Refinement
  std::size_t n_refine = 0;
  Canvas_mesher mesher(canvas, n_refine);
  mesher.refine();

  canvas.output_canvas_data_and_primal(canvas_str + "_tr");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}
