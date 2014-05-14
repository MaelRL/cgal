#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <CGAL/Anisotropic_mesher_3.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Anisotropic_tet_mesher_3.h>

#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
#include <Metric_field/External_metric_field.h>

#include <Domain/Constrain_surface_3_ellipse.h>
#include <Domain/Constrain_surface_3_cube.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>

#include <CGAL/IO/Star_set_output.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef CGAL::Timer                                          Timer;
typedef Stretched_Delaunay_3<K>                              Star;
typedef std::vector<Star*>                                   Star_vector;

int main(int argc, char** argv)
{
  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  int n = 1;

//geometry
  FT a = (argc > n) ? atof(argv[n++]) : 1.;
  FT b = (argc > n) ? atof(argv[n++]) : 1.;
  FT c = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;
  FT en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

//facet criteria
  FT approx = (argc > n) ? atof(argv[n++]) : 1e10;
  FT gamma = (argc > n) ? atof(argv[n++]) : 1.5;
  FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  FT f_r0 = (argc > n) ? atof(argv[n++]) : 1.0;

//cell criteria
  FT sliverity = (argc > n) ? atof(argv[n++]) : 0.2;
  FT c_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  FT c_r0 = (argc > n) ? atof(argv[n++]) : 1.0;

//misc
  FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma, beta, delta, nb,
                                                    max_times_to_try);

  timer.start();

  //----------- pick a domain! ----------
//  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);
  Constrain_surface_3_free_cube<K>* pdomain = new Constrain_surface_3_free_cube<K>(0.,0.,0.,1.,1.,1.);
//  Constrain_surface_3_polyhedral<K>* pdomain = new Constrain_surface_3_polyhedral<K>("../../data/Anisotropy_CMP/3DSurface/Fertility.off");

  //----------- pick a metric field! ----
//  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();
//  Hyperbolic_shock_metric_field<K>* metric_field = new Hyperbolic_shock_metric_field<K>(0.6);
  External_metric_field<K>* metric_field = new External_metric_field<K>(*pdomain, "../../data/Anisotropy_CMP/3DSurface/Fertility_Metric.txt");
//  Polyhedral_curvature_metric_field<K>* metric_field = new Polyhedral_curvature_metric_field<K>(pdomain);

  Starset<K> starset;

  //----------- pick a mesher! -----------
  Anisotropic_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);
//  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);
//  Anisotropic_tet_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  double elapsed_time = mesher.refine_mesh();
  std::cout << "elapsed time: " << elapsed_time << std::endl;
  mesher.report();

  std::ofstream out("bambimboum.mesh");
  output_medit(starset, out);
  std::ofstream out_facet("bambimboum_surf.mesh");
  output_surface_medit(starset, out_facet);
  std::ofstream out_off("bambimboum.off");
  output_off(starset, out_off);
  std::ofstream out_surface_off("bambimboum_surf.off");
  output_surface_off(starset, out_surface_off);

  delete criteria;
  delete pdomain;
  delete metric_field;

  return 0;
}

