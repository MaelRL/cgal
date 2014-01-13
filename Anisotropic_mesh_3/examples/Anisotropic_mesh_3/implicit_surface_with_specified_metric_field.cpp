
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <Domain/Constrain_surface_3_sphere.h>
//#include <Domain/Constrain_surface_3_ellipse.h>

//#include <Metric_field/Arctan_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
//#include <CGAL/Euclidean_metric_field.h>

//#define CGAL_PROFILE

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef Anisotropic_surface_mesher_3<K> Mesher;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 0.1;
  K::FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  K::FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  K::FT f_r0 = (argc > n) ? atof(argv[n++]) : 0.1;

//cell criteria
  K::FT sliverity = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT c_rho0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT c_r0 = (argc > n) ? atof(argv[n++]) : 0.;
  bool c_consistency = false;

//misc
  K::FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  K::FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, gamma0, f_rho0, f_r0,
                                                    sliverity, c_rho0, c_r0,
                                                    c_consistency, beta, delta,
                                                    max_times_to_try);

  Constrain_surface_3_sphere<K>* pdomain = new Constrain_surface_3_sphere<K>();

  //Constrain_surface_3_ellipse<K>* pdomain
  //  = new Constrain_surface_3_ellipse<K>(2., 2., 1.);

  //Arctan_metric_field<K>* metric_field = new Arctan_metric_field<K>(5.0, 0.2, 0.1);
  Hyperbolic_shock_metric_field<K>* metric_field =
    new Hyperbolic_shock_metric_field<K>(0.6, epsilon);

  /*
  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();
  if(argc > 1)
  {
    metric_field.a() = atof(argv[1]);
    metric_field.b() = atof(argv[2]);
    metric_field.c() = atof(argv[3]);
  }
  */

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  Mesher mesher(starset);
  mesher.refine_all();
  //mesher.output(std::string("sphere_arctan.off"));
  mesher.output(std::string("sphere_shock.off"));

  delete pdomain;
  return 0;
}
