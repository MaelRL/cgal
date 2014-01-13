#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Euclidean_metric_field.h>

#include <Domain/Constrain_surface_3_cyclide.h>
#include <CGAL/Criteria.h>
#include <CGAL/Surface_star_set_3.h>

#include "refinement_condition_is_between_planes.h"


using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Timer Timer;


std::string output_filename(const double& a,
                            const double& b,
                            const double& mu)
{
  std::ostringstream oss;
  oss << "cyclide_" << a << "_" << b << "_" << mu << ".off";
  return oss.str();
}


int main(int argc, char* argv[])
{
  Timer timer;
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT a = (argc > n) ? atof(argv[n++]) : 7.;
  K::FT b = (argc > n) ? atof(argv[n++]) : std::sqrt(48.);
  K::FT mu = (argc > n) ? atof(argv[n++]) : 3.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 1.0;
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
  
  Constrain_surface_3_cyclide<K>* pdomain = new Constrain_surface_3_cyclide<K>(a, b, mu);

  Implicit_curvature_metric_field<K>* metric_field = new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  std::string file = output_filename(a, b, mu);

  timer.start();
  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  timer.stop();
  
  starset.output(file.c_str()); // .off

  timer.start();
  starset.refine_all();
  timer.stop();

  starset.output_surface_star_set("cyclide_example.mesh");
  starset.output(file.c_str()); // .off

  delete pdomain;
  return 0;
}
