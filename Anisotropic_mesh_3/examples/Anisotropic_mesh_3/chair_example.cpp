#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_chair.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Timer Timer;

std::string output_filename(const double& a,
                            const double& b,
                            const double& k)
{
  std::ostringstream oss;
  oss << "chair_" << a << "_" << b << "_" << k << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  Timer timer;
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT a = (argc > n) ? atof(argv[n++]) : 0.8;
  K::FT b = (argc > n) ? atof(argv[n++]) : 0.4;
  K::FT k = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 1.0;
  K::FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  K::FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  K::FT f_r0 = (argc > n) ? atof(argv[n++]) : 1.0;

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

  timer.start();
  Constrain_surface_3_chair<K>* pdomain =
      new Constrain_surface_3_chair<K>(a, b, k);
  Implicit_curvature_metric_field<K>* metric_field =
      new Implicit_curvature_metric_field<K>(*pdomain, epsilon);

  typedef Anisotropic_surface_mesher_3<K> Mesher;

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);

  Mesher mesher(starset);
  mesher.refine_all();

  timer.stop();
  std::cerr << timer.time() << " seconds" << std::endl;

  std::string file = output_filename(a, b, k);
  mesher.output(file);

  delete pdomain;
  return 0;
}
