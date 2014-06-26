#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_cyclide.h>

#include "refinement_condition_is_between_planes.h"
#include <fstream>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

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
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT a = (argc > n) ? atof(argv[n++]) : 7.;
  K::FT b = (argc > n) ? atof(argv[n++]) : std::sqrt(48.);
  K::FT mu = (argc > n) ? atof(argv[n++]) : 3.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;
  K::FT en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 1.0;
  K::FT f_r0 = (argc > n) ? atof(argv[n++]) : 0.1;
  K::FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;

//cell criteria
  K::FT c_r0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT c_rho0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT sliverity = (argc > n) ? atof(argv[n++]) : 0.;

//misc
  K::FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  K::FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  K::FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta, nb,
                                                    max_times_to_try);
  
  Constrain_surface_3_cyclide<K>* pdomain = new Constrain_surface_3_cyclide<K>(a, b, mu);
  Implicit_curvature_metric_field<K>* metric_field =
      new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  mesher.refine_mesh();

  std::string file = output_filename(a, b, mu);
  output_surface_off(starset, file.c_str());

  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}
