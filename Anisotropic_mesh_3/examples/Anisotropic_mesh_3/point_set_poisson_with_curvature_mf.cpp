
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <Domain/Constrain_surface_3_pointset.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>


using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef Anisotropic_surface_mesher_3<K> Mesher;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
    
  Constrain_surface_3_pointset<K>* pdomain
    = new Constrain_surface_3_pointset<K>(argv[1]);

  int n = 2;
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
  const double en_factor = (argc > n) ? atof(argv[n++]) : 1.;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, gamma0, f_rho0, f_r0,
                                                    sliverity, c_rho0, c_r0,
                                                    c_consistency, beta, delta,
                                                    max_times_to_try);

  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  Mesher mesher(starset);
  mesher.refine_all();
  mesher.output("pointset_curvature_surface.off");
  
  delete pdomain;
  return 0;
}
