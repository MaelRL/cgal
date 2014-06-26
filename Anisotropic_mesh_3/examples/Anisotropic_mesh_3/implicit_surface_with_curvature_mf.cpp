#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <Domain/Constrain_surface_3_ellipse.h>
#include <Domain/Constrain_surface_3_sphere.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Metric_field/Spherical_metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef Anisotropic_surface_mesher_3<K> Mesher;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 0.1;
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
  const double en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta,
                                                    max_times_to_try);

  Constrain_surface_3_sphere<K>* pdomain = new Constrain_surface_3_sphere<K>();
  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);
  //Spherical_metric_field<K>* metric_field = new Spherical_metric_field<K>(epsilon);

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  mesher.refine_mesh();
  output_surface_off(starset, "sphere_curvature_surface.off");
  
  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}
