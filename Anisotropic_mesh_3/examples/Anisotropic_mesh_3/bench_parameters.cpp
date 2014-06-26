#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_torus.h>

#include "refinement_condition_is_between_planes.h"

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Timer Timer;

int main(int argc, char* argv[])
{
  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  std::ofstream fout("bench_beta_delta.txt");
  fout << "beta\t" << "delta\t" << "nbv\t" << "time(sec)" << std::endl;
  std::cout << "beta\t" << "delta\t" << "nbv\t" << "time(sec)" << std::endl;

  int n = 1;

//geometry
  K::FT R = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT r = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 0.1;
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
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 10;

  K::FT betas[] = {2., 3., 4., 5.};
  K::FT deltas[] = {0.1, 0.3, 0.5};

  Constrain_surface_3_torus<K>* pdomain = new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* metric_field =
      new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);

//  K::FT y = 1.1*(R + r);
//  K::Plane_3 plane1(0., 1., 0., -y);
//  K::Plane_3 plane2(0., 1., 0.,  y);
//  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
//  RCondition condition(plane1, plane2, true/*x>0 only*/);

  for(int i = 0; i < 4; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      std::cout << betas[i] << "\t" << deltas[j];
      Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                        c_r0, c_rho0, sliverity,
                                                        gamma0, betas[i], deltas[j], nb,
                                                        max_times_to_try);

      timer.reset();
      timer.start();

      Starset<K> starset;
      Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

      mesher.refine_mesh();
      timer.stop();

      fout << betas[i] << "\t"
           << deltas[j] << "\t"
           << starset.size() << "\t"
           <<  timer.time() << std::endl;

      delete criteria;
    }
  }

  delete metric_field;
  delete pdomain;

  return 0;
}
