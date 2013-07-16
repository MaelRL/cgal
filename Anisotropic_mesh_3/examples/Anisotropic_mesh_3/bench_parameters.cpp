#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_torus.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>

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

  K::FT R = 10.;
  K::FT r =  1.;
  K::FT epsilon = 0.1;
  K::FT r0 = 0.1;
  K::FT gamma0 = 1.5;
  int nb_initial_pts = 10;

  K::FT betas[] = {2., 3., 4., 5.};
  K::FT deltas[] = {0.1, 0.3, 0.5};

  Constrain_surface_3_torus<K>* pdomain =
    new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon);

  K::FT y = 1.1*(R + r);
  K::Plane_3 plane1(0., 1., 0., -y);
  K::Plane_3 plane2(0., 1., 0.,  y);
  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
  RCondition condition(plane1, plane2, true/*x>0 only*/);

  for(int i = 0; i < 4; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      std::cout << betas[i] << "\t" << deltas[j];
      Criteria_base<K>* criteria = new Criteria_base<K>(3.0, //radius_edge_ratio_
                                                        0.2,    //sliverity_
                                                        r0,     //circumradius_ 0.1
                                                        gamma0, //distortion_ 1.3
                                                        betas[i],   //beta_ 2.5
                                                        deltas[j]); //delta_ 0.3
      timer.reset();
      timer.start();
      Surface_star_set_3<K, RCondition> starset(criteria, metric_field, pdomain, nb_initial_pts, condition);
      //Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb_initial_pts);
      starset.refine_all();
      timer.stop();

      fout << betas[i] << "\t"
           << deltas[j] << "\t"
           << starset.number_of_stars() << "\t"
           <<  timer.time() << std::endl;
    }
  }

  delete pdomain;
  return 0;
}
