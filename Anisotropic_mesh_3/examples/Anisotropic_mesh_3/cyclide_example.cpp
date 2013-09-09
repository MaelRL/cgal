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
  std::ofstream fx("cyclide_info.txt");

  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  if(argc == 2)
  {
    std::cout << "cyclide_example parameters :" << std::endl;
    std::cout << " r0=1 gamma0=1.5 rho0=3 approx=0";
    std::cout << " nb=20 beta=2.5 delta=0.3 gamma1=2 y=5.5 xcondition=-1 " << std::endl;
    return 0;
  }

  K::FT a = (argc > 1) ? atof(argv[1]) : 7.;
  K::FT b = (argc > 2) ? atof(argv[2]) : std::sqrt(48);
  K::FT mu = (argc > 3) ? atof(argv[3]) : 3.;

  K::FT epsilon = (argc > 4) ? atof(argv[4]) : 0.1;
  K::FT r0 = (argc > 5) ? atof(argv[5]) : 0.;
  K::FT gamma0 = (argc > 6) ? atof(argv[6]) : 0.;
  K::FT rho0 = (argc > 7) ? atof(argv[7]) : 0.;
  K::FT approx = (argc > 8) ? atoi(argv[8]) : 0.1;

  int nb = (argc > 9) ? atof(argv[9]) : 20.;

  K::FT beta = (argc > 10) ? atof(argv[10]) : 2.5;
  K::FT delta = (argc > 11) ? atof(argv[11]) : 0.;

  Criteria_base<K>* criteria = new Criteria_base<K>(rho0, //radius_edge_ratio_
                                                    0.2,    //sliverity_
                                                    r0,     //circumradius_ 0.1
                                                    gamma0, //distortion_ 1.3
                                                    beta,   //beta_ 2.5
                                                    delta,  //delta_ 0.3
                                                    60,     //max_times_to_try_in_picking_region_
                                                    approx); //approximation

  fx << "Cyclide :" << std::endl;
  fx << "\ta = " << a << std::endl;
  fx << "\tb = " << b << std::endl;
  fx << "\tmu = " << mu << std::endl;
  fx << "\teps = " << epsilon << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;
  fx << "\trho0 = " << rho0 << std::endl;
  fx << "\tapprox = " << approx << std::endl;
  fx << "\tbeta = " << beta << std::endl;
  fx << "\tdelta = " << delta << std::endl;
  
  Constrain_surface_3_cyclide<K>* pdomain = new Constrain_surface_3_cyclide<K>(a, b, mu);

  Implicit_curvature_metric_field<K>* metric_field = new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  std::string file = output_filename(a, b, mu);

  fx << std::endl << "nbV" << "\t" << "time" << std::endl;
  timer.start();
  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  timer.stop();
  fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  
  starset.output(file.c_str());//off

  timer.start();
  starset.refine_all();
  timer.stop();
  fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;

  starset.output_surface_star_set("cyclide_example.mesh");//mesh
  starset.output(file.c_str());//off

  delete pdomain;
  return 0;
}
