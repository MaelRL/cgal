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

  K::FT a = (argc > 1) ? atof(argv[1]) : 0.8;
  K::FT b = (argc > 2) ? atof(argv[2]) : 0.4;
  K::FT k = (argc > 3) ? atof(argv[3]) : 1.;

  K::FT epsilon = (argc > 4) ? atof(argv[4]) : 1.0;

  K::FT r0 = (argc > 5) ? atof(argv[5]) : 1.0/*default*/;
  K::FT gamma0 = (argc > 6) ? atof(argv[6]) : 1.5;

  int nb = (argc > 7) ? atoi(argv[7]) : 20;
  
  Criteria_base<K>* criteria = new Criteria_base<K>(3.0, //radius_edge_ratio_
                                                    0.2, //sliverity_
                                                    r0, //circumradius_ 0.1
                                                    gamma0, //distortion_ 1.3
                                                    2.5, //beta_
                                                    0.3);//delta_

  timer.start();
  Constrain_surface_3_chair<K>* pdomain
    = new Constrain_surface_3_chair<K>(a, b, k);
  Implicit_curvature_metric_field<K>* metric_field =
     new Implicit_curvature_metric_field<K>(*pdomain, epsilon);

  typedef Anisotropic_surface_mesher_3<K> Mesher;

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);

  std::cout << "Max curvature : " << pdomain->global_max_curvature() << std::endl;
  std::cout << "Min curvature : " << pdomain->global_min_curvature() << std::endl;

  Mesher mesher(starset);
  mesher.refine_all();


  timer.stop();
  std::cerr << timer.time() << " seconds" << std::endl;

  std::string file = output_filename(a, b, k);
  mesher.output(file);

  delete pdomain;
  return 0;
}
