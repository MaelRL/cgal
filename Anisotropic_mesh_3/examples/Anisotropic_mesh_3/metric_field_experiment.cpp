#include <CGAL/Default_configuration.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_torus.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>



using namespace CGAL::Anisotropic_mesh_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;


int main(int argc, char* argv[])
{
  std::ofstream fx("experiment.txt");

#ifdef ANISO_USE_EIGEN
  std::cout << "Use Eigen" << std::endl;
#else
  std::cout << "Don't use Eigen" << std::endl;
#endif

  CGAL::default_random = CGAL::Random(0);

  K::FT R = (argc > 1) ? atof(argv[1]) : 10.;
  K::FT r = (argc > 2) ? atof(argv[2]) : 1.;
  K::FT epsilon = (argc > 3) ? atof(argv[3]) : (r / R);

  fx << "Torus :" << std::endl;
  fx << "\tR = " << R << std::endl;
  fx << "\tr = " << r << std::endl;
  fx << "\teps = " << epsilon << std::endl;

  Constrain_surface_3_torus<K>* pdomain
    = new Constrain_surface_3_torus<K>(R, r);

  Implicit_curvature_metric_field<K> metric_field(*pdomain, epsilon);

  //u turns around a small circle
  //v turns around a large circle
  
  fx << std::endl << "On small cirle [v = 45 deg] : " << std::endl;
  fx << "u\t" << "Distortions\t" << "Approximations" << std::endl;
  K::FT u = 0.;
  K::FT v = .25*CGAL_PI; //constant
  K::Point_3 p = pdomain->point_on_surface(u, v);
  Implicit_curvature_metric_field<K>::Metric m = metric_field.compute_metric(p);

  int n = 50; //n+1 points on half-the-small-circle
  for(int i = 1; i < n+1; ++i)
  {
    u = i * CGAL_PI / ((double)n);
    K::Point_3 pi = pdomain->point_on_surface(u, v);
    fx << (i * 180. / ((double)n)) << "\t";

    Implicit_curvature_metric_field<K>::Metric mi = metric_field.compute_metric(pi);
    K::FT dis = m.compute_distortion(mi);
    fx << dis << "\t\t";

    K::FT approx = pdomain->compute_sq_approximation(CGAL::midpoint(p, pi));
    fx << std::sqrt(approx) << std::endl;
  }

  fx << std::endl << "On large circle [u = 0 deg] : " << std::endl;
  fx << "v\t" << "Distortions\t" << "Approximations" << std::endl;
  u = 0.;//constant
  v = 0.;
  p = pdomain->point_on_surface(u, v);
  Implicit_curvature_metric_field<K>::Metric m2 = metric_field.compute_metric(p);
  for(int i = 1; i < n+1; ++i)
  {
    v = i * CGAL_PI / ((double)n);
    K::Point_3 pi = pdomain->point_on_surface(u, v);
    fx << (i * 180. / ((double)n)) << "\t";

    Implicit_curvature_metric_field<K>::Metric mi = metric_field.compute_metric(pi);
    K::FT dis = m2.compute_distortion(mi);
    fx << dis << "\t\t";

    K::FT approx = pdomain->compute_sq_approximation(CGAL::midpoint(p, pi));
    fx << std::sqrt(approx) << std::endl;
  }

  fx << std::endl << "On large circle [u = 90 deg] : " << std::endl;
  fx << "v\t" << "Distortions\t" << "Approximations" << std::endl;
  u = 0.5*CGAL_PI;//constant
  v = 0.;
  p = pdomain->point_on_surface(u, v);
  Implicit_curvature_metric_field<K>::Metric m3 = metric_field.compute_metric(p);
  for(int i = 1; i < n+1; ++i)
  {
    v = i * CGAL_PI / ((double)n);
    K::Point_3 pi = pdomain->point_on_surface(u, v);
    fx << (i * 180. / ((double)n)) << "\t";

    Implicit_curvature_metric_field<K>::Metric mi = metric_field.compute_metric(pi);
    K::FT dis = m3.compute_distortion(mi);
    fx << dis << "\t\t";

    K::FT approx = pdomain->compute_sq_approximation(CGAL::midpoint(p, pi));
    fx << std::sqrt(approx) << std::endl;
  }

  fx << std::endl << "On large circle [u = 180 deg] : " << std::endl;
  fx << "v\t" << "Distortions\t" << "Approximations" << std::endl;
  u = CGAL_PI;//constant;
  v = 0.;
  p = pdomain->point_on_surface(u, v);
  Implicit_curvature_metric_field<K>::Metric m4 = metric_field.compute_metric(p);
  for(int i = 1; i < n+1; ++i)
  {
    v = i * CGAL_PI / ((double)n);
    K::Point_3 pi = pdomain->point_on_surface(u, v);
    fx << (i * 180. / ((double)n)) << "\t";

    Implicit_curvature_metric_field<K>::Metric mi = metric_field.compute_metric(pi);
    K::FT dis = m4.compute_distortion(mi);
    fx << dis << "\t\t";

    K::FT approx = pdomain->compute_sq_approximation(CGAL::midpoint(p, pi));
    fx << std::sqrt(approx) << std::endl;
  }


  fx.close();
  delete pdomain;
  return 0;
}
