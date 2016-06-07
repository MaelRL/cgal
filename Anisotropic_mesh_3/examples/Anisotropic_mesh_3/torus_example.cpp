#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Torus_metric_field.h>
#include <Domain/Constrain_surface_3_torus.h>

#include "refinement_condition_is_between_planes.h"

#include <string>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Timer Timer;

std::string output_filename(const double& R,
                            const double& r)
{
  std::ostringstream oss;
  oss << "torus_" << R << "_" << r << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT R = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT r = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;
  K::FT en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

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

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta, nb,
                                                    max_times_to_try);

  Constrain_surface_3_torus<K>* pdomain = new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* metric_field = new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);
  //Torus_metric_field<K>* metric_field = new Torus_metric_field<K>(R, r, epsilon, en_factor);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();


//  K::FT y = (argc > n) ? atof(argv[n++]) : (0.5*(R + r));
//  int xcondition = (argc > n) ? atoi(argv[n++]) : -1;//default : no condition on x
//  K::Plane_3 plane1(0., 1., 0., 0.);
//  double angle = M_PI/16.;
//  K::Plane_3 plane2(K::Point_3(0.,0.,0.), K::Point_3(0.,0.,1.), K::Point_3(std::cos(angle), std::sin(angle),0.));
//  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
//  RCondition condition(plane1, plane2, 1);//(xcondition == 1));

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  mesher.refine_mesh();

  output_surface_medit(starset, "torus_example.mesh");
  std::string file = output_filename(R, r);
  output_surface_off(starset, file.c_str());

  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}
