#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_cylinder.h>
#include <Domain/Constrain_surface_3_torus.h>
#include <Metric_field/Cylinder_metric_field.h>

#include "refinement_condition_is_between_planes.h"

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT R = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT r = (argc > n) ? atof(argv[n++]) : 1.;

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

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta, nb,
                                                    max_times_to_try);

//  K::FT y = (argc > n) ? atof(argv[n++]) : (R + r + 1);
//  int xcondition = (argc > n) ? atoi(argv[n++]) : -1;//default : no condition on x
//  K::Plane_3 plane1(0., 1., 0., -y);
//  K::Plane_3 plane2(0., 1., 0., y);
//  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
//  RCondition condition(plane1, plane2, (xcondition == 1));

  Constrain_surface_3_torus<K>* pdomain = new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* curvature_mf = new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  Euclidean_metric_field<K>* euclid_mf = new Euclidean_metric_field<K>();

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> curvature_mesher(starset, pdomain, criteria, curvature_mf);
  curvature_mesher.refine_mesh();
  output_surface_off(starset, "mesh_curvature.off");

  starset.clear();
  Anisotropic_surface_mesher_3<K> euclidean_mesher(starset, pdomain, criteria, euclid_mf);
  euclidean_mesher.refine_mesh();
  output_surface_off(starset, "mesh_euclidean.off");

  delete curvature_mf;
  delete euclid_mf;
  delete pdomain;
  delete criteria;

  return 0;
}
