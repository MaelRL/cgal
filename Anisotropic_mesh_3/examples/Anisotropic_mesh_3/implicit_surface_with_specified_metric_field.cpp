
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <Domain/Constrain_surface_3_sphere.h>
//#include <Domain/Constrain_surface_3_ellipse.h>

//#include <Metric_field/Arctan_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
//#include <CGAL/Euclidean_metric_field.h>

//#define CGAL_PROFILE

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef Anisotropic_surface_mesher_3<K> Mesher;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);

  K::FT r0 = (argc > 1) ? atof(argv[1]) : 1.0/*default*/;
  K::FT gamma0 = (argc > 2) ? atof(argv[2]) : 1.3;
  K::FT epsilon = (argc > 3) ? atof(argv[3]) : 0.125;

  Criteria_base<K> criteria(3.0, //radius_edge_ratio_
                            0.2, //sliverity_
                            r0, //circumradius_ 0.1
                            gamma0, //distortion_ 1.3
                            2.5, //beta_
                            0.3);//delta_

  Constrain_surface_3_sphere<K>* pdomain = new Constrain_surface_3_sphere<K>();

  //Constrain_surface_3_ellipse<K>* pdomain
  //  = new Constrain_surface_3_ellipse<K>(2., 2., 1.);

  //Arctan_metric_field<K> metric_field(5.0, 0.2, 0.1);
  Hyperbolic_shock_metric_field<K> metric_field(0.6, epsilon);
  /*
  Euclidean_metric_field<K> metric_field;
  if(argc > 1)
  {
    metric_field.a() = atof(argv[1]);
    metric_field.b() = atof(argv[2]);
    metric_field.c() = atof(argv[3]);
  }
  */

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain);
  Mesher mesher(starset);
  mesher.refine_all();
  //mesher.output(std::string("sphere_arctan.off"));
  mesher.output(std::string("sphere_shock.off"));

  delete pdomain;
  return 0;
}
