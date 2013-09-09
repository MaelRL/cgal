
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

  Criteria_base<K>* criteria = new Criteria_base<K>(0., //radius_edge_ratio_
                                                    0.2, //sliverity_
                                                    0., //circumradius_
                                                    0., //distortion_
                                                    2.5, //beta_
                                                    0.3, //delta_
                                                    60,  //max tries in picking region
                                                    0.); //approximation

  const double epsilon = (argc > 1) ? atof(argv[1]) : 0.1;
  const double en_factor = (argc > 2) ? atof(argv[2]) : 1.;

  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain);
  Mesher mesher(starset);
  mesher.refine_all();
  mesher.output("pointset_curvature_surface.off");
  
  delete pdomain;
  return 0;
}
