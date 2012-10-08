
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Domain/Constrain_surface_3_ellipse.h>
#include <Domain/Constrain_surface_3_cylinder.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Metric_field/Arctan_metric_field.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>


using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef Anisotropic_surface_mesher_3<K> Mesher;

int main()
{
  CGAL::default_random = CGAL::Random(0);

  
//  Constrain_surface_3_ellipse<K>* pdomain 
//    = new Constrain_surface_3_ellipse<K>(1.0, 0.4, 0.6);

  Constrain_surface_3_cylinder<K>* pdomain
    = new Constrain_surface_3_cylinder<K>(1., 10.);

  

  Criteria_base<K> criteria(3.0, //radius_edge_ratio_
                            0.2, //sliverity_
                            0.1, //circumradius_
                            1.8, //distortion_
                            2.5, //beta_
                            0.3);//delta_
  Implicit_curvature_metric_field<K> metric_field(*pdomain, 1.0);
    

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain);
  Mesher mesher(starset);
  mesher.refine_all();
//  mesher.output("ellipse_curvature_surface.off");
  mesher.output("cylinder_curvature_surface.off");

  delete pdomain;
  return 0;
}
