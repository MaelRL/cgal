#include "StarSet_type.h"
#include "Scene_starset3_item.h"
#include "Anisotropic_meshing_thread.h"
#include "Anisotropic_mesh_function.h"
#include "Implicit_surface_type.h"

#include <CGAL/Bbox_3.h>

//#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Metric_field/Torus_metric_field.h>


Anisotropic_meshing_thread* cgal_code_anisotropic_mesh_3(const Implicit_surface* p_surface,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const int max_times_to_try_in_picking_region,
                                 const int dim,
                                 const int nb_initial_points)
{
  if( NULL == p_surface ) { return NULL; }
    
  const Implicit_surface* p_domain = p_surface->clone();
  
//  typedef CGAL::Anisotropic_mesh_3::Euclidean_metric_field<Kernel> Metric_field;
//  std::cout << "(Euclidean)." << std::endl;

//  typedef CGAL::Anisotropic_mesh_3::Implicit_curvature_metric_field<Kernel> Metric_field;
//  std::cout << "(Curvature metric field)." << std::endl;
  
  typedef Torus_metric_field<Kernel> Metric_field;
  

  Metric_field* metric_field = new Metric_field(0.7, 0.3, 1.);//*p_domain);;
  
  typedef Anisotropic_mesh_function<Implicit_surface, Metric_field> AMesh_function;
 
     
  Anisotropic_mesh_parameters param;
  param.radius_edge_ratio = radius_edge_ratio;
  param.sliverity = sliverity;
  param.circumradius = circumradius;
  param.distortion = distortion;
  param.beta = beta;
  param.delta = delta;
  param.max_times_to_try_in_picking_region = max_times_to_try_in_picking_region;
  param.dim = dim;

  // @TODO, WARNING: memory leak to be corrected later: criteria and
  // metric_field must be destroyed by somebody. The issue is that they
  // cannot be destroyed before the life end of the meshing thread.
  Criteria* criteria = new Criteria(param.radius_edge_ratio, param.sliverity, param.circumradius, 
                                    param.distortion, param.beta, param.delta, 
                                    param.max_times_to_try_in_picking_region);

  Scene_starset3_item* p_new_item 
    = new Scene_starset3_item(*criteria, *metric_field, p_domain, nb_initial_points);

  AMesh_function* p_mesh_function 
    = new AMesh_function(p_new_item->star_set(), param, criteria, metric_field);
  // The mesh function takes the ownership of 'criteria' and
  // 'metric_field', to release them at its destruction.

  return new Anisotropic_meshing_thread(p_mesh_function, p_new_item); 
}
