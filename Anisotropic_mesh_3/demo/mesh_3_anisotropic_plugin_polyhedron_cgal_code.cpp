#include <CGAL/Default_configuration.h>

#include "StarSet_type.h"
#include "Scene_starset3_item.h"
#include "Polyhedron_type.h" //contains Kernel
#include "Anisotropic_meshing_thread.h"
#include "Anisotropic_mesh_function.h"
#include "anisotropic_meshing_options.h"

#include <CGAL/Metric_field.h>
#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/Random.h>

Criteria* build_param_and_metric(const CGAL::Anisotropic_mesh_3::Constrain_surface_3_polyhedral<Kernel, Polyhedron>* const p_domain,
                                 Anisotropic_mesh_parameters & param,
                                 CGAL::Anisotropic_mesh_3::Metric_field<Kernel>* & mf,
                                 const double epsilon,
                                 const double approximation,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const std::size_t max_times_to_try_in_picking_region,
                                 const int dim,
                                 const Metric_options& metric,
                                 const double en_factor)
{
  param.approximation = approximation;
  param.radius_edge_ratio = radius_edge_ratio;
  param.sliverity = sliverity;
  param.circumradius = circumradius;
  param.distortion = distortion;
  param.beta = beta;
  param.delta = delta;
  param.max_times_to_try_in_picking_region = max_times_to_try_in_picking_region;
  param.dim = dim;

  //typedef CGAL::Anisotropic_mesh_3::Metric_field<Kernel> Metric_field;

  if(metric == EUCLIDEAN)
  {
    std::cout << "(Metric : Euclidean)." << std::endl;
    mf = new CGAL::Anisotropic_mesh_3::Euclidean_metric_field<Kernel>(1., 1., 1., epsilon, en_factor);
  }
  else if(metric == POLYHEDRON_CURVATURE)
  {
    std::cout << "(Metric : Polyhedral)" << std::endl;
    typedef CGAL::Anisotropic_mesh_3::Polyhedral_curvature_metric_field<Kernel> PMF;
    mf = new PMF(*p_domain, en_factor);
  }
  else if(metric == HYPERBOLIC_SHOCK)
  {
    std::cout << "(Hyperbolic shock metric field)." << std::endl;
    mf = new Hyperbolic_shock_metric_field<Kernel>(0.6, epsilon, en_factor);
  }

  // @TODO, WARNING: memory leak to be corrected later: criteria and
  // metric_field must be destroyed by somebody. The issue is that they
  // cannot be destroyed before the life end of the meshing thread.
  return new Criteria(param.radius_edge_ratio, param.sliverity, param.circumradius,
                                    param.distortion, param.beta, param.delta,
                                    param.max_times_to_try_in_picking_region, param.approximation);
}

Anisotropic_meshing_thread* cgal_code_anisotropic_mesh_3(const Polyhedron* p_poly,
                                 const double epsilon,
                                 const double approximation,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const std::size_t max_times_to_try_in_picking_region,
                                 const int dim,
                                 const int nb_initial_points,
                                 const Metric_options& metric,
                                 const bool pick_valid_causes_stop,
                                 const int pick_valid_max_failures,
                                 const double en_factor)
{
  CGAL::default_random = CGAL::Random(0);

  if( NULL == p_poly ) { return NULL; }
  typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_polyhedral<Kernel, Polyhedron>  Constrain_surface_polyhedral;
  const Constrain_surface_polyhedral* const p_domain = new Constrain_surface_polyhedral(*p_poly,epsilon);
  
    //ini
  Anisotropic_mesh_parameters param;
  Criteria* criteria = NULL;
  typedef CGAL::Anisotropic_mesh_3::Metric_field<Kernel> Metric_field;
  Metric_field* mf = NULL;

    //build
  criteria = build_param_and_metric(p_domain, param, mf, epsilon, approximation, radius_edge_ratio,
                                    sliverity, circumradius, distortion, beta, delta,
                                    max_times_to_try_in_picking_region, dim, metric, en_factor);

  Scene_starset3_item* p_new_item 
    = new Scene_starset3_item(criteria, mf, p_domain, nb_initial_points);

  typedef Anisotropic_mesh_function<Constrain_surface_polyhedral, Metric_field> AMesh_function;
  AMesh_function* p_mesh_function = new AMesh_function(p_new_item->star_set(), param, criteria,
                                                       mf, pick_valid_causes_stop,
                                                       pick_valid_max_failures);
  // The mesh function takes the ownership of 'criteria' and
  // 'metric_field', to release them at its destruction.

  return new Anisotropic_meshing_thread(p_mesh_function, p_new_item);
}
