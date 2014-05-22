#include <CGAL/Default_configuration.h>

#include "StarSet_type.h"
#include "Scene_starset3_item.h"
#include "Polyhedral_surface_type.h" //contains Kernel
#include "Anisotropic_meshing_thread.h"
#include "Anisotropic_mesh_function.h"
#include "anisotropic_meshing_options.h"

#include <CGAL/Metric_field.h>
#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/Random.h>

Criteria* build_criteria_and_metric(const Polyhedral_surface* const p_domain,
                                    CGAL::Anisotropic_mesh_3::Metric_field<Kernel>* & mf,
                                    const int dimension,
                                    const double approximation,
                                    const double facet_circumradius,
                                    const double facet_radius_edge_ratio,
                                    const double cell_circumradius,
                                    const double cell_radius_edge_ratio,
                                    const double sliverity,
                                    const double distortion,
                                    const double beta,
                                    const double delta,
                                    const int nb_initial_points,
                                    const std::size_t max_times_to_try_in_picking_region,
                                    const Metric_options& metric,
                                    const double epsilon,
                                    const double en_factor)
{
  if(metric == EUCLIDEAN)
  {
    std::cout << "(Metric : Euclidean)." << std::endl;
    mf = new CGAL::Anisotropic_mesh_3::Euclidean_metric_field<Kernel>(1., 1., 1., epsilon, en_factor);
  }
  else if(metric == POLYHEDRON_CURVATURE)
  {
    std::cout << "(Metric : Polyhedral)" << std::endl;
    typedef CGAL::Anisotropic_mesh_3::Polyhedral_curvature_metric_field<Kernel> PMF;
    mf = new PMF(*p_domain, epsilon, en_factor);
  }
  else if(metric == HYPERBOLIC_SHOCK)
  {
    std::cout << "(Hyperbolic shock metric field)." << std::endl;
    mf = new Hyperbolic_shock_metric_field<Kernel>(0.6, epsilon, en_factor);
  }
  std::cout << "epsilon/en_factor : " << epsilon << " " << en_factor << std::endl;

  // @TODO, WARNING: memory leak to be corrected later: criteria and
  // metric_field must be destroyed by somebody. The issue is that they
  // cannot be destroyed before the life end of the meshing thread.
  return new Criteria(approximation, facet_circumradius, facet_radius_edge_ratio,
                      cell_circumradius, cell_radius_edge_ratio, sliverity,
                      distortion, beta, delta, nb_initial_points,
                      max_times_to_try_in_picking_region, dimension);
}

Anisotropic_meshing_thread* cgal_code_anisotropic_mesh_3(const Polyhedron* p_poly,
                                                         const int dimension,
                                                         const double approximation,
                                                         const double facet_circumradius,
                                                         const double facet_radius_edge_ratio,
                                                         const double cell_circumradius,
                                                         const double cell_radius_edge_ratio,
                                                         const double sliverity,
                                                         const double distortion,
                                                         const double beta,
                                                         const double delta,
                                                         const int nb_initial_points,
                                                         const std::size_t max_times_to_try_in_picking_region,
                                                         const Metric_options& metric,
                                                         const double epsilon,
                                                         const double en_factor)
{
  CGAL::default_random = CGAL::Random(0);

  if( NULL == p_poly ) { return NULL; }
  const Polyhedral_surface* const p_domain = new Polyhedral_surface(*p_poly);

  Criteria* criteria = NULL;
  typedef CGAL::Anisotropic_mesh_3::Metric_field<Kernel> Metric_field;
  Metric_field* mf = NULL;

  criteria = build_criteria_and_metric(p_domain, mf, dimension, approximation, facet_circumradius,
                                       facet_radius_edge_ratio, cell_circumradius,
                                       cell_radius_edge_ratio, sliverity,
                                       distortion, beta, delta, nb_initial_points,
                                       max_times_to_try_in_picking_region,
                                       metric, epsilon, en_factor);

  Scene_starset3_item* p_new_item = new Scene_starset3_item(p_domain, mf, criteria);

  typedef Anisotropic_mesh_function<Polyhedral_surface, Metric_field> AMesh_function;
  AMesh_function* p_mesh_function = new AMesh_function(p_new_item->star_set());
  // The mesh function takes the ownership of 'criteria' and
  // 'metric_field', to release them at its destruction.

  return new Anisotropic_meshing_thread(p_mesh_function, p_new_item);
}
