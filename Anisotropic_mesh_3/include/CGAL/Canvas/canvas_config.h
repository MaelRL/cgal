#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_CONFIG_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_CONFIG_H

#define ANISO_GEO_FILTER_SEEDS_OUTSIDE_CANVAS_BBOX
#define verbosity 11

  // if one wants to still initialize vertices even if the seed is found to be in
  // an exterior cell (might want to use that when there are seeds on the border
  // of the domain since numerical issues might not find it in the domain)
#define ANISO_GEO_FORCE_SEED_INITIALIZATION

#define CGAL_MESH_3_VERBOSE

#include <limits>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

const double FT_inf = std::numeric_limits<double>::infinity();

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_CONFIG_H
