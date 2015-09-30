#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_ENUM_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_ENUM_H

namespace CGAL
{
namespace Anisotropic_mesh_3
{

enum FMM_state
{
  CHANGED = 0, // used only in the Konukoglu algorithm. 'CHANGED' implies 'KNOWN'
  KNOWN,
  TRIAL,
  FAR
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL,
  REBUILD_CHANGED, // used only in the Konukoglu algorithm
  REBUILD_BOTH // used only in the Konukoglu algorithm
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_ENUM_H
