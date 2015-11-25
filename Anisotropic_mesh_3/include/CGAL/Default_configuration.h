//#define CGAL_PROFILE

#ifndef CGAL_ANISOTROPIC_MESH_3_DEFAULT_CONFIGURATION_H
#define CGAL_ANISOTROPIC_MESH_3_DEFAULT_CONFIGURATION_H

// whether or not to check edge encroachment.
// if the surface has sharp edges and defined in the class,
// please uncomment the following macro definition.
//#define CHECK_EDGE_ENCROACHMENT

// whether or not to check surface encroachment.
// for volume mesher, if we do not check the encroachment of
// the surface domain, comment the following macro definition.
//#define CHECK_FACE_ENCROACHMENT

// for volume mesher, if checking encroachment, only
// check the encroachment of the neigbhoring facets.
// this can accelerate the process.
//#define ONLY_CONSIDER_NEIGHBORING_ENCROACHMENT

// show debug info.
#define CGAL_ANISOTROPIC_MESH_3_DEBUG_INFO

// do not pick points in the picking region
//#define NO_PICKING_REGION

// before creating a star, always check whether or
// not the new star is the same with existing stars.
//#define CHECK_EXISTING_VERTICES

// sometimes statistics are not needed/wanted
#define MAKE_STATISTICS

//#define VERBOSE_REFINEMENT 0

//#define USE_AABB_TREE_OF_BBOXES  1

#endif // CGAL_ANISOTROPIC_MESH_3_DEFAULT_CONFIGURATION_H
