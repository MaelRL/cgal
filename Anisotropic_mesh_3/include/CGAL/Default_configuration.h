// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Kan-Le Shi

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
#define CHECK_FACE_ENCROACHMENT

// if the following switcher is on, the insertion operation
// is gurenteed to be O(1) so that the entire algorithm can be
// O(n), where n is the number of vertices.
#define CREATE_STAR_ENHANCEMENT

// if the following switcher is on, the insertion operation
// is guaranteed to be O(1) so that the entire algorithm can be
// O(n), where n is the number of vertices.
// we cannot guarantee the correctness of the result if the
// following switcher is on.
// #define INSERT_POINT_ENHANCEMENT

// always on
#define ONLY_CONSIDER_SURFACE

// for volume mesher, if checking encroachment, only
// check the encroachment of the neigbhoring facets.
// this can accelarate the speed.
#define ONLY_CONSIDER_NEIGHBORING_ENCROACHMENT

// show debug info.
#define CGAL_ANISOTROPIC_MESH_3_DEBUG_INFO

// do not pick points in the picking region
//#define NO_PICKING_REGION

// before creating a star, always check whether or 
// not the new star is the same with existing stars.
//#define CHECK_EXISTING_VERTICES

// sometimes statistics are not needed/wanted
#define MAKE_STATISTICS

#define VERBOSE_REFINEMENT 0

#define USE_AABB_TREE_OF_BBOXES  1

#endif // CGAL_ANISOTROPIC_MESH_3_DEFAULT_CONFIGURATION_H