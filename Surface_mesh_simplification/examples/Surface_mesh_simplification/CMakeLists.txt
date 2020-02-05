# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Surface_mesh_simplification_Examples )

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

find_package( METIS )
if( METIS_FOUND )
  include_directories (${METIS_INCLUDE_DIRS})
else()
  message(STATUS "Examples that use the metis library will not be compiled.")
endif()

find_package( TBB )
if( NOT TBB_FOUND )
  message( STATUS "NOTICE: Intel TBB was not found. Examples will use sequential code." )
endif()

find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "edge_collapse_constrain_sharp_edges.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_enriched_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_linear_cell_complex.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_all_short_edges.cpp" )
create_single_source_cgal_program( "edge_collapse_bounded_normal_change.cpp" )
create_single_source_cgal_program( "edge_collapse_visitor_surface_mesh.cpp" )

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
if( Eigen3_FOUND )
  include( ${EIGEN3_USE_FILE} )
  create_single_source_cgal_program( "edge_collapse_garland_heckbert.cpp" )
else ()
  message(STATUS "Garland-Heckbert polices require the Eigen library, which has not been found; related examples will not be compiled.")
endif()

if(METIS_FOUND AND TBB_FOUND AND Eigen3_FOUND)
  create_single_source_cgal_program( "parallel_edge_collapse.cpp" )
  target_link_libraries( parallel_edge_collapse PRIVATE ${METIS_LIBRARIES} )
  CGAL_target_use_TBB( parallel_edge_collapse )
endif()

if(OpenMesh_FOUND)
  create_single_source_cgal_program( "edge_collapse_OpenMesh.cpp" )
  target_link_libraries( edge_collapse_OpenMesh PRIVATE ${OPENMESH_LIBRARIES} )
endif()
