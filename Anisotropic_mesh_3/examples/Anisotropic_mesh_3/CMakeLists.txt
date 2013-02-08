# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Anisotropic_mesh_3_examples)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories(../../include)
 
# Add macros
add_definitions(
	-DANISO_VERBOSE
#	-DANISO_DEBUG
	-DANISO_DEBUG_REFINEMENT
	-DANISO_DEBUG_METRIC
    -DANISO_APPROXIMATE_SPHERE_BBOX 
	-DANISO_BRUTE_FORCE_SIMULATE_INSERT_TO_STAR)

 
# Find CGAL 
find_package(CGAL COMPONENTS)


if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # use Eigen if it is installed
  find_package(Eigen3 ) #(requires 3.1.0-alpha1 or greater)
  if (NOT EIGEN3_FOUND)
  message ( STATUS "eigen not found" )
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()

  # Find LAPACK (optional), for curvatures estimation
  find_package(LAPACK) 
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
	
#   create_single_source_cgal_program("implicit_surface_with_specified_metric_field.cpp")
    create_single_source_cgal_program("implicit_surface_with_curvature_mf.cpp")
#   create_single_source_cgal_program("polyhedral_surface_with_curvature_mf.cpp")
#   create_single_source_cgal_program("torus_example.cpp")
#   create_single_source_cgal_program("cylinder_example.cpp")
#	create_single_source_cgal_program("bench_parameters.cpp")
#	create_single_source_cgal_program("chair_example.cpp")
#   create_single_source_cgal_program("debug_example.cpp")	
#	create_single_source_cgal_program("comparisons_example.cpp")	
#	create_single_source_cgal_program("metric_field_experiment.cpp")	
		
  else(LAPACK_FOUND)

    message(STATUS "NOTICE: Some of the executables in this directory need the LAPACK library and will not be compiled.")

  endif(LAPACK_FOUND)
  
else()  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
