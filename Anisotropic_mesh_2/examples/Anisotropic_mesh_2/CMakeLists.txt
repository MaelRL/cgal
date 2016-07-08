# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Anisotropic_mesh_2_examples)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories(../../include)
 
# Add macros
add_definitions(
	-DANISO_VERBOSE
#	-DANISO_DEBUG
#	-DANISO_DEBUG_REFINEMENT
#	-DANISO_DEBUG_METRIC
	-DANISO_OUTPUT_MESH_FOR_POLES
	-DCGAL_DEBUG_OUTPUT_VECTOR_FIELD
#	-DANISO_DUMP_INTERMEDIATE_MESHES
	-DANISO_APPROXIMATE_SPHERE_BBOX 
#	-DANISO_BRUTE_FORCE_SIMULATE_INSERT_TO_STAR
	)

 
# Find CGAL 
find_package(CGAL COMPONENTS)


if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  # Find Boost
  find_package(Boost COMPONENTS chrono timer)
  include_directories(${Boost_INCLUDE_DIRS})
  link_libraries(${Boost_LIBRARIES})

  include( CGAL_CreateSingleSourceCGALProgram )

  # Find OpenGL
  find_package(OpenGL)
  link_libraries(${OPENGL_LIBRARIES})
  
  # use Eigen if it is installed
  find_package(Eigen3 ) #(requires 3.1.0-alpha1 or greater)
  if (NOT EIGEN3_FOUND)
  message ( STATUS "eigen not found" )
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()

  create_single_source_cgal_program("mesher_example.cpp")
  create_single_source_cgal_program("histogram_from_external_mesh.cpp")
  create_single_source_cgal_program("grid_generator.cpp")
  create_single_source_cgal_program("geodesic_grid_generator.cpp")
  create_single_source_cgal_program("geodesic_grid_generator_2.cpp")
  create_single_source_cgal_program("geodesic_generator_on_tr.cpp")
  create_single_source_cgal_program("geodesic_generator_on_starset.cpp")
  create_single_source_cgal_program("Nash.cpp")

else()  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
