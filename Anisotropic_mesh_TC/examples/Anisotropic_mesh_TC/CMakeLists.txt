# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Anisotropic_mesh_TC_examples)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories(../../include)
 
# Add macros
add_definitions(
	-DANISO_VERBOSE
	)

# Find CGAL 
find_package(CGAL COMPONENTS)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

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

  create_single_source_cgal_program("aniso_brute.cpp")

else()  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
