# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Anisotropic_mesh_3_examples)

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
  find_package(Boost COMPONENTS chrono)
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

  # use tbb
  find_package( TBB )
  if(TBB_FOUND)
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif(TBB_FOUND)
 
  # Find LAPACK (optional), for curvatures estimation
  find_package(LAPACK) 
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )

  create_single_source_cgal_program("implicit_surface_with_specified_metric_field.cpp")
  create_single_source_cgal_program("implicit_surface_with_curvature_mf.cpp")
  create_single_source_cgal_program("ellipsoid_example.cpp")
  create_single_source_cgal_program("double_ellipsoid_example.cpp")
  create_single_source_cgal_program("polyhedral_surface_with_curvature_mf.cpp")
  create_single_source_cgal_program("torus_example.cpp")
  create_single_source_cgal_program("cylinder_example.cpp")
  create_single_source_cgal_program("bench_parameters.cpp")
  create_single_source_cgal_program("chair_example.cpp")
  create_single_source_cgal_program("debug_example.cpp")
  create_single_source_cgal_program("comparisons_example.cpp")
  create_single_source_cgal_program("metric_field_experiment.cpp")
  create_single_source_cgal_program("medial_axis_with_powercrust.cpp")
#  create_single_source_cgal_program("point_set_poisson_with_curvature_mf.cpp")
  create_single_source_cgal_program("compare_ellipsoids_iso_aniso.cpp")
  create_single_source_cgal_program("cell_mesher_example.cpp")
  create_single_source_cgal_program("cell_mesher_poly_example.cpp")
  create_single_source_cgal_program("mesher_lvl_example.cpp")
  create_single_source_cgal_program("small_angle_between_planes_example.cpp")
  create_single_source_cgal_program("histogram_from_external_mesh.cpp")
  create_single_source_cgal_program("Voronoi_painters.cpp")
  create_single_source_cgal_program("geodesic_Konukoglu_grid.cpp")
  create_single_source_cgal_program("geodesic_Campen_grid.cpp")
  create_single_source_cgal_program("geodesic_Campen_c2t3.cpp")
  create_single_source_cgal_program("geodesic_Campen_c3t3.cpp")
  create_single_source_cgal_program("geodesic_Campen_surf_starset.cpp")
  create_single_source_cgal_program("geodesic_Campen_starset.cpp")
  create_single_source_cgal_program("test_intersections.cpp")
  
  else(LAPACK_FOUND)
    message(STATUS "NOTICE: Some of the executables in this directory need the LAPACK library and will not be compiled.")
  endif(LAPACK_FOUND)
  
else()  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
