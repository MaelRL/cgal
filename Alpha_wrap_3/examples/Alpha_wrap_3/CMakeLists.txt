# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Alpha_wrap_3_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 REQUIRED) # (3.1.0 or greater)
include(CGAL_Eigen3_support)

# create a target per cppfile
create_single_source_cgal_program("triangle_mesh_wrap.cpp")
target_link_libraries(triangle_mesh_wrap PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("triangle_soup_wrap.cpp")
target_link_libraries(triangle_soup_wrap PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("point_set_wrap.cpp")
target_link_libraries(point_set_wrap PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("wrap_from_cavity.cpp")
target_link_libraries(wrap_from_cavity PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("mixed_inputs_wrap.cpp")
target_link_libraries(mixed_inputs_wrap PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("volumetric_wrap.cpp")
target_link_libraries(volumetric_wrap PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("sparse_wrapping_with_region_growing.cpp")
target_link_libraries(sparse_wrapping_with_region_growing PUBLIC CGAL::Eigen3_support)
