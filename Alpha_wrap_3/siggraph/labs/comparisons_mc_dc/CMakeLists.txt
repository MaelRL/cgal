# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( Isosurfacing )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package

find_package(OpenMP)


# Creating entries for target: Isosurfacing
# ############################

file(GLOB SOURCES 
  Isosurfacing_3/*.h 
  Isosurfacing_3/*/*.h 
  Isosurfacing_3/*/*/*.h 
)

file(GLOB EXECUTABLES Isosurfacing_3/*.cpp)

foreach(exec ${EXECUTABLES})
  get_filename_component(exec_name ${exec} NAME_WE)
  message(STATUS "Add executable ${exec_name}")
  
  add_executable( ${exec_name} ${SOURCES} ${exec} )

  # Link the executable to CGAL and third-party libraries
  target_link_libraries( ${exec_name} PRIVATE CGAL::CGAL )

  if(OpenMP_CXX_FOUND)
      target_link_libraries( ${exec_name} PRIVATE OpenMP::OpenMP_CXX )
  endif()
  
endforeach()

