# This is the CMake script for compiling the CGAL Mesh_3 demo implicit functions.

project( Mesh_3_constrained_surfaces_functions )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include directory of demo includes
include_directories( BEFORE ${Mesh_3_constrained_surfaces_functions_BINARY_DIR} ../include )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

if(WIN32)
 # Creates a new CMake option, turned ON by default
 option(ACTIVATE_MSVC_PRECOMPILED_HEADERS
   "Activate precompiled headers in MSVC"
   ON)
 
 # Macro to add precompiled headers for MSVC
 MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)
   IF(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS)
     GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
     SET(Sources ${${SourcesVar}})
 
     SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                 PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\"")
     SET_SOURCE_FILES_PROPERTIES(${Sources}
                                 PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledHeaders}\" /FI\"${PrecompiledHeader}\"")
     # Add precompiled header to SourcesVar
     LIST(APPEND ${SourcesVar} ${PrecompiledSource})
   ENDIF(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS)
 ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)
endif()


if(CGAL_Qt4_FOUND AND QT4_FOUND)
  include(${QT_USE_FILE})

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  ###########
  # PLUGINS #
  ###########
  include(AddFileDependencies)
  remove_definitions(-DQT_STATICPLUGIN)

  # polyhedron_demo_plugin is a macro copied from demo/Polyhedron/CMakeLists.txt
  macro(polyhedron_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      set(other_sources ${ARGN})
      set(option "")
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  set(CHAIR_SRC_FILES "Chair_function.cpp")
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" CHAIR_SRC_FILES)
  LIST(REMOVE_AT CHAIR_SRC_FILES 0) 
  polyhedron_demo_plugin(chair_constrained_surface_plugin Chair_function ${CHAIR_SRC_FILES})

  set(SPHERE_SRC_FILES "Sphere_function.cpp")
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" SPHERE_SRC_FILES)
  LIST(REMOVE_AT SPHERE_SRC_FILES 0) 
  polyhedron_demo_plugin(sphere_constrained_surface_plugin Sphere_function ${SPHERE_SRC_FILES})

  set(TANGLECUBE_SRC_FILES "Tanglecube_function.cpp")
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" TANGLECUBE_SRC_FILES)
  LIST(REMOVE_AT TANGLECUBE_SRC_FILES 0) 
  polyhedron_demo_plugin(tanglecube_constrained_surface_plugin Tanglecube_function ${TANGLECUBE_SRC_FILES})

  set(TORUS_SRC_FILES "Torus_function.cpp")
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" TORUS_SRC_FILES)
  LIST(REMOVE_AT TORUS_SRC_FILES 0) 
  polyhedron_demo_plugin(torus_constrained_surface_plugin Torus_function ${TORUS_SRC_FILES})

  set(ELLIPSOID_SRC_FILES "Ellipsoid_function.cpp")
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" ELLIPSOID_SRC_FILES)
  LIST(REMOVE_AT ELLIPSOID_SRC_FILES 0) 
  polyhedron_demo_plugin(ellipsoid_constrained_surface_plugin Ellipsoid_function ${ELLIPSOID_SRC_FILES})
  
else (CGAL_Qt4_FOUND AND QT4_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt4 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(MESH_3_MISSING_DEPS "Qt4, ${MESH_3_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${MESH_3_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND)
