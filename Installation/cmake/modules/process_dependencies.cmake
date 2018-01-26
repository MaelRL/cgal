foreach(n RANGE 5 ${CMAKE_ARGC})
  list(APPEND INPUT_FILES ${CMAKE_ARGV${n}})
endforeach()

if(NOT CGAL_PACKAGES_PREFIX)
  message(FATAL_ERROR
    "The variable `CGAL_PACKAGES_PREFIX` should be defined to the prefix of CGAL packages!")
endif()

#message("regexp: \\. ${CGAL_PACKAGES_PREFIX}/[^/]*/include/CGAL/.*h")
#TO BE CHECKED
get_filename_component(INSTALLATION "${CMAKE_BINARY_DIR}" DIRECTORY)
get_filename_component(BUILD_DIR "${INSTALLATION}" NAME)
foreach(INPUT_FILE ${INPUT_FILES})
  file(STRINGS ${INPUT_FILE} input)
  #message("input is : ${input}")
  foreach(line ${input})
    string(REGEX MATCHALL "^\\.* ${CGAL_PACKAGES_PREFIX}/[A-Za-z0-9_.-]*/include/CGAL/[A-Za-z0-9_/.-]*\\.h" header ${line})
    string(REGEX REPLACE "\\.* ${CGAL_PACKAGES_PREFIX}/" "" header "${header}")
    string(REGEX REPLACE "/.*" "" pkg "${header}")
    string(REPLACE "${BUILD_DIR}" " " pkg "${pkg}")
    if(header)
      list(APPEND headers ${header})
    endif()
    if(pkg)
      list(APPEND pkgs ${pkg})
    endif()
  endforeach()
endforeach()
if(headers)
  list(REMOVE_DUPLICATES headers)
  list(SORT headers)
endif()
if(pkgs)
  list(REMOVE_DUPLICATES pkgs)
  list(SORT pkgs)
endif()
if(OUTPUT_HEADERS_LIST)
  file(WRITE ${OUTPUT_HEADERS_LIST} "")
  foreach(header ${headers})
    file(APPEND ${OUTPUT_HEADERS_LIST} "${header}\n")
  endforeach()
endif()
if(OUTPUT_PACKAGES_LIST)
  file(WRITE ${OUTPUT_PACKAGES_LIST} "")
  foreach(pkg ${pkgs})
    file(APPEND ${OUTPUT_PACKAGES_LIST} "${pkg}\n")
  endforeach()
endif()
