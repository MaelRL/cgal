if (CGAL_target_use_TBB_included)
  return()
endif()
set(CGAL_target_use_TBB_included TRUE)

function(CGAL_target_use_TBB target)
  target_include_directories ( ${target} SYSTEM PRIVATE ${TBB_INCLUDE_DIRS} )
  target_link_libraries( ${target} PRIVATE ${TBB_LIBRARIES} )
  target_compile_options( ${target} PRIVATE -DNOMINMAX -DCGAL_LINKED_WITH_TBB )
endfunction()
