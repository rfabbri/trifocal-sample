set(trifocal_sources
  trifocal.cpp trifocal.h  # TODO(future): make into .hpp header-only lib
  trifocal-util.cpp trifocal-util.h
)
set(trifocal_app_sources
  trifocal-app.cpp trifocal-app.h
) 
set(trifocal_sample_sources
  trifocal-sample.cpp 
) 

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# -- MINUS ---------------------------------------------------------------------
# optional, header-only library can be used without this library
# library is useful when the function will be used in many places
# to avoid having to include Minus many times and slowdown the compiler
#
# For OpenMVG it is worth doing this so that MINUS is always compiled with
# highly optimized compiler flags, even when OpenMVG is debug
set(minus_sources
  minus/minus.hxx minus/minus.h 
#
  minus/chicago14a.h minus/chicago-default.h minus/chicago14a-default-data.h minus/chicago14a-io.h
  minus/chicago14a.hxx minus/chicago14a-default-data.hxx 
#
  minus/problem-defs.h
  
  minus/Templates/minus-chicago-alltypes-allformulations+double.cxx
  minus/Templates/minus_io_common+double.cxx
  minus/Templates/minus_util+double.cxx
)
# this library can be used if you want to precompile instantiations
# but header-only will also work, see README.md
add_library( minus ${minus_sources} )
# always build minus with specially tweaked optimizations.
# these were copied directly from minus original flags
set_target_properties(minus PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(minus PROPERTIES COMPILE_FLAGS "-O3 -march=native -ffast-math")
set_target_properties(minus PROPERTIES COMPILE_DEFINITIONS "NDEBUG")
# ------------------------------------------------------------------------------

add_library(openMVG_trifocal ${trifocal_sources})
add_library(openMVG_trifocal_app ${trifocal_app_sources})

target_link_libraries(openMVG_trifocal
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY} minus)

target_link_libraries(openMVG_trifocal_app openMVG_trifocal)

UNIT_TEST(openMVG trifocal_stub "openMVG_trifocal")
UNIT_TEST(openMVG trifocal "openMVG_trifocal_app")

# set(CMAKE_CXX_FLAGS "-g") # XXX temporary for build cycle
# set(CMAKE_C_FLAGS "-g") # XXX temporary for build cycle
# These work better: they add -g but they won't override previous flags
set_target_properties(openMVG_test_trifocal PROPERTIES COMPILE_FLAGS "-g")
#set_target_properties(openMVG_test_trifocal PROPERTIES COMPILE_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
set_target_properties(openMVG_test_trifocal_stub PROPERTIES COMPILE_FLAGS "-g")
set_target_properties(openMVG_trifocal PROPERTIES COMPILE_FLAGS "-g")
set_target_properties(openMVG_trifocal_app PROPERTIES COMPILE_FLAGS "-g")

add_executable(openMVG_sample_multiview_trifocal ${trifocal_sample_sources})
target_link_libraries(openMVG_sample_multiview_trifocal openMVG_trifocal_app)
target_compile_definitions(openMVG_sample_multiview_trifocal
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_trifocal PROPERTY FOLDER OpenMVG/Samples)
