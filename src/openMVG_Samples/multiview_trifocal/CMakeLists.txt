add_executable(openMVG_sample_multiview_trifocal trifocal.cpp)
target_link_libraries(openMVG_sample_multiview_trifocal
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})
target_compile_definitions(openMVG_sample_multiview_trifocal
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_trifocal PROPERTY FOLDER OpenMVG/Samples)