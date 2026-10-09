# Default to just the release library until we support TPL debug builds on LC
set(QHULL_BUILD_TYPE Release CACHE STRING "qhull build type")
if (QHULL_BUILD_TYPE STREQUAL "Debug")
  set(${lib_name}_libs libqhullstatic_rd.a)
else()
  set(${lib_name}_libs libqhullstatic_r.a)
endif()
