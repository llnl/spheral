
add_subdirectory(DeviceTestLib)

if(ENABLE_TESTS)
  if(ENABLE_DEV_BUILD)
    set(SPHERAL_DEPEND Spheral_DeviceTestLib)
  else()
    set(SPHERAL_DEPEND Spheral_CXX)
  endif()
  blt_add_executable( NAME spheral_offload_test
                      SOURCES spheral_offload_test.cc
                      DEPENDS_ON ${SPHERAL_DEPEND} ${spheral_blt_cxx_depends})
endif()
