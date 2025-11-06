
add_subdirectory(DeviceTestLib)

if(ENABLE_TESTS)
  if(ENABLE_DEV_BUILD)
    set(SPHERAL_DEPEND Spheral_DeviceTestLib)
  else()
    set(SPHERAL_DEPEND Spheral_CXX)
  endif()
  blt_add_executable(NAME spheral_offload_test
    SOURCES spheral_offload_test.cc
    DEPENDS_ON ${SPHERAL_DEPEND} ${SPHERAL_CXX_DEPENDS})
  get_property(SPHERAL_CXX_FLAGS GLOBAL PROPERTY SPHERAL_CXX_FLAGS)
  get_property(SPHERAL_LINK_FLAGS GLOBAL PROPERTY SPHERAL_LINK_FLAGS)
  target_compile_options(spheral_offload_test PUBLIC ${SPHERAL_CXX_FLAGS})
  target_link_options(spheral_offload_test PUBLIC ${SPHERAL_LINK_FLAGS})
endif()
