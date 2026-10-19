spheral_add_test(
  NAME artvisc_tests
  SOURCES artvisc_tests.cc
  DEPENDS_ON Spheral_CXX
  #DEBUG_LINKER
)

spheral_add_test(
  NAME tensorav_tests
  SOURCES tensorav_tests.cc
  DEPENDS_ON Spheral_CXX
  #DEBUG_LINKER
)
