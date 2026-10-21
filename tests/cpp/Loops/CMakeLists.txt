spheral_add_test(
  NAME forall_tests
  SOURCES forall_tests.cc
  DEPENDS_ON Spheral_CXX
  #DEBUG_LINKER
)

spheral_add_test(
  NAME sph_tests
  SOURCES sph_tests.cc
  DEPENDS_ON Spheral_CXX
  #DEBUG_LINKER
)
