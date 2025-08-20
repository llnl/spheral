spheral_add_test(
  NAME offload_tests
  SOURCES offload_tests.cc
  DEBUG_LINKER
)

spheral_add_test(
  NAME quadraticinterpolator_tests
  SOURCES quadraticinterpolator_tests.cc
  DEPENDS_ON Spheral_Utilities
  #DEBUG_LINKER
)

spheral_add_test(
  NAME chinterpolator_tests
  SOURCES chinterpolator_tests.cc
  DEPENDS_ON Spheral_Utilities
  #DEBUG_LINKER
)

spheral_add_test(
  NAME flexpointer_tests
  SOURCES flexpointer_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Geometry
  #DEBUG_LINKER
)
