spheral_add_test(
  NAME kernel_tests
  SOURCES kernel_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Kernel
  #DEBUG_LINKER
)

spheral_add_test(
  NAME tablekernel_tests
  SOURCES tablekernel_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Kernel
  #DEBUG_LINKER
)
