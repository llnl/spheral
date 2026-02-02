spheral_add_test(
  NAME nodepairlistview_test
  SOURCES nodepairlistview_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Neighbor chai
  #DEBUG_LINKER #In case it has unresolved symbols
)
