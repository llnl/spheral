spheral_add_test(
  NAME nodepairlistview_test
  SOURCES nodepairlistview_tests.cc
  DEPENDS_ON Spheral_CXX chai
  #DEBUG_LINKER #In case it has unresolved symbols
)

spheral_add_test(
  NAME pairwisefieldview_test
  SOURCES pairwisefieldview_tests.cc
  DEPENDS_ON Spheral_CXX chai
  #DEBUG_LINKER #In case it has unresolved symbols
)
