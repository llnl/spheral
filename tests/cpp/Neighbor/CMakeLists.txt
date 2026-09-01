spheral_add_test(
  NAME nodepairlistview_test
  SOURCES nodepairlistview_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Neighbor chai
  #DEBUG_LINKER #In case it has unresolved symbols
)

spheral_add_test(
  NAME pairwisefieldview_test
  SOURCES pairwisefieldview_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Neighbor chai
  #DEBUG_LINKER #In case it has unresolved symbols
)

spheral_add_test(
  NAME treeneighborview_test
  SOURCES treeneighborview_tests.cc
  DEPENDS_ON Spheral_Utilities Spheral_Geometry Spheral_Field Spheral_Distributed
             Spheral_DataOutput Spheral_NodeList Spheral_Neighbor Spheral_Threading chai
)
