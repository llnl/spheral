spheral_add_test(
  NAME forall_test
  SOURCES forall_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Material Spheral_NodeList Spheral_Neighbor Spheral_Utilities Spheral_Geometry Spheral_Distributed Spheral_DataOutput Spheral_FileIO Spheral_Kernel
  DEBUG_LINKER
)
