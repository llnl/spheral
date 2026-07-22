spheral_add_test(
  NAME forall_tests
  SOURCES forall_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Material Spheral_NodeList Spheral_Neighbor Spheral_Utilities Spheral_Geometry Spheral_Distributed Spheral_DataOutput Spheral_FileIO Spheral_Kernel Spheral_Threading
  #DEBUG_LINKER
)

spheral_add_test(
  NAME sph_tests
  SOURCES sph_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Material Spheral_NodeList Spheral_Neighbor Spheral_Utilities Spheral_Geometry Spheral_Distributed Spheral_DataOutput Spheral_FileIO Spheral_Kernel Spheral_Threading
  #DEBUG_LINKER
)
