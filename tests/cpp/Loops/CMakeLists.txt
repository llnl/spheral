spheral_add_test(
  NAME forall_test
  SOURCES forall_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Material Spheral_NodeList Spheral_Neighbor Spheral_Utilities Spheral_Geometry Spheral_Distributed Spheral_DataOutput Spheral_FileIO Spheral_Kernel
  DEBUG_LINKER
)

spheral_add_test(
  NAME sph_test
  SOURCES sph_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Material Spheral_NodeList Spheral_Neighbor Spheral_Utilities Spheral_Geometry Spheral_Distributed Spheral_DataOutput Spheral_FileIO Spheral_Kernel
  #DEBUG_LINKER
)

spheral_add_test(
  NAME novirt_test
  SOURCES novirt_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Material Spheral_NodeList Spheral_Neighbor Spheral_Utilities Spheral_Geometry Spheral_Distributed Spheral_DataOutput Spheral_FileIO Spheral_Kernel
  #DEBUG_LINKER
)
