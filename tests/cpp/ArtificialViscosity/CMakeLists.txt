spheral_add_test(
  NAME artvisc_tests
  SOURCES artvisc_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Geometry Spheral_NodeList Spheral_Utilities Spheral_Kernel Spheral_ArtificialViscosity Spheral_DataOutput Spheral_FileIO Spheral_Hydro Spheral_Distributed Spheral_Physics
  #DEBUG_LINKER
)
