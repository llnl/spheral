set(${lib_name}_libs libsiloh5.a)

if(APPLE)
  set(${lib_name}_libs libsiloh5.dylib)
endif()
Spheral_Handle_Ext(${lib_name} APPLE)
