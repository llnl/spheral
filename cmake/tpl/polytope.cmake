set(${lib_name}_libs libpolytope.a)

if(APPLE)
  set(${lib_name}_libs libpolytope.dylib)
endif()
Spheral_Handle_Ext(${lib_name} APPLE)
