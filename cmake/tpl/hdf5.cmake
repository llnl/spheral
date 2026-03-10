set(${lib_name}_libs libhdf5.dylib libhdf5_hl.dylib)

if(ENABLE_STATIC_TPL)
  string(REPLACE ".dylib" ".a;" ${lib_name}_libs ${${lib_name}_libs})
endif()
