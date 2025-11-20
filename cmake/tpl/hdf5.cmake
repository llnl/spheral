set(${lib_name}_libs libhdf5.so libhdf5_hl.so)

if(ENABLE_STATIC_TPL)
  string(REPLACE ".so" ".a;" ${lib_name}_libs ${${lib_name}_libs})
endif()
