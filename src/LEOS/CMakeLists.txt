include_directories(.)
set(LEOS_inst
  LEOS
)

set(LEOS_sources
  LEOS_bundle.cc
)

instantiate(LEOS_inst LEOS_sources)

set(LEOS_headers
  LEOS.hh
  LEOS_bundle.hh
)

spheral_add_obj_library(LEOS SPHERAL_OBJ_LIBS)
