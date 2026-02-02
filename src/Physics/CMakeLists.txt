include_directories(.)
set(Physics_inst
    GenericBodyForce
    GenericHydro
    Physics
   )

set(Physics_sources )

instantiate(Physics_inst Physics_sources)

set(Physics_headers
    GenericBodyForce.hh
    GenericHydro.hh
    Physics.hh
    )

spheral_add_obj_library(Physics SPHERAL_OBJ_LIBS)

