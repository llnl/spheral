include_directories(.)
set(DataOutput_sources
    RestartRegistrar.cc
    #RestartableObject.cc
    )

set(DataOutput_headers
    Restart.hh
    RestartHandle.hh
    RestartInline.hh
    RestartRegistrar.hh
    RestartRegistrarInline.hh
    RestartableObject.hh
    registerWithRestart.hh
    )

spheral_add_obj_library(DataOutput)

