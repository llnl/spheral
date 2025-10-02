include_directories(.)

set(ArtificialConduction_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND ArtificialConduction_sources
       ArtificialConductionInst.cc
       ArtificialConductionPolicyInst.cc
       )
endif()

set(ArtificialConduction_headers
    ArtificialConduction.hh
    ArtificialConductionPolicy.hh
    ) 

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND ArtificialConduction_headers
     ArtificialConduction.cc
     ArtificialConductionPolicy.cc
     )

spheral_add_obj_library(ArtificialConduction SPHERAL_OBJ_LIBS)

