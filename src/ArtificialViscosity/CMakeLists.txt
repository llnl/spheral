include_directories(.)

set(ArtificialViscosity_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND ArtificialViscosity_sources
       ArtificialViscosityHandleInst.cc
       CullenDehnenViscosityInst.cc
       IncrementCullenMultipliersInst.cc
       MorrisMonaghanReducingViscosityInst.cc
       MonaghanGingoldViscosityInst.cc
       LimitedMonaghanGingoldViscosityInst.cc
       FiniteVolumeViscosityInst.cc
       TensorMonaghanGingoldViscosityInst.cc
       TensorSVPHViscosityInst.cc
       TensorCRKSPHViscosityInst.cc
       )
endif()

# if(ENABLE_2D)
#   list(APPEND ArtificialViscosity_sources
#        MonaghanGingoldViscosityGSRZ.cc
#       )
# endif()

set(ArtificialViscosity_headers
    ArtificialViscosityHandle.hh
    ArtificialViscosityHandleInline.hh
    ArtificialViscosity.hh
    LimitedMonaghanGingoldViscosity.hh
    LimitedMonaghanGingoldViscosityRZ.hh
    CullenDehnenViscosity.hh
    FiniteVolumeViscosity.hh
    IncrementCullenMultipliers.hh
    MonaghanGingoldViscosity.hh
    MonaghanGingoldViscosityGSRZ.hh
    MonaghanGingoldViscosityRZ.hh
    MorrisMonaghanReducingViscosity.hh
    TensorCRKSPHViscosity.hh
    TensorMonaghanGingoldViscosity.hh
    TensorSVPHViscosity.hh
    )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND ArtificialViscosity_headers
     ArtificialViscosityHandle.cc
     CullenDehnenViscosity.cc
     IncrementCullenMultipliers.cc
     MorrisMonaghanReducingViscosity.cc
     MonaghanGingoldViscosity.cc
     LimitedMonaghanGingoldViscosity.cc
     FiniteVolumeViscosity.cc
     TensorMonaghanGingoldViscosity.cc
     TensorSVPHViscosity.cc
     TensorCRKSPHViscosity.cc
     )

spheral_add_obj_library(ArtificialViscosity SPHERAL_OBJ_LIBS)

