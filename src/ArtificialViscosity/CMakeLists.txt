include_directories(.)
set(ArtificialViscosity_inst
    ArtificialViscosityHandle
    CullenDehnenViscosity
    IncrementCullenMultipliers
    MorrisMonaghanReducingViscosity
    MonaghanGingoldViscosity
    LimitedMonaghanGingoldViscosity
    FiniteVolumeViscosity
    TensorMonaghanGingoldViscosity
    TensorSVPHViscosity
    TensorCRKSPHViscosity
   )


set(ArtificialViscosity_sources )

# if(ENABLE_2D)
#   list(APPEND ArtificialViscosity_sources
#        MonaghanGingoldViscosityGSRZ.cc
#       )
# endif()

instantiate(ArtificialViscosity_inst ArtificialViscosity_sources)

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

spheral_add_obj_library(ArtificialViscosity SPHERAL_OBJ_LIBS)

