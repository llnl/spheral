include_directories(.)
set(ArtificialViscosity_inst
  ArtificialViscosity
  CullenDehnenViscosity
  IncrementCullenMultipliers
  MorrisMonaghanReducingViscosity
  MonaghanGingoldViscosity
  LimitedMonaghanGingoldViscosity
  FiniteVolumeViscosity
  TensorMonaghanGingoldViscosity
  TensorSVPHViscosity
  TensorCRKSPHViscosity)


set(ArtificialViscosity_sources)

# if(ENABLE_2D)
#   list(APPEND ArtificialViscosity_sources
#        MonaghanGingoldViscosityGSRZ.cc
#       )
# endif()

instantiate(ArtificialViscosity_inst ArtificialViscosity_sources)

set(ArtificialViscosity_headers
  ArtificialViscosity.hh
  ArtificialViscosityBase.hh
  ArtificialViscosityBaseInline.hh
  ArtificialViscosityView.hh
  LimitedMonaghanGingoldViscosity.hh
  LimitedMonaghanGingoldViscosityView.hh
  CullenDehnenViscosity.hh
  FiniteVolumeViscosity.hh
  FiniteVolumeViscosityView.hh
  IncrementCullenMultipliers.hh
  MonaghanGingoldViscosity.hh
  MonaghanGingoldViscosityView.hh
  MorrisMonaghanReducingViscosity.hh
  TensorCRKSPHViscosity.hh
  TensorMonaghanGingoldViscosity.hh
  TensorMonaghanGingoldViscosityView.hh
  TensorSVPHViscosity.hh)

# The following are currently not compiled
# MonaghanGingoldViscosityGSRZ.hh
# MonaghanGingoldViscosityRZ.hh
# LimitedMonaghanGingoldViscosityRZ.hh

spheral_add_obj_library(ArtificialViscosity SPHERAL_OBJ_LIBS)

