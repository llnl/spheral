include_directories(.)
set(ArtificialViscosity_inst
  ArtificialViscosity
  PythonScalarArtificialViscosity
  CullenDehnenViscosity
  IncrementCullenMultipliers
  MorrisMonaghanReducingViscosity
  MonaghanGingoldViscosity
  LimitedMonaghanGingoldViscosity
  FiniteVolumeViscosity
  TensorMonaghanGingoldViscosity
  TensorSVPHViscosity
  TensorCRKSPHViscosity)


# if(ENABLE_2D)
#   list(APPEND ArtificialViscosity_sources
#        MonaghanGingoldViscosityGSRZ.cc
#       )
# endif()

instantiate(ArtificialViscosity_inst ArtificialViscosity_sources)

set(ArtificialViscosity_headers
  ArtificialViscosity.hh
  PythonScalarArtificialViscosity.hh
  ArtificialViscosityBase.hh
  ArtificialViscosityBaseInline.hh
  ArtificialViscosityView.hh
  LimitedMonaghanGingoldViscosity.hh
  LimitedMonaghanGingoldViscosityView.hh
  LimitedMonaghanGingoldViscosityViewInline.hh
  CullenDehnenViscosity.hh
  FiniteVolumeViscosity.hh
  FiniteVolumeViscosityView.hh
  FiniteVolumeViscosityViewInline.hh
  IncrementCullenMultipliers.hh
  MonaghanGingoldViscosity.hh
  MonaghanGingoldViscosityView.hh
  MonaghanGingoldViscosityViewInline.hh
  MorrisMonaghanReducingViscosity.hh
  TensorCRKSPHViscosity.hh
  TensorMonaghanGingoldViscosity.hh
  TensorMonaghanGingoldViscosityView.hh
  TensorMonaghanGingoldViscosityViewInline.hh
  TensorSVPHViscosity.hh)

# The following are currently not compiled
# MonaghanGingoldViscosityGSRZ.hh
# MonaghanGingoldViscosityRZ.hh
# LimitedMonaghanGingoldViscosityRZ.hh

spheral_add_obj_library(ArtificialViscosity SPHERAL_OBJ_LIBS)

