include_directories(.)
set(Hydro_inst
    EntropyPolicy
    CompatibleDifferenceSpecificThermalEnergyPolicy
    NonSymmetricSpecificThermalEnergyPolicy
    PressurePolicy
    SoundSpeedPolicy
    SpecificFromTotalThermalEnergyPolicy
    SpecificThermalEnergyPolicy
    SumVoronoiMassDensityPolicy
    VolumePolicy
    VoronoiMassDensityPolicy
    GammaPolicy
    computeSPHVolume
   )


set(Hydro_sources )

if (SPHERAL_ENABLE_2D)
  list(APPEND Hydro_sources RZNonSymmetricSpecificThermalEnergyPolicy.cc)
endif()

if (SPHERAL_ENABLE_1D AND SPHERAL_ENABLE_2D)
  list(APPEND Hydro_sources SphericalPositionPolicy.cc)
endif()

instantiate(Hydro_inst Hydro_sources)

set(Hydro_headers
    EntropyPolicy.hh
    GammaPolicy.hh
    HydroFieldNames.hh
    CompatibleDifferenceSpecificThermalEnergyPolicy.hh
    NonSymmetricSpecificThermalEnergyPolicy.hh
    PressurePolicy.hh
    RZNonSymmetricSpecificThermalEnergyPolicy.hh
    SoundSpeedPolicy.hh
    SpecificFromTotalThermalEnergyPolicy.hh
    SpecificThermalEnergyPolicy.hh
    SumVoronoiMassDensityPolicy.hh
    VolumePolicy.hh
    VoronoiMassDensityPolicy.hh
    computeSPHVolume.hh
    entropyWeightingFunction.hh
    )

spheral_add_obj_library(Hydro SPHERAL_OBJ_LIBS)

