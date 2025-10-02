include_directories(.)

set(GSPH_sources 
    GSPHFieldNames.cc)

if(ENABLE_INSTANTIATIONS)
  list(APPEND GSPH_sources
       computeSumVolumeInst.cc
       computeMFMDensityInst.cc
       initializeGradientsInst.cc
       GenericRiemannHydroInst.cc
       GSPHInst.cc
       MFMInst.cc
       MFVInst.cc
       Policies/MassFluxPolicyInst.cc
       Policies/MFVIncrementVelocityPolicyInst.cc
       Policies/MFVIncrementSpecificThermalEnergyPolicyInst.cc
       Policies/CompatibleMFVSpecificThermalEnergyPolicyInst.cc
       Limiters/LimiterBaseInst.cc
       Limiters/VanLeerLimiterInst.cc
       Limiters/SuperbeeLimiterInst.cc
       Limiters/MinModLimiterInst.cc
       Limiters/VanAlbaLimiterInst.cc
       Limiters/OspreLimiterInst.cc
       Limiters/BarthJespersenLimiterInst.cc
       WaveSpeeds/WaveSpeedBaseInst.cc
       WaveSpeeds/AcousticWaveSpeedInst.cc
       WaveSpeeds/DavisWaveSpeedInst.cc
       WaveSpeeds/EinfeldtWaveSpeedInst.cc
       RiemannSolvers/RiemannSolverBaseInst.cc
       RiemannSolvers/HLLCInst.cc
       RiemannSolvers/SecondOrderArtificialViscosityInst.cc)
endif()

set(GSPH_headers
    computeSumVolume.hh
    computeMFMDensity.hh
    initializeGradients.hh
    GSPHFieldNames.hh
    GenericRiemannHydro.hh
    GSPH.hh
    MFM.hh
    MFV.hh
    GenericRiemannHydroInline.hh
    GSPHInline.hh
    MFMInline.hh
    MFVInline.hh
    Policies/MassFluxPolicy.hh
    Policies/MFVIncrementVelocityPolicy.hh
    Policies/MFVIncrementSpecificThermalEnergyPolicy.hh
    Policies/CompatibleMFVSpecificThermalEnergyPolicy.hh
    Limiters/LimiterBase.hh
    Limiters/VanLeerLimiter.hh
    Limiters/SuperbeeLimiter.hh
    Limiters/MinModLimiter.hh
    Limiters/VanAlbaLimiter.hh
    Limiters/OspreLimiter.hh
    Limiters/BarthJespersenLimiter.hh
    WaveSpeeds/WaveSpeedBase.hh
    WaveSpeeds/AcousticWaveSpeed.hh
    WaveSpeeds/DavisWaveSpeed.hh
    WaveSpeeds/EinfeldtWaveSpeed.hh
    RiemannSolvers/RiemannSolverBase.hh
    RiemannSolvers/HLLC.hh
    RiemannSolvers/SecondOrderArtificialViscosity.hh)

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND GSPH_headers
     computeSumVolumeInst.cc
     computeMFMDensityInst.cc
     initializeGradientsInst.cc
     GenericRiemannHydroInst.cc
     GSPHInst.cc
     MFMInst.cc
     MFVInst.cc
     Policies/MassFluxPolicyInst.cc
     Policies/MFVIncrementVelocityPolicyInst.cc
     Policies/MFVIncrementSpecificThermalEnergyPolicyInst.cc
     Policies/CompatibleMFVSpecificThermalEnergyPolicyInst.cc
     Limiters/LimiterBaseInst.cc
     Limiters/VanLeerLimiterInst.cc
     Limiters/SuperbeeLimiterInst.cc
     Limiters/MinModLimiterInst.cc
     Limiters/VanAlbaLimiterInst.cc
     Limiters/OspreLimiterInst.cc
     Limiters/BarthJespersenLimiterInst.cc
     WaveSpeeds/WaveSpeedBaseInst.cc
     WaveSpeeds/AcousticWaveSpeedInst.cc
     WaveSpeeds/DavisWaveSpeedInst.cc
     WaveSpeeds/EinfeldtWaveSpeedInst.cc
     RiemannSolvers/RiemannSolverBaseInst.cc
     RiemannSolvers/HLLCInst.cc
     RiemannSolvers/SecondOrderArtificialViscosityInst.cc)

spheral_install_python_files(GSPHHydros.py)

spheral_add_obj_library(GSPH SPHERAL_OBJ_LIBS)
