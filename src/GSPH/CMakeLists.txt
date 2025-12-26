include_directories(.)

set(GSPH_inst
    computeSumVolume
    computeMFMDensity
    initializeGradients
    GenericRiemannHydro
    GSPH
    MFM
    MFV
    Policies/MassFluxPolicy
    Policies/MFVIncrementVelocityPolicy
    Policies/MFVIncrementSpecificThermalEnergyPolicy
    Policies/CompatibleMFVSpecificThermalEnergyPolicy
    Limiters/LimiterBase
    Limiters/VanLeerLimiter
    Limiters/SuperbeeLimiter
    Limiters/MinModLimiter
    Limiters/VanAlbaLimiter
    Limiters/OspreLimiter
    Limiters/BarthJespersenLimiter
    WaveSpeeds/WaveSpeedBase
    WaveSpeeds/AcousticWaveSpeed
    WaveSpeeds/DavisWaveSpeed
    WaveSpeeds/EinfeldtWaveSpeed
    RiemannSolvers/RiemannSolverBase
    RiemannSolvers/HLLC
    RiemannSolvers/SecondOrderArtificialViscosity)

set(GSPH_sources )

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

add_subdirectory(Limiters)
add_subdirectory(WaveSpeeds)
add_subdirectory(RiemannSolvers)
add_subdirectory(Policies)

instantiate(GSPH_inst GSPH_sources)

spheral_install_python_files(GSPHHydros.py)

spheral_add_obj_library(GSPH SPHERAL_OBJ_LIBS)

