include_directories(.)
set(RK_inst
    computeRKSumVolume
    computeOccupancyVolume
    computeHullVolumes
    computeHullVolume
    computeHVolumes
    HVolumePolicy
    ContinuityVolumePolicy
    computeRKVolumes
    RKCorrections
    RKUtilities012
    RKUtilities3
    RKUtilities4
    RKUtilities5
    RKUtilities6
    RKUtilities7
    ReproducingKernelMethods
    ReproducingKernel
    interpolateRK
    gradientRK
    hessianRK
    )

set(RK_sources )

if(SPHERAL_ENABLE_2D)
  list(APPEND RK_sources
    ContinuityVolumePolicyRZ.cc
    )
endif()

instantiate(RK_inst RK_sources)

set(RK_headers
    ContinuityVolumePolicy.hh
    ContinuityVolumePolicyRZ.hh
    HVolumePolicy.hh
    RKCoefficients.hh
    RKCorrectionParams.hh
    RKCorrections.hh
    RKFieldNames.hh
    RKUtilities.hh
    RKUtilitiesInline.hh
    ReproducingKernel.hh
    ReproducingKernelInline.hh
    ReproducingKernelMethods.hh
    ReproducingKernelMethodsInline.hh
    computeHVolumes.hh
    computeHullVolumes.hh
    computeOccupancyVolume.hh
    computeRKSumVolume.hh
    computeRKVolumes.hh
    gradientRK.hh
    hessianRK.hh
    interpolateRK.hh
    )

spheral_add_obj_library(RK SPHERAL_OBJ_LIBS)

