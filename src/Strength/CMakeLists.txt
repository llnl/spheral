include_directories(.)
set(Strength_inst
    BulkModulusPolicy
    DeviatoricStressPolicy
    MeltEnergyPolicy
    PlasticStrainPolicy
    ShearModulusPolicy
    YieldStrengthPolicy
   )

set(Strength_sources )

instantiate(Strength_inst Strength_sources)

set(Strength_headers
    BulkModulusPolicy.hh
    DeviatoricStressPolicy.hh
    MeltEnergyPolicy.hh
    PlasticStrainPolicy.hh
    ShearModulusPolicy.hh
    SolidFieldNames.hh
    YieldStrengthPolicy.hh
    effectiveKernelVolume.hh
    )

spheral_add_obj_library(Strength SPHERAL_OBJ_LIBS)

