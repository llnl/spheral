include_directories(.)
set(ExternalForce_inst
    ConstantAcceleration
    LinearAcceleration
    NFWPotential
    PointPotential
    )

set(ExternalForce_sources )
instantiate(ExternalForce_inst ExternalForce_sources)

set(ExternalForce_headers
    ConstantAcceleration.hh
    ConstantAccelerationInline.hh
    LinearAcceleration.hh
    LinearAccelerationInline.hh
    NFWPotential.hh
    NFWPotentialInline.hh
    PointPotential.hh
    PointPotentialInline.hh
    )

set(ExternalForce_includes
   ${SPHERAL_SOURCE_DIR}
   ${SPHERAL_SOURCE_DIR}/ExternalForce
   ${EIGEN_LIB_LOCATION}
   )


spheral_add_obj_library(ExternalForce)

