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

spheral_add_package(ExternalForce)

