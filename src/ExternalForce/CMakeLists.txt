include_directories(.)

set(ExternalForce_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND ExternalForce_sources
       ConstantAccelerationInst.cc
       LinearAccelerationInst.cc
       NFWPotentialInst.cc
       PointPotentialInst.cc
       )
endif()

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

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND ExternalForce_headers
     ConstantAcceleration.cc
     LinearAcceleration.cc
     NFWPotential.cc
     PointPotential.cc
     )

set(ExternalForce_includes
   ${SPHERAL_SOURCE_DIR}
   ${SPHERAL_SOURCE_DIR}/ExternalForce
   ${EIGEN_LIB_LOCATION}
   )

spheral_add_obj_library(ExternalForce SPHERAL_OBJ_LIBS)

