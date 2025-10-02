include_directories(.)

set(Gravity_sources
    ApproximatePolyhedralGravityModel.cc
    TreeInst.cc
    TreeGravityInst.cc
    PolyGravityInst.cc
    )

if(ENABLE_INSTANTIATIONS)
  list(APPEND Gravity_sources
       NBodyGravityInst.cc
       )
endif()

set(Gravity_headers
    CompatibleGravitationalVelocityPolicy.hh
    NBodyGravity.hh
    Tree.hh
    TreeGravity.hh
    TreeGravityInline.hh
    ApproximatePolyhedralGravityModel.hh
    ApproximatePolyhedralGravityModelInline.hh
    PolyGravity.hh
    )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND Gravity_headers
     NBodyGravity.cc
     )

spheral_add_obj_library(Gravity SPHERAL_OBJ_LIBS)

