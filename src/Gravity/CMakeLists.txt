include_directories(.)
set(Gravity_inst
  NBodyGravity
  )

set(Gravity_sources
  ApproximatePolyhedralGravityModel.cc
  TreeInst.cc
  TreeGravityInst.cc
  PolyGravityInst.cc
  )
instantiate(Gravity_inst Gravity_sources)

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

spheral_add_obj_library(Gravity)

