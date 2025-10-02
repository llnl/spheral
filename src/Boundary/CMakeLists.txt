include_directories(.)

set(Boundary_sources 
    FacetedVolumeBoundary.cc
    )

if(ENABLE_1D)
  list(APPEND Boundary_sources
       SphericalOriginBoundary.cc
       )
endif()

if(ENABLE_2D)
  list(APPEND Boundary_sources
       AxisBoundaryRZ.cc
       )
endif()

if(ENABLE_3D)
  list(APPEND Boundary_sources
       CylindricalBoundary.cc  
       SphericalBoundary.cc
       ConstantZVelocityBoundaryInst.cc
       )
endif()

if(ENABLE_INSTANTIATIONS)
  list(APPEND Boundary_sources
       BoundaryInst.cc
       CRKSPHVoidBoundaryInst.cc
       PlanarBoundaryInst.cc
       ReflectingBoundaryInst.cc
       RigidBoundaryInst.cc
       PeriodicBoundaryInst.cc
       ConstantVelocityBoundaryInst.cc
       ConstantXVelocityBoundaryInst.cc
       ConstantYVelocityBoundaryInst.cc
       ConstantRVelocityBoundaryInst.cc
       ConstantBoundaryInst.cc
       findNodesTouchingThroughPlanesInst.cc
       InflowOutflowBoundaryInst.cc
       )
endif()

set(Boundary_headers
    AxisBoundaryRZ.hh
    SphericalOriginBoundary.hh
    Boundary.hh
    BoundaryInline.hh
    CRKSPHVoidBoundary.hh
    ConstantBoundary.hh
    ConstantBoundaryInline.hh
    ConstantBoundaryUtilities.hh
    ConstantRVelocityBoundary.hh
    ConstantVelocityBoundary.hh
    ConstantVelocityBoundaryInline.hh
    ConstantXVelocityBoundary.hh
    ConstantYVelocityBoundary.hh
    ConstantZVelocityBoundary.hh
    CylindricalBoundary.hh
    CylindricalBoundaryInline.hh
    FacetedVolumeBoundary.hh
    FacetedVolumeBoundaryInline.hh
    InflowOutflowBoundary.hh
    InflowOutflowBoundaryInline.hh
    PeriodicBoundary.hh
    PeriodicPlanarBoundary.hh
    PlanarBoundary.hh
    PlanarBoundaryInline.hh
    PlanarSymmetryBoundary.hh
    ReflectingBoundary.hh
    ReflectingBoundaryInline.hh
    RigidBoundary.hh
    SphericalBoundary.hh
    SphericalBoundaryInline.hh
    ThreeOneDimensionBoundary.hh
    findNodesTouchingThroughPlanes.hh
    mapPositionThroughPlanes.hh
    )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND Boundary_headers
     Boundary.cc
     CRKSPHVoidBoundary.cc
     PlanarBoundary.cc
     ReflectingBoundary.cc
     RigidBoundary.cc
     PeriodicBoundary.cc
     ConstantVelocityBoundary.cc
     ConstantXVelocityBoundary.cc
     ConstantYVelocityBoundary.cc
     ConstantRVelocityBoundary.cc
     ConstantBoundary.cc
     findNodesTouchingThroughPlanes.cc
     InflowOutflowBoundary.cc
     )

if(ENABLE_3D)
  list(APPEND Boundary_headers
       ConstantZVelocityBoundary.cc
       )
endif()

spheral_add_obj_library(Boundary SPHERAL_OBJ_LIBS)

