include_directories(.)
set(Boundary_inst
    Boundary
    CRKSPHVoidBoundary
    PlanarBoundary
    ReflectingBoundary
    RigidBoundary
    PeriodicBoundary
    ConstantVelocityBoundary
    ConstantXVelocityBoundary
    ConstantYVelocityBoundary
    ConstantRVelocityBoundary
    ConstantBoundary
    findNodesTouchingThroughPlanes
    InflowOutflowBoundary
   )


set(Boundary_sources 
  FacetedVolumeBoundary.cc
  )

if(SPHERAL_ENABLE_1D)
  list(APPEND Boundary_sources
       SphericalOriginBoundary.cc
       )
endif()

if(SPHERAL_ENABLE_2D)
  list(APPEND Boundary_sources
       AxisBoundaryRZ.cc
       )
endif()

if(SPHERAL_ENABLE_3D)
  list(APPEND Boundary_sources
      CylindricalBoundary.cc  
      SphericalBoundary.cc
      ConstantZVelocityBoundaryInst.cc
      )
endif()

instantiate(Boundary_inst Boundary_sources)

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

spheral_add_obj_library(Boundary SPHERAL_OBJ_LIBS)

