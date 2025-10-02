include_directories(.)

set(DEM_sources
    DEMFieldNames.cc
    ContactStorageLocation.cc
    DEMDimension.cc)

if(ENABLE_INSTANTIATIONS)
  list(APPEND DEM_sources
       DEMBaseInst.cc
       SolidBoundary/SolidBoundaryBaseInst.cc
       SolidBoundary/DEMBoundaryPolicyInst.cc
       SolidBoundary/InfinitePlaneSolidBoundaryInst.cc
       SolidBoundary/CircularPlaneSolidBoundaryInst.cc
       SolidBoundary/RectangularPlaneSolidBoundaryInst.cc
       SolidBoundary/CylinderSolidBoundaryInst.cc
       SolidBoundary/SphereSolidBoundaryInst.cc
       SolidBoundary/ClippedSphereSolidBoundaryInst.cc
       ReplacePairFieldListInst.cc
       IncrementPairFieldListInst.cc
       ReplaceAndIncrementPairFieldListInst.cc
       LinearSpringDEMInst.cc)
endif()

set(DEM_headers
    SolidBoundary/SolidBoundaryBase.hh
    SolidBoundary/DEMBoundaryPolicy.hh
    SolidBoundary/InfinitePlaneSolidBoundary.hh
    SolidBoundary/CircularPlaneSolidBoundary.hh
    SolidBoundary/RectangularPlaneSolidBoundary.hh
    SolidBoundary/CylinderSolidBoundary.hh
    SolidBoundary/SphereSolidBoundary.hh
    SolidBoundary/ClippedSphereSolidBoundary.hh
    ReplacePairFieldList.hh
    ReplaceAndIncrementPairFieldList.hh
    IncrementPairFieldList.hh
    LinearSpringDEM.hh
    DEMBase.hh
    DEMDimension.hh
    DEMFieldNames.hh
    ContactStorageLocation.hh
    setUniqueNodeIDs.hh
    setUniqueNodeIDsInline.hh)

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND DEM_headers
     DEMBase.cc
     SolidBoundary/SolidBoundaryBase.cc
     SolidBoundary/DEMBoundaryPolicy.cc
     SolidBoundary/InfinitePlaneSolidBoundary.cc
     SolidBoundary/CircularPlaneSolidBoundary.cc
     SolidBoundary/RectangularPlaneSolidBoundary.cc
     SolidBoundary/CylinderSolidBoundary.cc
     SolidBoundary/SphereSolidBoundary.cc
     SolidBoundary/ClippedSphereSolidBoundary.cc
     ReplacePairFieldList.cc
     IncrementPairFieldList.cc
     ReplaceAndIncrementPairFieldList.cc
     LinearSpringDEM.cc)

spheral_install_python_files(DEM.py)

spheral_add_obj_library(DEM SPHERAL_OBJ_LIBS)
