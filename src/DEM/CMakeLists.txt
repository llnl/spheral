include_directories(.)

set(DEM_inst 
    DEMBase
    SolidBoundary/SolidBoundaryBase
    SolidBoundary/DEMBoundaryPolicy
    SolidBoundary/InfinitePlaneSolidBoundary
    SolidBoundary/CircularPlaneSolidBoundary
    SolidBoundary/RectangularPlaneSolidBoundary
    SolidBoundary/CylinderSolidBoundary
    SolidBoundary/SphereSolidBoundary
    SolidBoundary/ClippedSphereSolidBoundary
    ReplacePairFieldList
    IncrementPairFieldList
    ReplaceAndIncrementPairFieldList
    LinearSpringDEM)

set(DEM_sources
    ContactStorageLocation.cc
    DEMDimension.cc)

instantiate(DEM_inst DEM_sources)

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

add_subdirectory(SolidBoundary)

spheral_install_python_files(DEM.py)

spheral_add_obj_library(DEM SPHERAL_OBJ_LIBS)

