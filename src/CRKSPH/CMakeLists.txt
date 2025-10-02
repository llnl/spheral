include_directories(.)

set(CRKSPH_sources
    centerOfMass.cc
   )

if(ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHRZ.cc
       SolidCRKSPHRZ.cc
      )
endif()

if(ENABLE_INSTANTIATIONS)
  list(APPEND CRKSPH_sources
       computeCRKSPHSumMassDensityInst.cc
       computeSolidCRKSPHSumMassDensityInst.cc
       detectSurfaceInst.cc
       editMultimaterialSurfaceTopologyInst.cc
       zerothOrderSurfaceCorrectionsInst.cc
       CRKSPHBaseInst.cc
       CRKSPHInst.cc
       SolidCRKSPHInst.cc
      )
endif()

set(CRKSPH_headers
    CRKSPHBase.hh
    CRKSPH.hh
    CRKSPHRZ.hh
    CRKSPHVariant.hh
    SolidCRKSPH.hh
    SolidCRKSPHRZ.hh
    centerOfMass.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHSumMassDensity.hh
    computeHullSumMassDensity.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    detectSurface.hh
    editMultimaterialSurfaceTopology.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND Boundary_headers
     computeCRKSPHSumMassDensity.cc
     computeSolidCRKSPHSumMassDensity.cc
     detectSurface.cc
     editMultimaterialSurfaceTopology.cc
     zerothOrderSurfaceCorrections.cc
     CRKSPHBase.cc
     CRKSPH.cc
     SolidCRKSPH.cc
     )

spheral_add_obj_library(CRKSPH SPHERAL_OBJ_LIBS)

spheral_install_python_files(CRKSPHHydros.py)

