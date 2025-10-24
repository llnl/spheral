include_directories(.)
set(CRKSPH_inst
    computeCRKSPHSumMassDensity
    computeSolidCRKSPHSumMassDensity
    detectSurface
    editMultimaterialSurfaceTopology
    zerothOrderSurfaceCorrections
    CRKSPHBase
    CRKSPH
    SolidCRKSPH
   )


set(CRKSPH_sources
    centerOfMass.cc
   )

if(SPHERAL_ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHRZ.cc
       SolidCRKSPHRZ.cc
      )
endif()

if(SPHERAL_ENABLE_3D)
  list(APPEND CRKSPH_sources
      )
endif()

instantiate(CRKSPH_inst CRKSPH_sources)

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

spheral_add_obj_library(CRKSPH SPHERAL_OBJ_LIBS)

spheral_install_python_files(CRKSPHHydros.py)

