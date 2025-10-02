include_directories(.)

set(FSISPH_sources
    FSIFieldNames.cc
    )

if(ENABLE_INSTANTIATIONS)
  list(APPEND FSISPH_sources
       SolidFSISPHInst.cc
       computeFSISPHSumMassDensityInst.cc
       computeHWeightedFSISPHSumMassDensityInst.cc
       computeInterfacePressureCorrectedSumMassDensityInst.cc
       SlideSurfaceInst.cc
       )
endif()

set(FSISPH_headers
    SolidFSISPH.hh
    SolidFSISPHInline.hh
    computeFSISPHSumMassDensity.hh
    computeHWeightedFSISPHSumMassDensity.hh
    computeInterfacePressureCorrectedSumMassDensity.hh
    FSIFieldNames.hh
    SlideSurface.hh
    SlideSurfaceInline.hh
    )

list(APPEND FSISPH_sources
     SolidFSISPH.cc
     computeFSISPHSumMassDensity.cc
     computeHWeightedFSISPHSumMassDensity.cc
     computeInterfacePressureCorrectedSumMassDensity.cc
     SlideSurface.cc
     )

spheral_install_python_files(
  FSISPHHydros.py
  SlideSurfaces.py
  )

spheral_add_obj_library(FSISPH SPHERAL_OBJ_LIBS)

