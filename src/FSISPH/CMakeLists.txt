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

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND FSISPH_headers
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

