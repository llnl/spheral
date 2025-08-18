include_directories(.)

set(FSISPH_inst
    SolidFSISPH
    computeFSISPHSumMassDensity
    computeHWeightedFSISPHSumMassDensity
    computeInterfacePressureCorrectedSumMassDensity
    SlideSurface
    )


set(FSISPH_sources
    SolidFSISPHRZ.cc
    FSIFieldNames.cc )


instantiate(FSISPH_inst FSISPH_sources)

set(FSISPH_headers
    SolidFSISPH.hh
    SolidFSISPHInline.hh
    SolidFSISPHRZ.hh
    computeFSISPHSumMassDensity.hh
    computeHWeightedFSISPHSumMassDensity.hh
    computeInterfacePressureCorrectedSumMassDensity.hh
    FSIFieldNames.hh
    SlideSurface.hh
    SlideSurfaceInline.hh
    )

spheral_install_python_files(
  FSISPHHydros.py
  SlideSurfaces.py
  )

spheral_add_obj_library(FSISPH SPHERAL_OBJ_LIBS)

