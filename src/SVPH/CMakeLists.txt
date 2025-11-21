include_directories(.)
set(SVPH_inst
    computeSVPHCorrections
    computeSVPHCorrectionsOnFaces
    computeSumVoronoiCellMassDensityFromFaces
    sampleFieldListSVPH
    gradientFieldListSVPH
    SVPHFacetedHydroBase
    SVPHCorrectionsPolicy
    SVPHMassDensityPolicy
    CompatibleFaceSpecificThermalEnergyPolicy
    SpecificThermalEnergyVolumePolicy
    CellPressurePolicy
   )


set(SVPH_sources )

instantiate(SVPH_inst SVPH_sources)

set(SVPH_headers
    CellPressurePolicy.hh
    CompatibleFaceSpecificThermalEnergyPolicy.hh
    SVPHCorrectionsPolicy.hh
    SVPHFacetedHydroBase.hh
    SVPHFacetedHydroBaseInline.hh
    SVPHFieldNames.hh
    SVPHHydroBase.hh
    SVPHHydroBaseInline.hh
    SVPHMassDensityPolicy.hh
    SpecificThermalEnergyVolumePolicy.hh
    computeSVPHCorrections.hh
    computeSVPHCorrectionsOnFaces.hh
    computeSumVoronoiCellMassDensityFromFaces.hh
    gradientFieldListSVPH.hh
    sampleFieldListSVPH.hh
    )

spheral_install_python_files(SVPHHydros.py)

spheral_add_obj_library(SVPH SPHERAL_OBJ_LIBS)

