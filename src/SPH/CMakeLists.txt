include_directories(.)
set(SPH_inst
  SPHBase
  SPH
  PSPH
  SolidSPH
  computeSPHOmegaGradhCorrection
  computeSPHSumMassDensity
  correctSPHSumMassDensity
  computePSPHCorrections
  computeSumVoronoiCellMassDensity
)

set(SPH_sources )

if(ENABLE_1D)
  list(APPEND SPH_sources
       SphericalSPHInst.cc
       SphericalSPH.cc
       SolidSphericalSPH.cc
      )
endif()

if(ENABLE_2D)
  list(APPEND SPH_sources
       SPHRZ.cc
       SolidSPHRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND SPH_sources
       computeHydrostaticEquilibriumPressure.cc
      )
endif()

instantiate(SPH_inst SPH_sources)

set(SPH_headers
  SPH.hh
  PSPH.hh
  SPHBase.hh
  SPHRZ.hh
  SphericalSPH.hh
  SolidSPH.hh
  SolidSPHRZ.hh
  SolidSphericalSPH.hh
  computeHydrostaticEquilibriumPressure.hh
  computePSPHCorrections.hh
  computeSPHOmegaGradhCorrection.hh
  computeSPHSumMassDensity.hh
  computeSumVoronoiCellMassDensity.hh
  correctSPHSumMassDensity.hh
)

spheral_install_python_files(
  SPHHydros.py
  PSPHHydros.py
  SPHUtilities.py
)

spheral_add_obj_library(SPH SPHERAL_OBJ_LIBS)

