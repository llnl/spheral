set(spheralversion "v${CMAKE_PROJECT_VERSION}")
if (EXISTS ${CMAKE_SOURCE_DIR}/.git)
  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git describe --tags --abbrev=0
    OUTPUT_VARIABLE tagcomm
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git rev-parse --short HEAD
    OUTPUT_VARIABLE githash
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE abbref
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git branch --show-current
    OUTPUT_VARIABLE gitbranch
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  set(SPHERAL_VERSION_LIST ${tagcomm} ${githash} ${abbref})
  string(REPLACE ";" " " spheralversion "${SPHERAL_VERSION_LIST}")
endif()
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Spheral_banner.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/Spheral_banner.py)

# Configure SpheralConfigs.py.in
set(SPHERAL_COMP_CONFIG)
blt_list_append( TO SPHERAL_COMP_CONFIG ELEMENTS "LEOS" IF ENABLE_LEOS)
blt_list_append( TO SPHERAL_COMP_CONFIG ELEMENTS "ANEOS" IF ENABLE_ANEOS)
blt_list_append( TO SPHERAL_COMP_CONFIG ELEMENTS "OPENSUBDIV" IF ENABLE_OPENSUBDIV)
blt_list_append( TO SPHERAL_COMP_CONFIG ELEMENTS "TIMERS" IF SPHERAL_ENABLE_TIMERS)
blt_list_append( TO SPHERAL_COMP_CONFIG ELEMENTS "MPI" IF ENABLE_MPI)
blt_list_append( TO SPHERAL_COMP_CONFIG ELEMENTS "SHARED" IF ENABLE_SHARED)
set(SPHERAL_HYDRO_IMPORTS "SPHHydros;PSPHHydros;CRKSPHHydros;DEM;SPHUtilities")
blt_list_append( TO SPHERAL_HYDRO_IMPORTS ELEMENTS "FSISPHHydros;SlideSurfaces" IF SPHERAL_ENABLE_FSISPH)
blt_list_append( TO SPHERAL_HYDRO_IMPORTS ELEMENTS "GSPHHydros" IF SPHERAL_ENABLE_GSPH)
blt_list_append( TO SPHERAL_HYDRO_IMPORTS ELEMENTS "SVPHHydros" IF SPHERAL_ENABLE_SVPH)

set(SPHERAL_DIMS 1)
blt_list_append( TO SPHERAL_DIMS ELEMENTS 2 IF ENABLE_2D)
blt_list_append( TO SPHERAL_DIMS ELEMENTS 3 IF ENABLE_3D)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/SpheralConfigs.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralConfigs.py)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/spheralDimensions.py
  ${CMAKE_CURRENT_BINARY_DIR}/spheralDimensions.py)

spheral_install_python_files(
  Spheral.py
  Spheral1d.py
  Spheral2d.py
  Spheral3d.py
  SphericalSpheral.py
  SpheralRZ.py
  SolidSpheral1d.py
  SolidSpheral2d.py
  SolidSpheral3d.py
  SolidSpheralRZ.py
  SpheralTestUtilities.py
  SpheralTimer.py
  SpheralConservation.py
  SpheralController.py
  SpheralOptionParser.py
  SpheralVisitDump.py
  SpheralVoronoiSiloDump.py
  SpheralPolytopeSiloDump.py
  Spheral1dVizDump.py
  SpheralMatplotlib.py
  SpheralTimingParser.py
  findLastRestart.py
  Pnorm.py
  filearraycmp.py
  numericalIntegration.py
  newtonRaphson.py
  bisectFunction.py
  NodeHistory.py
  AverageStrain.py
  EulerianTracerHistory.py
  SpheralHadesDump.py
  siloPointmeshDump.py
  SpheralPointmeshSiloDump.py
  resampleNodeList.py
  overlayNodeList.py
  writeSiloQuadMesh.py
  buildOrderedDict.py
  CaptureStdout.py
  ${CMAKE_CURRENT_BINARY_DIR}/Spheral_banner.py
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralConfigs.py
  ${CMAKE_CURRENT_BINARY_DIR}/spheralDimensions.py
  )
