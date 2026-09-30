################################
# Build spheral packages 
################################
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/InstantiateCXX.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralAddLibs.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralInstallPythonFiles.cmake)

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

# Set the name of the current package we are creating
set(SPHERAL_ACTIVE_CXX_TARGET Spheral_CXX)
# Create the target for this overarching package
spheral_initialize_cxx_target("spheralCXX.cc")

# These are internal Spheral packages
list(APPEND _packages
  ArtificialViscosity
  Boundary
  CRKSPH
  CXXTests
  DEM
  Damage
  DataBase
  DataOutput
  Distributed
  Threading
  Field
  FieldOperations
  FileIO
  Geometry
  Hydro
  Integrator
  Kernel
  KernelIntegrator
  Material
  Mesh
  Neighbor
  NodeGenerators
  NodeList
  Physics
  Porosity
  RK
  VoronoiCells
  SPH
  SolidMaterial
  Strength
  Utilities
  SmoothingScale
)

if (SPHERAL_ENABLE_ARTIFICIAL_CONDUCTION)
  list(APPEND _packages ArtificialConduction)
endif()

if (SPHERAL_ENABLE_EXTERNAL_FORCE)
  list(APPEND _packages ExternalForce)
endif()

if (SPHERAL_ENABLE_FSISPH)
  list(APPEND _packages FSISPH)
endif()

if (SPHERAL_ENABLE_GRAVITY)
  list(APPEND _packages Gravity)
endif()

if (SPHERAL_ENABLE_GSPH)
  list(APPEND _packages GSPH)
endif()

if (SPHERAL_ENABLE_SVPH)
  list(APPEND _packages SVPH)
endif()

if (SPHERAL_ENABLE_SUNDIALS)
  list(APPEND _packages Solvers)
endif()

if (SPHERAL_ENABLE_LEOS)
  list(APPEND _packages LEOS)
endif()

if(SPHERAL_ENABLE_PYTHON)
  list(APPEND _packages
    SimulationControl
    PYB11)
endif()

configure_file(config.hh.in
  ${PROJECT_BINARY_DIR}/src/config.hh)
install(FILES ${PROJECT_BINARY_DIR}/src/config.hh
  DESTINATION include)

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()

# This calls LLNLSpheralInstallObjs.cmake
if(EXISTS ${EXTERNAL_SPHERAL_OBJ_CMAKE})
  include(${EXTERNAL_SPHERAL_OBJ_CMAKE})
endif()

spheral_install_cxx_library()
