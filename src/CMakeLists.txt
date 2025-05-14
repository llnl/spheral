################################
# Build spheral packages 
################################
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/InstantiateCXX.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralAddLibs.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralInstallPythonFiles.cmake)

include_directories(.)
include_directories(${PROJECT_BINARY_DIR}/src)

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

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

if (ENABLE_SUNDIALS)
  list(APPEND _packages Solvers)
endif()

if (ENABLE_LEOS)
  list(APPEND _packages LEOS)
endif()

if(NOT ENABLE_CXXONLY)
  list(APPEND _packages
    SimulationControl
    PYB11)
endif()

configure_file(config.hh.in
  ${PROJECT_BINARY_DIR}/src/config.hh)

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()

# For dev builds, we must call LLNLSpheralInstallObjs.cmake
# to ensure LLNLSpheral libraries are added to the Spheral_CXX
# target
if(NOT ENABLE_DEV_BUILD)
  set(CXX_sources spheralCXX.cc)
endif()
# Retrieve the global list populated in spheral_obj_add_library
get_property(SPHERAL_OBJ_LIBS GLOBAL PROPERTY SPHERAL_OBJ_LIBS)
# Must use quotes when passing lists as inputs for functions
spheral_add_cxx_library(CXX "${SPHERAL_OBJ_LIBS}")

# This calls LLNLSpheralInstallObjs.cmake
if(EXISTS ${EXTERNAL_SPHERAL_OBJ_CMAKE})
  include(${EXTERNAL_SPHERAL_OBJ_CMAKE})
endif()
