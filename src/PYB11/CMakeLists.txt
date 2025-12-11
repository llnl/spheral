# When Python targets are created, they add the module name to SPHERAL_MODULE_LIST
set (_python_packages
  CXXTypes
  CXXTests
  Geometry
  Silo
  DataOutput
  NodeList
  FieldView
  FieldListView
  Field
  FieldList
  Kernel
  Neighbor
  Material
  FileIO
  Utilities
  RK
  VoronoiCells
  DataBase
  Boundary
  Physics
  Hydro
  Integrator
  NodeGenerators
  FieldOperations
  SPH
  CRKSPH
  DEM
  ArtificialViscosity
  Mesh
  Damage
  SolidMaterial
  Strength
  Porosity
  KernelIntegrator
  SmoothingScale
  polytope
  OpenMP)

if(ENABLE_MPI)
  list(APPEND _python_packages Distributed)
endif()

if (SPHERAL_ENABLE_ANEOS)
  list(APPEND _python_packages ANEOS)
endif()

if (SPHERAL_ENABLE_HELMHOLTZ)
  list(APPEND _python_packages Helmholtz)
endif()

if (SPHERAL_ENABLE_ARTIFICIAL_CONDUCTION)
  list(APPEND _python_packages ArtificialConduction)
endif()

if (SPHERAL_ENABLE_EXTERNAL_FORCE)
  list(APPEND _python_packages ExternalForce)
endif()

if (SPHERAL_ENABLE_FSISPH)
  list(APPEND _python_packages FSISPH)
endif()

if (SPHERAL_ENABLE_GRAVITY)
  list(APPEND _python_packages Gravity)
endif()

if (SPHERAL_ENABLE_GSPH)
  list(APPEND _python_packages GSPH)
endif()

if (SPHERAL_ENABLE_SVPH)
  list(APPEND _python_packages SVPH)
endif()

if (SPHERAL_ENABLE_SUNDIALS)
  list(PREPEND _python_packages Solvers)
endif()

if (SPHERAL_ENABLE_LEOS)
  list(APPEND _python_packages LEOS)
endif()

foreach(_python_package ${_python_packages})
  add_subdirectory(${_python_package})
endforeach()

# Add the master module that includes all our compiled submodules
get_property(SPHERAL_MODULE_LIST GLOBAL PROPERTY SPHERAL_MODULE_LIST)
set(CXX_MODULE_LIBS ${SPHERAL_MODULE_LIST})
list(TRANSFORM CXX_MODULE_LIBS PREPEND "Spheral")

spheral_add_pybind11_library(CompiledModules DUMMY_MODULE_LIST
  DEPENDS        ${CXX_MODULE_LIBS}
  IS_SUBMODULE   OFF
  SUBMODULES     ${CXX_MODULE_LIBS})

# This global list is filled in each spheral_add_pybind11_library call
string(REPLACE ";" " " PYTHONPKGS "${SPHERAL_MODULE_LIST}")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/SpheralCompiledPackages.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralCompiledPackages.py)

set(PYB11_python_targets
  SpheralPickle.py
  PolyhedronFileUtilities.py
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralCompiledPackages.py)

spheral_install_python_files(${PYB11_python_targets})

install(FILES       Spheral.pth
        DESTINATION ${SPHERAL_SITE_PACKAGES_PATH})
