include(ExternalProject)

#-------------------------------------------------------------------------------
# Configure CMake
#-------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

if (NOT SPHERAL_CMAKE_MODULE_PATH)
  set(SPHERAL_CMAKE_MODULE_PATH "${SPHERAL_ROOT_DIR}/cmake")
endif()
list(APPEND CMAKE_MODULE_PATH "${SPHERAL_CMAKE_MODULE_PATH}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#-------------------------------------------------------------------------------
# Add Spheral CMake Macros for tests and executables
#-------------------------------------------------------------------------------
include(SpheralMacros)

#-------------------------------------------------------------------------------
# Set Compiler Flags / Options
#-------------------------------------------------------------------------------
include(Compilers)

#-------------------------------------------------------------------------------
# Configure and Include blt
#-------------------------------------------------------------------------------

# Need to define Python paths here as BLT finds it's own Python package.
set(Python_EXECUTABLE ${python_DIR}/bin/python3)
set(Python3_EXECUTABLE ${python_DIR}/bin/python3)

set(ENABLE_MPI ON CACHE BOOL "")
set(ENABLE_OPENMP ON CACHE BOOL "")
set(BLT_DOCS_TARGET_NAME "blt_docs" CACHE STRING "")

if(NOT SPHERAL_BLT_DIR)
  set (SPHERAL_BLT_REL_DIR "${SPHERAL_ROOT_DIR}/cmake/blt" CACHE PATH "")
  get_filename_component(SPHERAL_BLT_DIR "${SPHERAL_BLT_REL_DIR}" ABSOLUTE)
endif()

if (NOT EXISTS "${SPHERAL_BLT_DIR}/SetupBLT.cmake")
  message(FATAL_ERROR
            "${SPHERAL_BLT_DIR} is not present.\n"
            "call cmake with -DSPHERAL_BLT_DIR=/your/installation/of/blt\n")
endif()

include(${SPHERAL_BLT_DIR}/SetupBLT.cmake)

#-------------------------------------------------------------------------------
# Include standard build system logic and options / definitions
#-------------------------------------------------------------------------------
# TODO: Prefix Spheral options
set(ENABLE_CXXONLY OFF CACHE BOOL "enable C++ only build without python bindings")
set(ENABLE_1D ON CACHE BOOL "enable 1d")
set(ENABLE_2D ON CACHE BOOL "enable 2d")
set(ENABLE_3D ON CACHE BOOL "enable 3d")
set(ENABLE_INSTANTIATIONS ON CACHE BOOL "enable instantiations")
set(SPHERAL_ENABLE_TIMERS OFF CACHE BOOL "enable timer")
set(ENABLE_ANEOS ON CACHE BOOL "enable the ANEOS equation of state package")
set(ENABLE_OPENSUBDIV ON CACHE BOOL "enable the Opensubdiv Pixar extension for refining polyhedra")
set(ENABLE_HELMHOLTZ ON CACHE BOOL "enable the Helmholtz equation of state package")

option(SPHERAL_ENABLE_ARTIFICIAL_CONDUCTION "Enable the artificial conduction package" ON)
option(SPHERAL_ENABLE_EXTERNAL_FORCE "Enable the external force package" ON)
option(SPHERAL_ENABLE_FSISPH "Enable the FSISPH package" ON)
option(SPHERAL_ENABLE_GRAVITY "Enable the gravity package" ON)
option(SPHERAL_ENABLE_GSPH "Enable the GSPH package" ON)
option(SPHERAL_ENABLE_SVPH "Enable the SVPH package" ON)
option(SPHERAL_ENABLE_GLOBALDT_REDUCTION "Enable global allreduce for the time step" ON)
option(SPHERAL_ENABLE_LONGCSDT "Enable longitudinal sound speed time step constraint" ON)

option(SPHERAL_ENABLE_LOGGER "Enable debug log printing" OFF)
option(ENABLE_DEV_BUILD "Build separate internal C++ libraries for faster code development" OFF)
option(ENABLE_STATIC_CXXONLY "build only static libs" OFF)
option(ENABLE_SHARED "Building C++ libs shared" ON)

if(ENABLE_STATIC_CXXONLY)
  set(ENABLE_CXXONLY ON)
  set(ENABLE_SHARED OFF)
endif()

if(ENABLE_MPI)
  set(BLT_MPI_COMPILE_FLAGS -DUSE_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DOMPI_SKIP_MPICXX)
  list(APPEND SPHERAL_CXX_DEPENDS mpi)
endif()

if(ENABLE_OPENMP)
  list(APPEND SPHERAL_CXX_DEPENDS openmp)
endif()

if(ENABLE_CUDA)
  set(CMAKE_CUDA_FLAGS  "${CMAKE_CUDA_FLAGS} -arch=${CUDA_ARCH} --extended-lambda -Xcudafe --display_error_number")
  #set(CMAKE_CUDA_FLAGS  "${CMAKE_CUDA_FLAGS} -arch=${CUDA_ARCH} --expt-relaxed-constexpr --extended-lambda -Xcudafe --display_error_number")
  set(CMAKE_CUDA_STANDARD 17)
  list(APPEND SPHERAL_CXX_DEPENDS cuda)
  set(SPHERAL_ENABLE_CUDA ON)
endif()

if(ENABLE_HIP)
  list(APPEND SPHERAL_CXX_DEPENDS blt::hip)
  list(APPEND SPHERAL_CXX_DEPENDS blt::hip_runtime)
  set(SPHERAL_ENABLE_HIP ON)
endif()

#-------------------------------------------------------------------------------#
# Set a default build type if none was specified
#-------------------------------------------------------------------------------#
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build (debug, release, etc)." FORCE)

  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#-------------------------------------------------------------------------------
# Should we build sphinx documentation
#-------------------------------------------------------------------------------
set(ENABLE_DOCS OFF CACHE BOOL "enable sphinx Spheral documentation")

#-------------------------------------------------------------------------------
# Locate third party libraries
#-------------------------------------------------------------------------------
include(${SPHERAL_ROOT_DIR}/cmake/InstallTPLs.cmake)

include(${SPHERAL_ROOT_DIR}/cmake/CMakeDefinitions.cmake)

#-------------------------------------------------------------------------------
# Set full rpath information by default
#-------------------------------------------------------------------------------
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#-------------------------------------------------------------------------------
# Set global variables used for dependencies
#-------------------------------------------------------------------------------
# List of external dependencies
set_property(GLOBAL PROPERTY SPHERAL_BLT_DEPENDS "${SPHERAL_BLT_DEPENDS}")
# List of compiler dependencies
set_property(GLOBAL PROPERTY SPHERAL_CXX_DEPENDS "${SPHERAL_CXX_DEPENDS}")

#-------------------------------------------------------------------------------
# Prepare to build the src
#-------------------------------------------------------------------------------
configure_file(${SPHERAL_ROOT_DIR}/src/config.hh.in
  ${PROJECT_BINARY_DIR}/src/config.hh)
include_directories(${PROJECT_BINARY_DIR}/src)

add_subdirectory(${SPHERAL_ROOT_DIR}/src)

#-------------------------------------------------------------------------------
# Add the documentation
#-------------------------------------------------------------------------------
if (NOT ENABLE_CXXONLY)
  add_subdirectory(${SPHERAL_ROOT_DIR}/docs)
endif()

#-------------------------------------------------------------------------------
# Build C++ tests and install tests to install directory
#-------------------------------------------------------------------------------
if (ENABLE_TESTS)
  add_subdirectory(${SPHERAL_ROOT_DIR}/tests)

  spheral_install_python_tests(${SPHERAL_ROOT_DIR}/tests/ ${SPHERAL_TEST_INSTALL_PREFIX})
  # Always install performance.py in the top of the testing script
  install(FILES ${SPHERAL_ROOT_DIR}/tests/performance.py
    DESTINATION ${SPHERAL_TEST_INSTALL_PREFIX})
endif()

include(${SPHERAL_ROOT_DIR}/cmake/SpheralConfig.cmake)
