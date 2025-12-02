include(CMakeDependentOption)
#-------------------------------------------------------------------------------
# Set Spheral CMake options
#-------------------------------------------------------------------------------
option(SPHERAL_ENABLE_PYTHON "Build Spheral python libraries" ON)

option(SPHERAL_ENABLE_1D "Enable 1D" ON)
option(SPHERAL_ENABLE_2D "Enable 2D" ON)
option(SPHERAL_ENABLE_3D "Enable 3D" ON)
option(SPHERAL_ENABLE_TIMERS "Enable Caliper timers" OFF)
option(SPHERAL_ENABLE_TESTS "Enable tests" ON)
option(SPHERAL_COMBINE_INSTANTIATIONS "Combine instantiations" ON)

option(SPHERAL_ENABLE_ANEOS "Enable the ANEOS equation of state package" ON)
option(SPHERAL_ENABLE_OPENSUBDIV "Enable the Opensubdiv Pixar extension for refining polyhedra" ON)
option(SPHERAL_ENABLE_HELMHOLTZ "Enable the Helmholtz equation of state package" ON)

option(SPHERAL_ENABLE_ARTIFICIAL_CONDUCTION "Enable the artificial conduction package" ON)
option(SPHERAL_ENABLE_EXTERNAL_FORCE "Enable the external force package" ON)
option(SPHERAL_ENABLE_FSISPH "Enable the FSISPH package" ON)
option(SPHERAL_ENABLE_GRAVITY "Enable the gravity package" ON)
option(SPHERAL_ENABLE_GSPH "Enable the GSPH package" ON)
option(SPHERAL_ENABLE_SVPH "Enable the SVPH package" ON)
option(SPHERAL_ENABLE_GLOBALDT_REDUCTION "Enable global allreduce for the time step" ON)
option(SPHERAL_ENABLE_LONGCSDT "Enable longitudinal sound speed time step constraint" ON)
cmake_dependent_option(SPHERAL_ENABLE_SUNDIALS "Enable use of SUNDIALS" ON ENABLE_MPI OFF)
option(SPHERAL_ENABLE_LEOS "Enable use of LEOS" OFF)

option(SPHERAL_NETWORK_CONNECTED "Enable use of network. Disable if using a build cache" ON)
option(SPHERAL_ENABLE_LOGGER "Enable debug log printing" OFF)
option(ENABLE_DEV_BUILD "Build separate internal C++ libraries for faster code development" OFF)
cmake_dependent_option(SPHERAL_UNIFIED_MEMORY "Configure Spheral to use unified memory on GPU machines" OFF SPHERAL_GPU_ENABLED OFF)

# Default is to build shared when python is enabled and build static if python is disabled
set(DEFAULT_STATIC ON)
if(SPHERAL_ENABLE_PYTHON)
  set(DEFAULT_STATIC OFF)
endif()
option(SPHERAL_ENABLE_STATIC "Building static C++ libraries" ${DEFAULT_STATIC})
cmake_dependent_option(SPHERAL_ENABLE_SHARED "Building shared C++ libraries" ON "NOT SPHERAL_ENABLE_STATIC" OFF)

#-------------------------------------------------------------------------------
# Experimental/unused options
#-------------------------------------------------------------------------------

cmake_dependent_option(SPHERAL_ENABLE_RDC "Build using the RDC flag. Currently broken (DO NOT USE)" OFF SPHERAL_GPU_ENABLED OFF)

if (SPHERAL_ENABLE_RDC AND SPHERAL_ENABLE_SHARED)
  message(FATAL_ERROR "Must use SPHERAL_ENABLE_STATIC if enabling SPHERAL_ENABLE_RDC (once it is working)")
endif()

#-------------------------------------------------------------------------------
# Should we build sphinx documentation
#-------------------------------------------------------------------------------
cmake_dependent_option(SPHERAL_ENABLE_DOCS "Enable sphinx Spheral documentation" OFF SPHERAL_ENABLE_PYTHON OFF)

#-------------------------------------------------------------------------------
# For using static TPLs
#-------------------------------------------------------------------------------
option(ENABLE_STATIC_TPL "Assume TPLs are built statically" OFF)

#-------------------------------------------------------------------------------
# Debug options
#-------------------------------------------------------------------------------
option(SPHERAL_ENABLE_BOUNDCHECKING "Check bounds on STL types (expensive, GNU only)" OFF)
option(SPHERAL_ENABLE_NAN_EXCEPTIONS "Raise an excpetion when a NAN occurs (GNU only)" OFF)
set(DEFAULT_DBC_MODE "None")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEFAULT_DBC_MODE "All")
endif()
set(SPHERAL_DBC_MODE ${DEFAULT_DBC_MODE} CACHE STRING "Design-by-contract flag for debugging. Options: 'All', 'Pre', 'None'")
