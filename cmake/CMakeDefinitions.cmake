#-----------------------------------------------------------------------------------
# Definitions to be added as compile flags for spheral 
#-----------------------------------------------------------------------------------

# If we're building debug default DBC to All
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("-- building Debug")
  add_compile_definitions("DEBUG=1")
  if (NOT DBC_MODE)
    set(DBC_MODE "All")
  endif()
else()
  add_compile_definitions("DEBUG=0")
endif()

# The DBC flag
if (DBC_MODE STREQUAL "All")
  message("-- DBC (design by contract) set to All")
  add_compile_definitions("DBC_COMPILE_ALL")
elseif (DBC_MODE STREQUAL "Pre")
  message("-- DBC (design by contract) set to Pre")
  add_compile_definitions("DBC_COMPILE_PRE")
else()
  message("-- DBC (design by contract) off")
endif()

# Bound checking option -- very expensive at run time
if (SPHERAL_ENABLE_BOUNDCHECKING)
  message("-- bound checking enabled")
  add_compile_definitions(_GLIBCXX_DEBUG=1)
else()
  message("-- bound checking disabled")
endif()

set(_comp_flags
  ENABLE_MPI
  ENABLE_HIP
  ENABLE_CUDA
  ENABLE_NAN_EXCEPTIONS
  ENABLE_OPENSUBDIV
  ENABLE_GLOBALDT_REDUCTION
  ENABLE_LONGCSDT
  ENABLE_PYTHON
  ENABLE_TIMERS
  ENABLE_LOGGER)

foreach(_comp ${_comp_flags})
  if(SPHERAL_${_comp})
    add_compile_definitions(SPHERAL_${_comp})
  endif()
endforeach()

# NAN handling (Gnu only)
if (SPHERAL_ENABLE_NAN_EXCEPTIONS)
  message("-- Enabling NAN floating point exceptions (only applicable to GNU compilers")
endif()

# Default Polytope options (currently undefined until polytope is fixed)
#add_definitions(-DUSE_TETGEN)
#add_definitions(-DUSE_TRIANGLE)
#add_definitions(-DUSE_POLYTOPE)

# Choose the dimensions we build
if (SPHERAL_ENABLE_1D)
  add_compile_definitions(SPHERAL1D)
endif()
if (SPHERAL_ENABLE_2D)
  add_compile_definitions(SPHERAL2D)
endif()
if (SPHERAL_ENABLE_3D)
  add_compile_definitions(SPHERAL3D)
endif()
