#-----------------------------------------------------------------------------------
# Definitions to be added as compile flags for spheral 
#-----------------------------------------------------------------------------------

set(SPHERAL_COMPILE_DEFS )

# If we're building debug default DBC to All
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("-- building Debug")
  list(APPEND SPHERAL_COMPILE_DEFS "DEBUG=1")
else()
  list(APPEND SPHERAL_COMPILE_DEFS "DEBUG=0")
endif()

# The DBC flag
if (SPHERAL_DBC_MODE STREQUAL "All")
  message("-- DBC (design by contract) set to All")
  list(APPEND SPHERAL_COMPILE_DEFS "DBC_COMPILE_ALL")
elseif (SPHERAL_DBC_MODE STREQUAL "Pre")
  message("-- DBC (design by contract) set to Pre")
  list(APPEND SPHERAL_COMPILE_DEFS "DBC_COMPILE_PRE")
else()
  message("-- DBC (design by contract) off")
endif()

# Bound checking option -- very expensive at run time
if (SPHERAL_ENABLE_BOUNDCHECKING)
  message("-- bound checking enabled")
  list(APPEND SPHERAL_COMPILE_DEFS _GLIBCXX_DEBUG=1)
else()
  message("-- bound checking disabled")
endif()

# Add SPHERAL_* definitions based on the corresponding option
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
    list(APPEND SPHERAL_COMPILE_DEFS SPHERAL_${_comp})
  endif()
endforeach()

# NAN handling (Gnu only)
if (SPHERAL_ENABLE_NAN_EXCEPTIONS)
  message("-- Enabling NAN floating point exceptions (only applicable to GNU compilers")
endif()

# Default Polytope options (currently undefined until polytope is fixed)
#list(APPEND SPHERAL_COMPILE_DEFS USE_TETGEN)
#list(APPEND SPHERAL_COMPILE_DEFS USE_TRIANGLE)
#list(APPEND SPHERAL_COMPILE_DEFS USE_POLYTOPE)

# Choose the dimensions we build
if (SPHERAL_ENABLE_1D)
  list(APPEND SPHERAL_COMPILE_DEFS SPHERAL1D)
endif()
if (SPHERAL_ENABLE_2D)
  list(APPEND SPHERAL_COMPILE_DEFS SPHERAL2D)
endif()
if (SPHERAL_ENABLE_3D)
  list(APPEND SPHERAL_COMPILE_DEFS SPHERAL3D)
endif()

#-------------------------------------------------------------------------------#
# Check if std::span is available
#-------------------------------------------------------------------------------#
include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(span HAS_STD_SPAN_HEADER)

if(SPHERAL_UNIFIED_MEMORY)
  message("-- Enabling unified memory for GPU architectures")
  list(APPEND SPHERAL_COMPILE_DEFS SPHERAL_UNIFIED_MEMORY)
  if(HAS_STD_SPAN_HEADER AND CMAKE_CXX_STANDARD GREATER_EQUAL 20)
    message(STATUS "std::span header found")
    list(APPEND SPHERAL_COMPILE_DEFS SPHERAL_USE_STD_SPAN)
  else()
    message(STATUS "std::span not available -- falling back to boost::span")
  endif()
else()
  message("-- Assuming non-unified memory for GPU architectures")
endif()

set_property(GLOBAL PROPERTY SPHERAL_COMPILE_DEFS "${SPHERAL_COMPILE_DEFS}")
