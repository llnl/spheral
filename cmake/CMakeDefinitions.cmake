#-----------------------------------------------------------------------------------
# Definitions to be added as compile flags for spheral 
#-----------------------------------------------------------------------------------

option(ENABLE_BOUNDCHECKING "Check bounds on STL types (expensive, Gnu only)" OFF)
option(ENABLE_NAN_EXCEPTIONS "Raise an excpetion when a NAN occurs (Gnu only)" OFF)

# If we're building debug default DBC to All
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("-- building Debug")
  add_definitions("-DDEBUG=1")
  if (NOT DBC_MODE)
    set(DBC_MODE "All")
  endif()
else()
  add_definitions("-DDEBUG=0")
endif()

# The DBC flag
if (DBC_MODE STREQUAL "All")
  message("-- DBC (design by contract) set to All")
  add_definitions("-DDBC_COMPILE_ALL")
elseif (DBC_MODE STREQUAL "Pre")
  message("-- DBC (design by contract) set to Pre")
  add_definitions("-DDBC_COMPILE_PRE")
else()
  message("-- DBC (design by contract) off")
endif()

# Bound checking option -- very expensive at run time
if (ENABLE_BOUNDCHECKING)
  message("-- bound checking enabled")
  add_definitions(-D_GLIBCXX_DEBUG=1)
else()
  message("-- bound checking disabled")
endif()

# NAN handling (Gnu only)
if (ENABLE_NAN_EXCEPTIONS)
  message("-- Enabling NAN floating point exceptions (only applicable to Gnu compilers")
  add_definitions(-DENABLE_NAN_EXCEPTIONS)
endif()

# CXXONLY
if (ENABLE_CXXONLY)
  add_definitions(-DCXXONLY=1)
endif()

# SPHERAL_ENABLE_GLOBALDT_REDUCTION
if (SPHERAL_ENABLE_GLOBALDT_REDUCTION)
  add_definitions(-DGLOBALDT_REDUCTION)
endif()

# SPHERAL_ENABLE_LONGCSDT
if (SPHERAL_ENABLE_LONGCSDT)
  add_definitions(-DLONGCSDT)
endif()

# Default Polytope options
add_definitions(-DUSE_TETGEN=0)
add_definitions(-DUSE_TRIANGLE=0)
add_definitions(-DNOPOLYTOPE=1)

# Are we using Opensubdiv?
if (ENABLE_OPENSUBDIV)
  add_definitions(-DENABLE_OPENSUBDIV)
endif()

# Choose the dimensions we build
if (ENABLE_1D)
  add_definitions(-DSPHERAL1D=1)
endif()
if (ENABLE_2D)
  add_definitions(-DSPHERAL2D=1)
endif()
if (ENABLE_3D)
  add_definitions(-DNOR3D=1)
  add_definitions(-DSPHERAL3D=1)
endif()

if (SPHERAL_ENABLE_TIMERS)
  add_definitions(-DSPHERAL_ENABLE_TIMERS=1)
endif()

if (ENABLE_MPI)
  add_definitions(-DUSE_MPI=1)
endif()

#-------------------------------------------------------------------------------#
# Check if std::span is available
#-------------------------------------------------------------------------------#
include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(span HAS_STD_SPAN_HEADER)

if(HAS_STD_SPAN_HEADER AND CMAKE_CXX_STANDARD GREATER_EQUAL 20)
    message(STATUS "std::span header found")
    add_definitions("-DSPHERAL_USE_STD_SPAN")
else()
    message(STATUS "std::span not available -- falling back to boost::span")
endif()

