#-----------------------------------------------------------------------------------
# Definitions to be added as compile flags for spheral 
#-----------------------------------------------------------------------------------

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
if (SPHERAL_ENABLE_BOUNDCHECKING)
  message("-- bound checking enabled")
  add_definitions(-D_GLIBCXX_DEBUG=1)
else()
  message("-- bound checking disabled")
endif()

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
  add_definitions(-DSPHERAL1D=1)
endif()
if (SPHERAL_ENABLE_2D)
  add_definitions(-DSPHERAL2D=1)
endif()
if (SPHERAL_ENABLE_3D)
  add_definitions(-DSPHERAL3D=1)
endif()
