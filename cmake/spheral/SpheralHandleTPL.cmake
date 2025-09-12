
#----------------------------------------------------------------------------------------
#                                   Spheral_Handle_TPL
#----------------------------------------------------------------------------------------

# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# <lib_name>_DIR        : REQUIRED : The installation location of the TPL
# <lib_name>_INCLUDES   : OPTIONAL : Specific includes for the TPL

# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# <lib_name>     : REQUIRED : The name of the target TPL
# TPL_CMAKE_DIR  : REQUIRED : Directory containing files for each TPL
#                             listing their library names

# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# <lib_name> : Exportable target for the TPL
#----------------------------------------------------------------------------------------

function(Spheral_Handle_TPL lib_name TPL_CMAKE_DIR)

  # Make shortcut variable for directory assigned to ${lib_name}_DIR
  set(lib_dir "${${lib_name}_DIR}")
  if (NOT ${lib_name}_DIR)
    message(FATAL_ERROR "${lib_name}_DIR not set.")
  else()
    message("${lib_name}_DIR set.")
    message("Searching ${lib_name} for : ${${lib_name}_DIR}")
  endif()

  # Find libraries
  set(${lib_name}_libs "")
  # Library names to be set in <tpl>.cmake file
  include(${TPL_CMAKE_DIR}/${lib_name}.cmake)
  # If library names are given, find them
  set(${lib_name}_LIBRARIES )
  foreach(libpath ${${lib_name}_libs})
    find_library(${libpath}_clib NAMES ${libpath}
      PATHS ${lib_dir}/lib
      REQUIRED
      NO_CACHE
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH)
    list(APPEND ${lib_name}_LIBRARIES ${${libpath}_clib})
    message("Importing libraries for ${${libpath}_clib}")
    # find_library treats output as a standard variable from 3.21+
    # We get different behavior on earlier CMake versions.
    if(${CMAKE_VERSION} VERSION_LESS "3.21.0")
      unset(${libpath}_clib CACHE)
    else()
      unset(${libpath}_clib)
    endif()
  endforeach()
  # Find includes by assuming they are explicitly provided as ${lib_name}_INCLUDES
  set(${lib_name}_INCLUDE_DIR ${${lib_name}_INCLUDES})

  # If include directories are not explicity set, look for an include directory
  string(COMPARE EQUAL "${${lib_name}_INCLUDE_DIR}" "" incl_test)
  if(incl_test) # Includes are not explicitly provided
    # Look for an include directory but it isn't required
    if(EXISTS "${lib_dir}/include/")
      set(${lib_name}_INCLUDE_DIR "${lib_dir}/include")
      message("Importing includes for ${lib_name}")
      message("")
    endif()
  else() # Includes are explicitly provided
    # Check to be sure they exist
    if(EXISTS ${${lib_name}_INCLUDE_DIR})
      message("Importing includes for ${lib_name}")
      message("")
    else()
      message(FATAL_ERROR "Include directories for ${lib_name} not found")
    endif()
  endif()

  blt_import_library(NAME ${lib_name}
    TREAT_INCLUDES_AS_SYSTEM ON
    INCLUDES ${${lib_name}_INCLUDE_DIR}
    LIBRARIES ${${lib_name}_LIBRARIES}
    LINK_FLAGS ${${lib_name}_LINK_FLAGS}
    EXPORTABLE ON)
  get_target_property(_is_imported ${lib_name} IMPORTED)
  if(NOT ${_is_imported})
    install(TARGETS ${lib_name}
      EXPORT spheral_cxx-targets
      DESTINATION lib/cmake)
  endif()
  set_target_properties(${lib_name} PROPERTIES EXPORT_NAME spheral::${lib_name})
endfunction()
