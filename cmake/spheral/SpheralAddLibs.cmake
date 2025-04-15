#----------------------------------------------------------------------------------------
#                                   spheral_add_obj_library
#----------------------------------------------------------------------------------------
# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# SPHERAL_BLT_DEPENDS    : REQUIRED : List of external dependencies
# SPHERAL_CXX_DEPENDS    : REQUIRED : List of compiler dependencies
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# SPHERAL_SUBMOD_DEPENDS : REQUIRED : List of submodule dependencies
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# package_name  : REQUIRED : Desired package name
# obj_list_name : REQUIRED : The NAME of the global variable that is the list of
#                            internal target libraries (not the list itself)
# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# Spheral_<package_name> : Target for a given spheral package
# <obj_list_name> : List of internal Spheral target objects, appended with target name
#----------------------------------------------------------------------------------------
function(spheral_add_obj_library package_name obj_list_name)
  # Assumes global variable SPHERAL_BLT_DEPENDS exists and is filled with external dependencies
  get_property(SPHERAL_BLT_DEPENDS GLOBAL PROPERTY SPHERAL_BLT_DEPENDS)
  # Assumes global variable SPHERAL_CXX_DEPENDS exists and is filled with compiler dependencies
  get_property(SPHERAL_CXX_DEPENDS GLOBAL PROPERTY SPHERAL_CXX_DEPENDS)
  # For including files in submodules, currently unused
  get_property(SPHERAL_SUBMOD_INCLUDES GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES)

  if(ENABLE_DEV_BUILD)
    blt_add_library(NAME Spheral_${package_name}
      HEADERS     ${${package_name}_headers}
      SOURCES     ${${package_name}_sources}
      DEPENDS_ON  ${SPHERAL_CXX_DEPENDS} ${SPHERAL_BLT_DEPENDS} 
      SHARED      TRUE)
  else()
    blt_add_library(NAME Spheral_${package_name}
      HEADERS     ${${package_name}_headers}
      SOURCES     ${${package_name}_sources}
      DEPENDS_ON  ${SPHERAL_CXX_DEPENDS} ${SPHERAL_BLT_DEPENDS} 
      OBJECT      TRUE)
  endif()
  target_include_directories(Spheral_${package_name} SYSTEM PUBLIC ${SPHERAL_SUBMOD_INCLUDES})
  # Install the headers
  install(FILES ${${package_name}_headers}
    DESTINATION include/${package_name})
  if(ENABLE_DEV_BUILD)
    # Export target name is either spheral_cxx-targets or spheral_llnlcxx-targets
    if (${obj_list_name} MATCHES "LLNL")
      set(export_target_name spheral_llnlcxx-targets)
    else()
      set(export_target_name spheral_cxx-targets)
    endif()
    install(TARGETS Spheral_${package_name}
      EXPORT ${export_target_name}
      DESTINATION lib)
  endif()
  # Append Spheral_${package_name} to the global object list
  # For example, SPHERAL_OBJ_LIBS or LLNLSPHERAL_OBJ_LIBS
  set_property(GLOBAL APPEND PROPERTY ${obj_list_name} Spheral_${package_name})

endfunction()

#----------------------------------------------------------------------------------------
#                                   spheral_add_cxx_library
#----------------------------------------------------------------------------------------
# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# SPHERAL_BLT_DEPENDS    : REQUIRED : List of external dependencies
# SPHERAL_CXX_DEPENDS    : REQUIRED : List of compiler dependencies
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# SPHERAL_SUBMOD_DEPENDS : REQUIRED : List of submodule dependencies
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# package_name   : REQUIRED : Desired package name (either CXX or LLNLCXX)
# _cxx_obj_list  : REQUIRED : List of internal targets to include
# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# Spheral_<package_name> : Exportable target for interal package name library
#----------------------------------------------------------------------------------------
function(spheral_add_cxx_library package_name _cxx_obj_list)
  # Assumes global variable SPHERAL_BLT_DEPENDS exists and is filled with external dependencies
  get_property(SPHERAL_BLT_DEPENDS GLOBAL PROPERTY SPHERAL_BLT_DEPENDS)
  # Assumes global variable spheral_cxx_depends exists and is filled with compiler dependencies
  get_property(SPHERAL_CXX_DEPENDS GLOBAL PROPERTY SPHERAL_CXX_DEPENDS)
  # For including files in submodules, currently unused
  get_property(SPHERAL_SUBMOD_INCLUDES GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES)
  # Convert package name to lower-case for export target name
  string(TOLOWER ${package_name} lower_case_package)
  set(export_target_name spheral_${lower_case_package}-targets)

  if(ENABLE_DEV_BUILD)
    add_library(Spheral_${package_name} INTERFACE)
    target_link_libraries(Spheral_${package_name} INTERFACE ${_cxx_obj_list})
  else()
    # Build static or shared spheral C++ library
    blt_add_library(NAME Spheral_${package_name}
      HEADERS     ${${package_name}_headers}
      SOURCES     ${${package_name}_sources}
      DEPENDS_ON  ${_cxx_obj_list} ${SPHERAL_CXX_DEPENDS} ${SPHERAL_BLT_DEPENDS}
      SHARED      ${ENABLE_SHARED})
  endif()
  target_include_directories(Spheral_${package_name} SYSTEM PRIVATE ${SPHERAL_SUBMOD_INCLUDES})
  if(ENABLE_CUDA)
    set_target_properties(Spheral_${package_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  ## This cleans up library targets created with object libs. It is turned off as it triggers
  ## a failure on Werror and pedantic builds.
  #set(_properties COMPILE_DEFINITIONS LINK_LIBRARIES LINK_OPTIONS INTERFACE_LINK_OPTIONS COMPILE_OPTIONS INTERFACE_COMPILE_OPTIONS)
  #foreach(_prop ${_properties})
  #  get_target_property(temp_prop Spheral_${package_name} ${_prop})
  #  list(REMOVE_DUPLICATES temp_prop)
  #  set_target_properties(Spheral_${package_name} PROPERTIES ${_prop} "${temp_prop}")
  #endforeach()

  #set_target_properties(Spheral_${package_name} PROPERTIES INTERFACE_LINK_LIBRARIES "")

  # Install Spheral C++ target and set it as an exportable CMake target
  install(TARGETS Spheral_${package_name}
    DESTINATION   lib
    EXPORT        ${export_target_name})
  # Export Spheral target
  install(EXPORT ${export_target_name} DESTINATION lib/cmake)

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(Spheral_${package_name} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endfunction()

#----------------------------------------------------------------------------------------
#                                   spheral_add_pybin11_library_package
#----------------------------------------------------------------------------------------
# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# SPHERAL_BLT_DEPENDS    : REQUIRED : List of external dependencies
# EXTRA_PYB11_SPHERAL_ENV_VARS : OPTIONAL : Additional directories containing python filed, used by LLNLSpheral
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# SPHERAL_SUBMOD_DEPENDS : REQUIRED : List of submodule dependencies
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# package_name     : REQUIRED : Desired package name
# module_list_name : REQUIRED : The NAME of the global variable that is the list of
#                               Spheral python modules (not the list itself)
# INCLUDES : OPTIONAL : Target specific includes
# DEPENDS : OPTIONAL : Target specific dependencies
# SOURCE : OPTIONAL : Target specific sources
# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# Spheral<package_name> : Target for a given Spheral python module
# Spheral<package_name>_src : Target for the PYB11Generated source code for a given Spheral module
# <module_list_name> : List of Spheral python modules, appended with current module name
#----------------------------------------------------------------------------------------

function(spheral_add_pybind11_library package_name module_list_name)

  # Define our arguments
  set(options )
  set(oneValueArgs )
  set(multiValueArgs INCLUDES SOURCES DEPENDS)
  cmake_parse_arguments(${package_name} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  # message("** ${package_name}_INCLUDES: ${${package_name}_INCLUDES}")
  # message("** ${package_name}_SOURCES: ${${package_name}_SOURCES}")
  # message("** ${package_name}_DEPENDS: ${${package_name}_DEPENDS}")

  # List directories in which spheral .py files can be found.
  set(PYTHON_ENV 
      ${EXTRA_PYB11_SPHERAL_ENV_VARS}
      "${SPHERAL_ROOT_DIR}/src/PYB11:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/${PYB11_MODULE_NAME}:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/polytope:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Distributed:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/OpenMP:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/CXXTypes:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Geometry:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/PolyClipper:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Silo:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/DataOutput:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/NodeList:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Field:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldList:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Kernel:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Neighbor:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Material:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FileIO:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/DataBase:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Boundary:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Physics:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Hydro:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ExternalForce:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Gravity:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Integrator:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Utilities:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/NodeGenerators:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldOperations:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SPH:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/RK:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/CRKSPH:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ArtificialViscosity:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SVPH:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Mesh:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Damage:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SolidMaterial:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Strength:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ArtificialConduction:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/KernelIntegrator:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Solvers:"
      "${CMAKE_BINARY_DIR}/src/SimulationControl"
      )

  # Format list into a one line shell friendly format
  STRING(REPLACE ";" "<->" PYTHON_ENV_STR ${PYTHON_ENV})

  string(JOIN ":" PYTHON_ENV_STR ${PYTHON_ENV_STR} ${SPACK_PYTHONPATH})

  # Get the TPL dependencies
  get_property(SPHERAL_BLT_DEPENDS GLOBAL PROPERTY SPHERAL_BLT_DEPENDS)
  list(APPEND SPHERAL_DEPENDS Spheral_CXX ${${package_name}_DEPENDS})

  set(MODULE_NAME Spheral${package_name})
  PYB11Generator_add_module(${package_name}
                            MODULE          ${MODULE_NAME}
                            SOURCE          ${package_name}_PYB11.py
                            DEPENDS         ${SPHERAL_CXX_DEPENDS} ${EXTRA_BLT_DEPENDS} ${SPHERAL_DEPENDS}
                            PYTHONPATH      ${PYTHON_ENV_STR}
                            INCLUDES        ${CMAKE_CURRENT_SOURCE_DIR} ${${package_name}_INCLUDES} ${PYBIND11_ROOT_DIR}/include
                            COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:${SPHERAL_PYB11_TARGET_FLAGS}>"
                            USE_BLT         ON
                            EXTRA_SOURCE    ${${package_name}_SOURCES}
                            INSTALL         OFF
                            VIRTUAL_ENV     python_build_env
                            )
  target_include_directories(${MODULE_NAME} SYSTEM PRIVATE ${SPHERAL_EXTERN_INCLUDES})

  install(TARGETS     ${MODULE_NAME}
          DESTINATION ${SPHERAL_SITE_PACKAGES_PATH}/Spheral
          )
  set_property(GLOBAL APPEND PROPERTY ${module_list_name} ${package_name})
  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(${MODULE_NAME} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

endfunction()
