#----------------------------------------------------------------------------------------
#                                   spheral_initialize_cxx_target
#----------------------------------------------------------------------------------------
# Create the monolithic CXX target
# Priot to calling this, set SPHERAL_CURRENT_LIB_TARGET to either Spheral_CXX or Spheral_LLNLCXX
function(spheral_initialize_cxx_target source_files)
  set(_main_target ${SPHERAL_CURRENT_LIB_TARGET})
  if(ENABLE_DEV_BUILD)
    add_library(${_main_target} INTERFACE)
  else()
    blt_add_library(NAME ${_main_target}
      SOURCES ${source_files}
      DEFINES ${SPHERAL_COMPILE_DEFS}
      DEPENDS_ON ${SPHERAL_CXX_DEPENDS} ${SPHERAL_BLT_DEPENDS}
      SHARED ${SPHERAL_ENABLE_SHARED})
    target_compile_options(${_main_target} PRIVATE ${SPHERAL_CXX_FLAGS})
    target_link_options(${_main_target} PRIVATE ${SPHERAL_LINK_FLAGS})
  endif()
  if(ENABLE_CUDA AND SPHERAL_ENABLE_RDC)
    set_target_properties(${_main_target} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()
endfunction()

# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# SPHERAL_BLT_DEPENDS    : REQUIRED : List of external dependencies
# SPHERAL_CXX_DEPENDS    : REQUIRED : List of compiler dependencies
# SPHERAL_COMPILE_DEFS   : REQUIRED : List of compiler definitions
# SPHERAL_CXX_FLAGS      : REQUIRED : List of C++ compiler options
# SPHERAL_LINK_FLAGS     : REQUIRED : List of link options
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# SPHERAL_CURRENT_LIB_TARGET : REQUIRED : Package target, either Spheral_CXX or Spheral_LLNLCXX
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# package_name  : REQUIRED : Desired package name
# obj_list_name : REQUIRED : The NAME of the global variable that is the list of
#                            internal target libraries in a development build
# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# Spheral_<package_name> : Target for a given spheral package
# <obj_list_name> : List of internal Spheral targets, appended with target name in a
#                   development build
#----------------------------------------------------------------------------------------
function(spheral_add_obj_library package_name)
  # Main package target, either Spheral_CXX or Spheral_LLNLCXX
  set(_main_target ${SPHERAL_CURRENT_LIB_TARGET})
  if(ENABLE_DEV_BUILD)
    blt_add_library(NAME Spheral_${package_name}
      HEADERS     ${${package_name}_headers}
      SOURCES     ${${package_name}_sources}
      DEFINES     ${SPHERAL_COMPILE_DEFS}
      DEPENDS_ON  ${SPHERAL_CXX_DEPENDS} ${SPHERAL_BLT_DEPENDS} 
      SHARED      TRUE)
    target_link_options(Spheral_${package_name} PUBLIC ${SPHERAL_LINK_FLAGS})
    target_link_libraries(${_main_target} INTERFACE Spheral_${package_name})
    target_compile_options(Spheral_${package_name} PRIVATE ${SPHERAL_CXX_FLAGS})
  else()
    target_include_directories(${_main_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    set(package_sources)
    foreach(source IN LISTS ${package_name}_sources)
      get_filename_component(source_path "${source}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
      list(APPEND package_sources "${source_path}")
    endforeach()
    target_sources(${_main_target} PRIVATE ${package_sources})
    # Grab the CXX files
    blt_split_source_list_by_language(
      SOURCES  ${package_sources}
      CXX_LIST package_cxx_sources)
    if(ENABLE_HIP)
      set_property(SOURCE ${package_cxx_sources}
        TARGET_DIRECTORY ${_main_target}
        PROPERTY LANGUAGE HIP)
    endif()
    if(ENABLE_CUDA)
      set_property(SOURCE ${package_cxx_sources}
        TARGET_DIRECTORY ${_main_target}
        PROPERTY LANGUAGE CUDA)
    endif()
  endif()
  # Install the headers
  install(FILES ${${package_name}_headers}
    DESTINATION include/${package_name})
  if(ENABLE_DEV_BUILD)
    # Export target name is either spheral_cxx-targets or spheral_llnlcxx-targets
    if (${_main_target} MATCHES "LLNL")
      set(export_target_name spheral_llnlcxx-targets)
    else()
      set(export_target_name spheral_cxx-targets)
    endif()
    install(TARGETS Spheral_${package_name}
      EXPORT ${export_target_name}
      DESTINATION lib)
  endif()
endfunction()

#----------------------------------------------------------------------------------------
#                                   spheral_install_cxx_library
#----------------------------------------------------------------------------------------
function(spheral_install_cxx_library package_name)
  set(_main_target ${SPHERAL_CURRENT_LIB_TARGET})
  string(TOLOWER ${package_name} lower_case_package)
  set(export_target_name spheral_${lower_case_package}-targets)
  install(TARGETS ${_main_target}
    DESTINATION   lib
    EXPORT        ${export_target_name})

  # Export Spheral target
  install(EXPORT ${export_target_name} DESTINATION lib/cmake)
endfunction()

#----------------------------------------------------------------------------------------
#                                   spheral_add_pybind11_library_package
#----------------------------------------------------------------------------------------
# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# SPHERAL_BLT_DEPENDS    : REQUIRED : List of external dependencies
# EXTRA_PYB11_SPHERAL_ENV_VARS : OPTIONAL : Additional directories containing python filed, used by LLNLSpheral
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# package_name     : REQUIRED : Desired package name
# module_list_name : REQUIRED : The NAME of the global variable that is the list of
#                               Spheral python modules (not the list itself)
# INCLUDES       : OPTIONAL : Target specific includes
# DEPENDS        : OPTIONAL : Target specific dependencies
# SOURCE         : OPTIONAL : Target specific sources
# MULTIPLE_FILES : OPTIONAL : Generate multiple pybind11 output files to parallelize compilation
# IS_SUBMODULE   : OPTIONAL : (default ON) Compile as a submodule of SpheralCompiledPackages
# SUBMODULES     : OPTIONAL : (default "") List of submodules of this module
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
  set(oneValueArgs MULTIPLE_FILES IS_SUBMODULE)
  set(multiValueArgs INCLUDES SOURCES DEPENDS SUBMODULES)
  cmake_parse_arguments(${package_name} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  # message("** ${package_name}_INCLUDES: ${${package_name}_INCLUDES}")
  # message("** ${package_name}_SOURCES: ${${package_name}_SOURCES}")
  # message("** ${package_name}_DEPENDS: ${${package_name}_DEPENDS}")

  if (NOT DEFINED ${package_name}_IS_SUBMODULE)
    set(${package_name}_IS_SUBMODULE "ON")
  endif()
  if (NOT DEFINED ${package_name}_SUBMODULES)
    set(${package_name}_SUBMODULES "")
  endif()

  # List directories in which spheral .py files can be found.
  set(PYTHON_ENV 
      ${EXTRA_PYB11_SPHERAL_ENV_VARS}
      "${SPHERAL_ROOT_DIR}/src/PYB11"
      "${SPHERAL_ROOT_DIR}/src/PYB11/${PYB11_MODULE_NAME}"
      "${SPHERAL_ROOT_DIR}/src/PYB11/polytope"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Distributed"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Threading"
      "${SPHERAL_ROOT_DIR}/src/PYB11/OpenMP"
      "${SPHERAL_ROOT_DIR}/src/PYB11/CXXTypes"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Geometry"
      "${SPHERAL_ROOT_DIR}/src/PYB11/PolyClipper"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Silo"
      "${SPHERAL_ROOT_DIR}/src/PYB11/DataOutput"
      "${SPHERAL_ROOT_DIR}/src/PYB11/NodeList"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldView"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldListView"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Field"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldList"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Kernel"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Neighbor"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Material"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FileIO"
      "${SPHERAL_ROOT_DIR}/src/PYB11/DataBase"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Boundary"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Physics"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Hydro"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ExternalForce"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Gravity"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Integrator"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Utilities"
      "${SPHERAL_ROOT_DIR}/src/PYB11/NodeGenerators"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldOperations"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SPH"
      "${SPHERAL_ROOT_DIR}/src/PYB11/RK"
      "${SPHERAL_ROOT_DIR}/src/PYB11/CRKSPH"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ArtificialViscosity"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SVPH"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Mesh"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Damage"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SolidMaterial"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Strength"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ArtificialConduction"
      "${SPHERAL_ROOT_DIR}/src/PYB11/KernelIntegrator"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Solvers"
      "${CMAKE_BINARY_DIR}/src/SimulationControl"
      )

  # Format python environment lists into a one line shell friendly format
  list(APPEND PYTHON_ENV ${PYTHON_ENV} ${SPACK_PYTHONPATH})
  blt_list_remove_duplicates(TO PYTHON_ENV)
  list(JOIN PYTHON_ENV ":" PYTHON_ENV_STR)

  # Get the TPL dependencies
  get_property(SPHERAL_PYB11_TARGET_FLAGS GLOBAL PROPERTY SPHERAL_PYB11_TARGET_FLAGS)
  list(APPEND SPHERAL_DEPENDS Spheral_CXX ${${package_name}_DEPENDS})

  set(MODULE_NAME Spheral${package_name})
  PYB11Generator_add_module(${package_name}
    MODULE          ${MODULE_NAME}
    SOURCE          ${package_name}_PYB11.py
    DEPENDS         ${SPHERAL_CXX_DEPENDS} ${SPHERAL_BLT_DEPENDS} ${EXTRA_BLT_DEPENDS} ${SPHERAL_DEPENDS}
    DEFINES         ${SPHERAL_COMPILE_DEFS}
    INCLUDES        ${CMAKE_CURRENT_SOURCE_DIR} ${${package_name}_INCLUDES} ${PYBIND11_ROOT_DIR}/include
    COMPILE_OPTIONS ${SPHERAL_PYB11_TARGET_FLAGS}
    USE_BLT         ON
    EXTRA_SOURCE    ${${package_name}_SOURCES}
    INSTALL         OFF # ${SPHERAL_SITE_PACKAGES_PATH}/Spheral
    VIRTUAL_ENV     python_build_env
    MULTIPLE_FILES  ${${package_name}_MULTIPLE_FILES}
    PYTHONPATH      ${PYTHON_ENV_STR}
    IS_SUBMODULE    ${${package_name}_IS_SUBMODULE}
    SUBMODULES      ${${package_name}_SUBMODULES})

  target_include_directories(${MODULE_NAME} SYSTEM PRIVATE ${SPHERAL_EXTERN_INCLUDES})

  add_dependencies(${MODULE_NAME} generate_spheralDimensions)

  if (NOT ${${package_name}_IS_SUBMODULE})
    add_custom_command(TARGET ${MODULE_NAME}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_BINARY_DIR}/lib/${MODULE_NAME}.so
      ${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/Spheral/${MODULE_NAME}.so)
  endif()

  install(TARGETS     ${MODULE_NAME}
          DESTINATION ${SPHERAL_SITE_PACKAGES_PATH}/Spheral)

  set_property(GLOBAL APPEND PROPERTY ${module_list_name} ${package_name})
  get_property(SPHERAL_LINK_FLAGS GLOBAL PROPERTY SPHERAL_LINK_FLAGS)
  target_link_options(Spheral${package_name} PUBLIC ${SPHERAL_LINK_FLAGS})

endfunction()
