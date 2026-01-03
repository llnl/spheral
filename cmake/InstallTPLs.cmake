#-----------------------------------------------------------------------------------
# Define the Third Party Libs to be used here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Set the location of the <tpl>.cmake files
set(TPL_SPHERAL_CMAKE_DIR ${SPHERAL_ROOT_DIR}/cmake/tpl)

# Initialize TPL options
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleExt.cmake)
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralPRT.cmake)

#-----------------------------------------------------------------------------------
# Submodules
#-----------------------------------------------------------------------------------

if (SPHERAL_ENABLE_PYTHON)
  # Find the appropriate Python
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  set(SPHERAL_SITE_PACKAGES_PATH "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages" )
  list(APPEND SPHERAL_CXX_DEPENDS Python3::Python)

  # Set the PYB11Generator path
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${SPHERAL_ROOT_DIR}/extern/PYB11Generator" CACHE PATH "")
  endif()
  # Set the pybind11 path
  if (NOT PYBIND11_ROOT_DIR)
    set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
    set(PYBIND11_NOPYTHON TRUE)
  endif()
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)
  list(APPEND SPHERAL_CXX_DEPENDS pybind11_headers)
  install(TARGETS pybind11_headers
    EXPORT spheral_cxx-targets
    DESTINATION lib/cmake)
  set_target_properties(pybind11_headers PROPERTIES EXPORT_NAME spheral::pybind11_headers)

  # Install Spheral Python Build Dependencies to a python virtual env in the build tree.

  # Need to set up the build env here so the python library targets can depend on
  # python_build_env.
  set(BUILD_REQ_LIST ${SPHERAL_ROOT_DIR}/scripts/build-requirements.txt)
  list(APPEND BUILD_REQ_LIST ${SPHERAL_BINARY_DIR}/scripts/runtime-requirements.txt)
  if(SPHERAL_ENABLE_DOCS)
    list(APPEND BUILD_REQ_LIST ${SPHERAL_ROOT_DIR}/scripts/docs-requirements.txt)
  endif()

  Spheral_Python_Env(python_build_env
    REQUIREMENTS ${BUILD_REQ_LIST}
    PREFIX ${CMAKE_BINARY_DIR}
  )
endif()

# This is currently unfilled in spheral
set_property(GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES "${SPHERAL_SUBMOD_INCLUDES}")

# PolyClipper
if (NOT polyclipper_DIR)
  # If no PolyClipper is specified, build it as an internal target
  set(polyclipper_DIR "${SPHERAL_ROOT_DIR}/extern/PolyClipper")
  # Must set this so PolyClipper doesn't include unnecessary python scripts
  if(NOT SPHERAL_ENABLE_PYTHON)
    set(POLYCLIPPER_ENABLE_PYTHON OFF)
  endif()
  set(POLYCLIPPER_MODULE_GEN OFF)
  set(POLYCLIPPER_ENABLE_DOCS OFF)
  set(POLYCLIPPER_INSTALL_DIR "PolyClipper/include")
  add_subdirectory(${polyclipper_DIR} ${CMAKE_CURRENT_BINARY_DIR}/PolyClipper)
  list(APPEND SPHERAL_BLT_DEPENDS PolyClipperAPI)
  install(TARGETS PolyClipperAPI
    EXPORT spheral_cxx-targets
    DESTINATION lib/cmake)
  set_target_properties(PolyClipperAPI PROPERTIES EXPORT_NAME spheral::PolyClipperAPI)
  message("Found PolyClipper External Package.")
else()
  list(APPEND SPHERAL_EXTERN_LIBS polyclipper)
endif()

#-----------------------------------------------------------------------------------
# Find pre-compiled TPLs
#-----------------------------------------------------------------------------------

# Any targets that used find package must be added to these lists
set(SPHERAL_FP_TPLS )
set(SPHERAL_FP_DIRS )

# Use find_package to get axom (which brings in fmt) and patch fmt
find_package(axom REQUIRED NO_DEFAULT_PATH PATHS ${axom_DIR}/lib/cmake)
list(APPEND SPHERAL_BLT_DEPENDS axom)
list(APPEND SPHERAL_FP_TPLS axom)
list(APPEND SPHERAL_FP_DIRS ${axom_DIR}/lib/cmake)

# This is a hack to handle transitive issues that come
# from using object libraries with newer version of axom
foreach(_comp ${AXOM_COMPONENTS_ENABLED})
  get_target_property(axom_deps axom::${_comp} INTERFACE_LINK_LIBRARIES)
  # strip cuda out so we have control over when cuda is enabled
  list(REMOVE_DUPLICATES axom_deps)
  list(REMOVE_ITEM axom_deps cuda)
  list(APPEND SPHERAL_BLT_DEPENDS ${axom_deps})
endforeach()

message("-----------------------------------------------------------------------------")
# Use find_package to get adiak
find_package(adiak REQUIRED NO_DEFAULT_PATH PATHS ${adiak_DIR}/lib/cmake/adiak)
if(adiak_FOUND)
  list(APPEND SPHERAL_BLT_DEPENDS adiak::adiak)
  list(APPEND SPHERAL_FP_TPLS adiak)
  list(APPEND SPHERAL_FP_DIRS ${adiak_DIR})
  message("Found Adiak External Package.")
endif()

message("-----------------------------------------------------------------------------")
# Use find_package to get polytope
find_package(polytope NO_DEFAULT_PATH PATHS ${polytope_DIR}/lib/cmake)
if(POLYTOPE_FOUND)
  list(APPEND SPHERAL_BLT_DEPENDS polytope)
  list(APPEND SPHERAL_FP_TPLS polytope)
  list(APPEND SPHERAL_FP_DIRS ${polytope_DIR})
  # Install Polytope python library to our site-packages
  if (SPHERAL_ENABLE_PYTHON)
    install(FILES ${POLYTOPE_INSTALL_PREFIX}/${POLYTOPE_SITE_PACKAGES_PATH}/polytope.so
      DESTINATION ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/polytope/
    )
    if (NOT EXISTS ${POLYTOPE_INSTALL_PREFIX}/${POLYTOPE_SITE_PACKAGES_PATH}/polytope.so)
      message(FATAL_ERROR
        "${POLYTOPE_INSTALL_PREFIX}/${POLYTOPE_SITE_PACKAGES_PATH}/polytope.so not found")
    endif()
  endif()
  message("Found Polytope External Package.")
else()
  list(APPEND SPHERAL_EXTERN_LIBS polytope)
endif()

message("-----------------------------------------------------------------------------")
# Use find_package to get caliper
if (SPHERAL_ENABLE_TIMERS)
  # Save caliper_DIR because it gets overwritten by find_package
  if(NOT CONFIG_CALIPER_DIR)
    # Only save if it does not exists already
    set(CONFIG_CALIPER_DIR "${caliper_DIR}" CACHE PATH "Configuration Caliper directory")
  endif()
  find_package(caliper REQUIRED NO_DEFAULT_PATH PATHS ${caliper_DIR}/share/cmake/caliper)
  if(caliper_FOUND)
    list(APPEND SPHERAL_BLT_DEPENDS caliper)
    list(APPEND SPHERAL_FP_TPLS caliper)
    list(APPEND SPHERAL_FP_DIRS ${caliper_DIR})
    message("Found Caliper External Package.")
  endif()
endif()

message("-----------------------------------------------------------------------------")
# HDF5
# This is a hack to allow other codes to use old versions of hdf5
# Ideally, if(NOT ENABLE_STATIC_TPL) would be replaced and the
# find_package call would be moved outside of the if statement:
#
# find_package(hdf5 NO_DEFAULT_PATH PATHS ${hdf5_DIR})
# if (hdf5_FOUND)


if(NOT ENABLE_STATIC_TPL)
  find_package(hdf5 REQUIRED NO_DEFAULT_PATH PATHS ${hdf5_DIR})
  message("Found HDF5 External Package.")
  list(APPEND SPHERAL_FP_TPLS hdf5)
  list(APPEND SPHERAL_FP_DIRS ${hdf5_DIR})
  if(ENABLE_STATIC_TPL)
    list(APPEND SPHERAL_BLT_DEPENDS hdf5-static hdf5_hl-static)
  else()
    list(APPEND SPHERAL_BLT_DEPENDS hdf5-shared hdf5_hl-shared)
  endif()
else()
  list(APPEND SPHERAL_EXTERN_LIBS hdf5)
endif()

message("-----------------------------------------------------------------------------")
find_package(RAJA REQUIRED NO_DEFAULT_PATH PATHS ${raja_DIR})
if (RAJA_FOUND)
  message("Found RAJA External Package.")
endif()

message("-----------------------------------------------------------------------------")
find_package(umpire REQUIRED NO_DEFAULT_PATH PATHS ${umpire_DIR})
if (umpire_FOUND)
  message("Found umpire External Package.")
endif()

message("-----------------------------------------------------------------------------")
# Chai
find_package(chai REQUIRED NO_DEFAULT_PATH PATHS ${chai_DIR})
if(chai_FOUND)
  message("Found chai External Package.")
endif()

list(APPEND SPHERAL_BLT_DEPENDS chai camp RAJA umpire)
list(APPEND SPHERAL_FP_TPLS chai RAJA umpire)
list(APPEND SPHERAL_FP_DIRS ${chai_DIR} ${raja_DIR} ${umpire_DIR})

message("-----------------------------------------------------------------------------")
# Use find_package to get Sundials
if (SPHERAL_ENABLE_SUNDIALS)
  set(SUNDIALS_DIR "${sundials_DIR}")
  find_package(SUNDIALS REQUIRED NO_DEFAULT_PATH
    COMPONENTS kinsol nvecparallel nvecmpiplusx nvecserial
    PATHS ${sundials_DIR}/lib64/cmake/sundials ${sundials_DIR}/lib/cmake/sundials)
  if(SUNDIALS_FOUND)
    set(SUNDIAL_LIBS kinsol nvecparallel nvecmpiplusx nvecserial)
    foreach(_lib ${SUNDIAL_LIBS})
      list(APPEND SPHERAL_BLT_DEPENDS SUNDIALS::${_lib}_static)
    endforeach()
    list(APPEND SPHERAL_FP_TPLS SUNDIALS)
    list(APPEND SPHERAL_FP_DIRS ${sundials_DIR})
    message("Found SUNDIALS External Package.")
  endif()
endif()

set_property(GLOBAL PROPERTY SPHERAL_FP_TPLS ${SPHERAL_FP_TPLS})
set_property(GLOBAL PROPERTY SPHERAL_FP_DIRS ${SPHERAL_FP_DIRS})

message("-----------------------------------------------------------------------------")
if (singularity_eos_DIR)
  find_package(ports-of-call REQUIRED NO_DEFAULT_PATH PATHS ${ports_of_call_DIR})
  if (ports-of-call_FOUND)
    message("Found ports-of-call External Package")
    list(APPEND SPHERAL_BLT_DEPENDS ports-of-call::ports-of-call)
    list(APPEND SPHERAL_FP_TPLS ports-of-call::ports-of-call)
    list(APPEND SPHERAL_FP_DIRS ${ports_of_call_DIR})
    blt_convert_to_system_includes(TARGET ports-of-call::ports-of-call)
  endif()

  find_package(mpark_variant REQUIRED NO_DEFAULT_PATH PATHS ${mpark_variant_DIR})
  if (mpark_variant_FOUND)
    message("Found mpark_variant External Package")
    list(APPEND SPHERAL_BLT_DEPENDS mpark_variant)
    list(APPEND SPHERAL_FP_TPLS mpark_variant)
    list(APPEND SPHERAL_FP_DIRS ${mpark_variant_DIR})
    blt_convert_to_system_includes(TARGET mpark_variant)
  endif()

  find_package(Kokkos REQUIRED NO_DEFAULT_PATH PATHS ${kokkos_DIR})
  if (Kokkos_FOUND)
    message("Found kokkos External Package")
    list(APPEND SPHERAL_FP_DIRS ${ports_of_call_DIR})
    foreach(_tar IN ITEMS Kokkos::kokkos Kokkos::kokkoscore)
      list(APPEND SPHERAL_BLT_DEPENDS ${_tar})
      list(APPEND SPHERAL_FP_TPLS ${_tar})
      blt_convert_to_system_includes(TARGET ${_tar})
    endforeach()
  endif()

  find_package(KokkosKernels REQUIRED NO_DEFAULT_PATH PATHS ${kokkos_kernels_DIR})
  if (KokkosKernels_FOUND)
    message("Found KokkosKernels External Package")
    list(APPEND SPHERAL_BLT_DEPENDS Kokkos::kokkoskernels)
    list(APPEND SPHERAL_FP_TPLS Kokkos::kokkoskernels)
    list(APPEND SPHERAL_FP_DIRS ${kokkos_kernels_DIR})
    blt_convert_to_system_includes(TARGET Kokkos::kokkoskernels)
  endif()

  find_package(spiner REQUIRED NO_DEFAULT_PATH PATHS ${spiner_DIR})
  if (spiner_FOUND)
    message("Found spiner External Package")
    list(APPEND SPHERAL_BLT_DEPENDS spiner::spiner)
    list(APPEND SPHERAL_FP_TPLS spiner::spiner)
    list(APPEND SPHERAL_FP_DIRS ${spiner_DIR})
    blt_convert_to_system_includes(TARGET spiner::spiner)
  endif()

  find_package(singularity-eos REQUIRED NO_DEFAULT_PATH PATHS ${singularity_eos_DIR})
  if(singularity-eos_FOUND)
    message("Found Singularity-EOS External Package")
    list(APPEND SPHERAL_FP_DIRS ${singularity_eos_DIR})
    foreach(_tar IN ITEMS singularity-eos::singularity-eos singularity-eos::eospac-wrapper)
      list(APPEND SPHERAL_BLT_DEPENDS ${_tar})
      list(APPEND SPHERAL_FP_TPLS ${_tar})
      blt_convert_to_system_includes(TARGET ${_tar})
    endforeach()
  endif()
  
  list(APPEND SPHERAL_EXTERN_LIBS eospac)
endif()

message("-----------------------------------------------------------------------------")
# In case we start using find_package on Silo, we should save the silo_DIR path
set(CONFIG_SILO_DIR "${silo_DIR}" CACHE PATH "Configuration Silo directory")
# TPLs that must be imported
list(APPEND SPHERAL_EXTERN_LIBS boost eigen qhull silo)

blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS leos IF SPHERAL_ENABLE_LEOS)
blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS aneos IF SPHERAL_ENABLE_ANEOS)
blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS opensubdiv IF SPHERAL_ENABLE_OPENSUBDIV)

# Create and install target library for each external library
foreach(lib ${SPHERAL_EXTERN_LIBS})
  if(NOT TARGET ${lib})
    Spheral_Handle_TPL(${lib} ${TPL_SPHERAL_CMAKE_DIR})
  endif()
  list(APPEND SPHERAL_BLT_DEPENDS ${lib})
endforeach()
blt_convert_to_system_includes(TARGETS "${SPHERAL_BLT_DEPENDS}")
# Note: SPHERAL_BLT_DEPENDS is made global after this in SetupSpheral.cmake

# This calls LLNLSpheralInstallTPLs.cmake
if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()

if (SPHERAL_ENABLE_PYTHON)
  configure_file(
    ${POLYTOPE_INSTALL_PREFIX}/${SPHERAL_SITE_PACKAGES_PATH}/polytope/polytope.so
    ${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/polytope/polytope.so
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    COPYONLY)

  install(FILES ${POLYTOPE_INSTALL_PREFIX}/${SPHERAL_SITE_PACKAGES_PATH}/polytope/polytope.so
    DESTINATION ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/polytope/
  )
endif()
