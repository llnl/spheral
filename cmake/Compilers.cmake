message("-- C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")

#-------------------------------------------------------------------------------
# Set compiler options
#-------------------------------------------------------------------------------

option(ENABLE_WARNINGS "show compiler warnings" ON)
option(ENABLE_WARNINGS_AS_ERRORS "make warnings errors" OFF)

option(ENABLE_UNUSED_VARIABLE_WARNINGS "show unused variable compiler warnings" ON)
option(ENABLE_UNUSED_PARAMETER_WARNINGS "show unused parameter warnings" OFF)
option(ENABLE_MISSING_INCLUDE_DIR_WARNINGS "Warn for missing include directories" ON)

set(LANG_STR "CXX")
if (ENABLE_HIP)
  set(LANG_STR "HIP")
  # CMake sets the wrong optimization flag for debug mode with HIP
  # But this causes the asan to fail to compile for some reason
  string(REPLACE "-O" "" CMAKE_HIP_FLAGS_DEBUG "${CMAKE_HIP_FLAGS_DEBUG}")
  set(CMAKE_HIP_FLAGS_DEBUG "${CMAKE_HIP_FLAGS_DEBUG} -g -O0")
endif()


set(CXX_COMPILE_FLAGS "")
if (ENABLE_WARNINGS)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    list(APPEND CXX_COMPILE_FLAGS -fdiagnostics-show-option -Wno-unused-command-line-argument -Wno-c++17-extensions)
    if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 18.0.0)
      list(APPEND CXX_COMPILE_FLAGS -Wno-deprecated-declarations -Wno-gnu-zero-variadic-macro-arguments)
      if(CMAKE_CXX_COMPILER_VERSION LESS 20.0.0)
        list(APPEND CXX_COMPILE_FLAGS -Wno-enum-constexpr-conversion)
        # We build some Fortran code from outside sources (like the Helmholtz EOS) that
        # cause building errors if the compiler is too picky...
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-missing-include-dirs")
      endif()
    endif()
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()
message("-- Compiler warnings ${ENABLE_WARNINGS}")

if (ENABLE_WARNINGS_AS_ERRORS)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list(APPEND CXX_COMPILE_FLAGS /W4 /WX)
  else()
    list(APPEND CXX_COMPILE_FLAGS -Wall -Wextra -pedantic -Werror -Wl,--fatal-warnings)
  endif()
  message("-- Treating warnings as errors")
endif()

if (NOT ENABLE_UNUSED_VARIABLE_WARNINGS)
  list(APPEND CXX_COMPILE_FLAGS -Wno-unused-variable)
endif()
message("-- Compiler unused variable warnings ${ENABLE_UNUSED_VARIABLE_WARNINGS}")


if (NOT ENABLE_UNUSED_PARAMETER_WARNINGS)
  list(APPEND CXX_COMPILE_FLAGS -Wno-unused-parameter)
endif()
message("-- Compiler unused parameter warnings ${ENABLE_UNUSED_PARAMETER_WARNINGS}")


if (NOT ENABLE_MISSING_INCLUDE_DIR_WARNINGS)
  list(APPEND CXX_COMPILE_FLAGS -Wno-missing-include-dirs)
endif()
message("-- Compiler missing include dir warnings ${ENABLE_MISSING_INCLUDE_DIR_WARNINGS}")

set(CUDA_WARNING_FLAGS -Xcudafe=\"--diag_suppress=esa_on_defaulted_function_ignored\")

if (SPHERAL_ENABLE_RDC AND ENABLE_HIP)
  list(APPEND CXX_COMPILE_FLAGS -fgpu-rdc)
endif()

if (SPHERAL_ENABLE_ASAN)
  list(APPEND CXX_COMPILE_FLAGS -fsanitize=address)
endif()

set_property(GLOBAL PROPERTY SPHERAL_CXX_FLAGS "${SPHERAL_CXX_FLAGS}"
  "$<$<COMPILE_LANGUAGE:${LANG_STR}>:${CXX_COMPILE_FLAGS}>")
message("-- Using CXX compile flags ${CXX_COMPILE_FLAGS}")

# Currently unused
set_property(GLOBAL PROPERTY SPHERAL_CUDA_FLAGS
  "$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_WARNING_FLAGS}>")

#-------------------------------------------------------------------------------
# Set link options
#-------------------------------------------------------------------------------

set(CXX_LINK_FLAGS "")

if (SPHERAL_ENABLE_RDC AND ENABLE_HIP)
  list(APPEND CXX_LINK_FLAGS -fgpu-rdc)
endif()

if (SPHERAL_ENABLE_ASAN)
  list(APPEND CXX_LINK_FLAGS -fsanitize=address)
  get_property(SPHERAL_ENV_LINES GLOBAL PROPERTY SPHERAL_ENV_LINES)
  message("------------------------Configuring ASAN------------------------------------")
  message("-- Found ASAN libraries at ${ASAN_LIBRARIES}")
  # Modify the hip arch if necessary
  if (ENABLE_HIP)
    list(APPEND SPHERAL_ENV_LINES "export HSA_XNACK=1")
    string(FIND "${CMAKE_HIP_ARCHITECTURES}" "xnack" idx)
    if (idx EQUAL -1)
      set(CMAKE_HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}:xnack+" CACHE STRING "" FORCE)
      message("-- Adding xnack+ to CMAKE_HIP_ARCHITECTURES, new value ${CMAKE_HIP_ARCHITECTURES}")
    endif()
  endif()
  set_property(GLOBAL PROPERTY SPHERAL_ENV_LINES "${SPHERAL_ENV_LINES}")
  message("----------------------------------------------------------------------------")
endif()

# Remove ompstub from fortran link libraries
if("ompstub" IN_LIST CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
  list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ompstub")
endif()

set_property(GLOBAL PROPERTY SPHERAL_LINK_FLAGS "${CXX_LINK_FLAGS}")
message("-- Using link flags ${CXX_LINK_FLAGS}")

#-------------------------------------------------------------------------------
# PYB11 Target Flags
#-------------------------------------------------------------------------------
set(SPHERAL_PYB11_FLAGS ${CXX_COMPILE_FLAGS})
list(APPEND SPHERAL_PYB11_FLAGS
  -O1 # This is necessary, do not remove
  -Wno-unused-local-typedefs
  -Wno-overloaded-virtual)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  list(APPEND SPHERAL_PYB11_FLAGS
    -Wno-self-assign-overloaded
    -Wno-inconsistent-missing-override
    -Wno-delete-non-abstract-non-virtual-dtor
    -Wno-delete-abstract-non-virtual-dtor)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND SPHERAL_PYB11_FLAGS
    -Wno-pedantic
    -fno-var-tracking-assignments)
endif()

set_property(GLOBAL PROPERTY SPHERAL_PYB11_TARGET_FLAGS
  "$<$<COMPILE_LANGUAGE:${LANG_STR}>:${SPHERAL_PYB11_FLAGS}>")

#-------------------------------------------------------------------------------
# Compiler specific flags
#-------------------------------------------------------------------------------
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS -wd11074,11076,654)
  set(SPHERAL_PYB11_TARGET_FLAGS )
endif()
