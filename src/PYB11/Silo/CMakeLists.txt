# Ignore -Wterminate warning that is thrown from DBC.hh VERIFY2

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(SpheralSilo SPHERAL_MODULE_LIST.cc PROPERTIES COMPILE_FLAGS -Wno-terminate)
endif()

spheral_add_pybind11_library(Silo SPHERAL_MODULE_LIST
                             INCLUDES ${SPHERAL_ROOT_DIR}/src/PYB11/Silo SPHERAL_MODULE_LIST)
