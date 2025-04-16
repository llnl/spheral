configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Integrator_PYB11.py
  ${CMAKE_CURRENT_BINARY_DIR}/Integrator_PYB11.py)

spheral_add_pybind11_library(Integrator SPHERAL_MODULE_LIST)
