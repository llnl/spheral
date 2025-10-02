include_directories(.)

set(DeviceTestLib_sources 
    DeviceTest.cc)

set(DeviceTestLib_headers
    DeviceTest.hh)

spheral_install_python_files(
    )

spheral_add_obj_library(DeviceTestLib SPHERAL_OBJ_LIBS)

