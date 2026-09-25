include_directories(.)
set(DeviceTestLib_inst
    )


set(DeviceTestLib_sources 
    DeviceTest.cc)

instantiate(DeviceTestLib_inst DeviceTestLib_sources)

set(DeviceTestLib_headers
    DeviceTest.hh)

spheral_install_python_files(
    )

spheral_add_obj_library(DeviceTestLib)

