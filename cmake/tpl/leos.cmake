set(leos_libs libleos_C.a libleos.a liblip-cpp.a libyaml-cpp.a libleospact.a)
# If we ever support debug TPL builds we'll need the following since the name of the libarary changes
# if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#   list(APPEND leos_libs libyaml-cppd.a)
# else()
#   list(APPEND leos_libs libyaml-cpp.a)
# endif()
