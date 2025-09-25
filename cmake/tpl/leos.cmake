# If LEOS is built as a debug build, the libyaml name changes
file(GLOB FULL_YAML_LIB "${leos_DIR}/lib/libyaml-cpp*")
get_filename_component(YAML_LIB "${FULL_YAML_LIB}" NAME)
set(leos_libs libleos.a liblip-cpp.a ${YAML_LIB})
