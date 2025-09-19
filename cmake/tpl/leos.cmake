# If LEOS is built as a debug build, the libyaml name changes
file(GLOB YAML_LIB "${leos_DIR}/lib/libyaml-cpp*")
set(leos_libs libleos.a liblip-cpp.a ${YAML_LIB})
