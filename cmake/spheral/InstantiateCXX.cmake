#-----------------------------------------------------------------------------------
# Instantiate
#     - Generate and add C++ source files to the project based on dimensional
#       requirements of the build
#
# _inst_var   : *name* of list variable containing base names to instantiate
#               The file ${_inst_var}Inst.cc.py must exist.
#               If instantiation is disabled, ${_inst_var}.cc will be added
#               to the source files, if it exists.
# _source_var : *name* of list variable to append source files to.
# Uses SPHERAL_CURRENT_LIB_TARGET to determine which monolithic CXX target is being created.
#-----------------------------------------------------------------------------------

function(instantiate _inst_var _source_var)
  set(_tmp_source)

  # Create our list of dimension to instantiate
  set(_dims )
  if(SPHERAL_ENABLE_1D)
    list(APPEND _dims 1)
   endif()

  if(SPHERAL_ENABLE_2D)
    list(APPEND _dims 2)
  endif()

  if(SPHERAL_ENABLE_3D)
    list(APPEND _dims 3)
  endif()

  # Iterate over each Instantation file
  foreach(_inst ${${_inst_var}})
    if(SPHERAL_COMBINE_INSTANTIATIONS)
      # Generate a C++ file with the format: <Name>Inst.cc
      set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
      set(_inst_file ${_inst}Inst.cc)

      # Generate the C++ file
      # Uses BLT's python for instantiations to work when building CXX_ONLY as well as with python
      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
        DEPENDS ${_inst_py}
        COMMAND ${Python3_EXECUTABLE} ${SPHERAL_ROOT_DIR}/src/helpers/InstantiationGenerator.py ${_inst_py} ${_inst_file} ${_dims}
        BYPRODUCTS ${_inst_file}
        COMMENT "Generating instantiation ${_inst_file}...")

      # Add the instantiation files to the sources
      list(APPEND _tmp_source ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
    else()
      # Generate a C++ file for each dimension with the format: <Name>Inst<N>d.cc
      foreach(_dim ${_dims})
        set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
        set(_inst_file ${_inst}Inst${_dim}d.cc)

        # Generate the C++ file
        # Uses BLT's python for instantiations to work when building CXX_ONLY as well as with python
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
          DEPENDS ${_inst_py}
          COMMAND ${Python3_EXECUTABLE} ${SPHERAL_ROOT_DIR}/src/helpers/InstantiationGenerator.py ${_inst_py} ${_inst_file} ${_dim}
          BYPRODUCTS ${_inst_file}
          COMMENT "Generating instantiation ${_inst_file}...")

        # Add the instantiation files to the sources
        list(APPEND _tmp_source ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
      endforeach()
    endif()
  endforeach()

  # CMake associates a custom command with targets created in the directory
  # where the command was declared. Keep generated instantiation sources
  # owned by this package directory and make the monolithic library depend
  # on the package-local target below.
  if(NOT ENABLE_DEV_BUILD)
    # Extract package name from input variables
    string(REGEX REPLACE "_.*" "" package_name "${_inst_var}")
    set(target_name Spheral_${package_name}_generated_sources)
    add_custom_target(${target_name}
      DEPENDS ${_tmp_source})
    # Make monolithic CXX target depend on it
    add_dependencies(${SPHERAL_CURRENT_LIB_TARGET} ${target_name})
  endif()

  set(${_source_var} ${${_source_var}} ${_tmp_source} PARENT_SCOPE)

endfunction()
