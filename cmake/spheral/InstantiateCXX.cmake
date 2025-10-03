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
#-----------------------------------------------------------------------------------

function(instantiate _inst_var _source_var)
  set(_tmp_source)

  # Create our list of dimension to instantiate
  set(_dims )

  if(ENABLE_1D)
     list(APPEND _dims 1)
  endif()

  if(ENABLE_2D)
     list(APPEND _dims 2)
  endif()

  if(ENABLE_3D)
     list(APPEND _dims 3)
  endif()

  # Iterate over each Instantation file
  foreach(_inst ${${_inst_var}})
    if(ENABLE_INSTANTIATIONS)
      if(SPHERAL_COMBINE_INSTANTIATIONS)
        # Generate a C++ file with the format: <Name>Inst.cc
        set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
        set(_inst_file ${_inst}Inst.cc)

        string(REPLACE ";" " " _dims_string "${_dims}")

        # Generate the C++ file
        # Uses BLT's python for instantiations to work when building CXX_ONLY as well as with python
        add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
                           DEPENDS ${_inst_py}
                           COMMAND ${Python3_EXECUTABLE} ${SPHERAL_ROOT_DIR}/src/helpers/InstantiationGenerator.py ${_inst_py} ${_inst_file} ${_dims_string}
                           BYPRODUCTS ${_inst_file}
                           COMMENT "Generating instantiation ${_inst_file}..."
                           )

        # Add the instantiation files to the sources
        list(APPEND _tmp_source ${_inst_file})
      else()
        # Generate a C++ file for each dimension with the format: <Name>Inst<N>d.cc
        foreach(_dim ${_dims})
          set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
          set(_inst_file ${_inst}Inst${_dim}d.cc)

          # Generate the C++ file
          # Uses BLT's python for instantiations to work when building CXX_ONLY as well as with python
          add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
                             DEPENDS ${_inst_py}
                             COMMAND ${Python3_EXECUTABLE} ${SPHERAL_ROOT_DIR}/src/helpers/InstantiationGenerator.py ${_inst_py} ${_inst_file} ${_dim}
                             BYPRODUCTS ${_inst_file}
                             COMMENT "Generating instantiation ${_inst_file}..."
                             )

          # Add the instantiation files to the sources
          list(APPEND _tmp_source ${_inst_file})
        endforeach()
      endif()
    else()

      # If the base source file exists, add it to the source files
      # These may not exist for header-only implementations
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}.cc)
        list(APPEND _tmp_source ${_inst}.cc)
      endif()

    endif()

  endforeach()

  set(${_source_var} ${${_source_var}} ${_tmp_source} PARENT_SCOPE)

endfunction()
