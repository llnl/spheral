#-----------------------------------------------------------------------------------
# spheral_install_python_files
#     - Copies Python files to the Spheral install path, and byte compiles them
#
# The list of python files should be passed as the arguments
#
# Note, if ENABLE_CXXONLY is set, this function does nothing
#-----------------------------------------------------------------------------------


function(spheral_install_python_files)

  #  file(COPY ${_file} DESTINATION ${CMAKE_BINARY_DIR}/${SPHERAL_SITE_PACKAGES_PATH}/Spheral/)
  foreach(_file ${ARGV})
    get_filename_component(_filename ${_file} NAME)
    configure_file(
      ${_file}
      ${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/Spheral/${_filename}
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    COPYONLY)
  endforeach()

  if (NOT ENABLE_CXXONLY)
    install(FILES ${ARGV}
      DESTINATION ${SPHERAL_SITE_PACKAGES_PATH}/Spheral)
    install(CODE "execute_process( \
    COMMAND ${PYTHON_EXE} -m compileall ${SPHERAL_SITE_PACKAGES_PATH}/Spheral \
            WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX})")
  endif()

endfunction()

#----------------------------------------------------------------------------------------
#                                   spheral_instalL_python_tests
#----------------------------------------------------------------------------------------
# ----------------------
# INPUT VARIABLES
# ----------------------
# test_dir  : REQUIRED : Source directory of tests to install
# test_dest : REQUIRED : Destination for tests
function(spheral_install_python_tests test_dir test_dest)
  install(DIRECTORY ${test_dir}
    USE_SOURCE_PERMISSIONS
    DESTINATION "${test_dest}"
    PATTERN "*CMakeLists.txt*" EXCLUDE
    PATTERN "*.cmake" EXCLUDE
    PATTERN "*.in" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
    PATTERN "performance.py" EXCLUDE
    PATTERN "*~" EXCLUDE)
  # performance.py must be installed in the top test directory
endfunction()
