if (NOT ENABLE_CXXONLY)

  # We use the spack generated build time python path to
  # pick up all of the python modules needed for our runtime.
  # Since some modules need to install from directories other than
  # site-packages, we strip away site-packages and pack and
  # copy the full contents of the python library prefix to
  # our virtual env in spheral-setup-venv.sh
  string(REGEX REPLACE "lib\/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}\/site-packages\/?[A-Za-z]*:" ";" VIRTUALENV_PYTHONPATH_COPY "${SPACK_PYTHONPATH}:")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-setup-venv.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/atstest.in"
    "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
  )

  install(FILES 
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
      "${CMAKE_CURRENT_SOURCE_DIR}/spheral_ats.py"
      "${CMAKE_CURRENT_SOURCE_DIR}/devtools/performance_analysis.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/runtime-requirements.txt.in"
    "${CMAKE_CURRENT_BINARY_DIR}/runtime-requirements.txt"
  )

  # Install Runtime python dependencies into the Install directory.
  Spheral_Python_Env(python_runtime_env 
    REQUIREMENTS runtime-requirements.txt 
    PREFIX ${CMAKE_INSTALL_PREFIX}
  )

  if(NOT POLYTOPE_SITE_PACKAGES_PATH)
    install(FILES ${polytope_DIR}/${SPHERAL_SITE_PACKAGES_PATH}/polytope/polytope.so
      DESTINATION ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/polytope/
    )
  endif()

  install(CODE "execute_process( \
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh)"
  )
endif()
