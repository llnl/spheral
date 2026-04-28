if (SPHERAL_ENABLE_PYTHON)

  # Our runtime-requirements file needs to be configured with the ATS
  # submodule path.
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/runtime-requirements.txt.in"
    "${SPHERAL_BINARY_DIR}/scripts/runtime-requirements.txt"
  )

  # ***********
  # Build Stage
  # ***********

  # Configure the build time scripts to run a build dir spheral environment.
  get_property(SPHERAL_ENV_LINES GLOBAL PROPERTY SPHERAL_ENV_LINES)
  string(REPLACE ";" "\n" SPHERAL_ENV_LINES "${SPHERAL_ENV_LINES}")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in"
    "${CMAKE_BINARY_DIR}/bin/spheral"
    @ONLY
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  # *************
  # Install Stage
  # *************

  # Set up the target for building the install python virtual environment.
  Spheral_Python_Env(python_runtime_env
    REQUIREMENTS ${SPHERAL_BINARY_DIR}/scripts/runtime-requirements.txt
    PREFIX ${CMAKE_INSTALL_PREFIX}
  )

  # Ensure devtool scripts are available on install.
  install(FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/spheral_ats.py"
      "${CMAKE_CURRENT_SOURCE_DIR}/devtools/performance_analysis.py"
    DESTINATION scripts
  )

  install(FILES
    "${CMAKE_BINARY_DIR}/bin/spheral"
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/atstest.sh"
    RENAME "spheral-ats"
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  # Configure scripts for the install stage.
  install( CODE "
  message(\"-- Setting up spheral python virtual environment.\")

  execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target python_runtime_env)
  ")

  # byte compile all of the venv site-packages.
  install(CODE "
  message(\"-- Byte compiling the Spheral Virtual Environment ...\")
  execute_process(
    COMMAND ${CMAKE_INSTALL_PREFIX}/bin/spheral -m compileall -q ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}
  )
  message(\"-- Byte compilation complete.\")
  ")

  # This is unused so no need to do it at the moment
  # message("-- Creating local Silo module path file.")
  # file(WRITE "${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/silo_lib.pth"
  #   "@CONFIG_SILO_DIR@/lib"
  # )
endif()
