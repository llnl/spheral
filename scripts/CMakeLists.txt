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
  set(SPHERAL_EXEC_DIR ${CMAKE_BINARY_DIR})
  set(SPHERAL_SCRIPT_DIR ${SPHERAL_ROOT_DIR})
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in"
    "${CMAKE_BINARY_DIR}/bin/spheral"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/atstest.in"
    "${CMAKE_BINARY_DIR}/bin/spheral-ats"
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

  # Move the contents of our site-packages dir into the venv.
  install(DIRECTORY ${CMAKE_INSTALL_PREFIX}/${SPHERAL_SITE_PACKAGES_PATH}/
          DESTINATION ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH})

  # Ensure devtool scripts are available on install.
  install(FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/spheral_ats.py"
      "${CMAKE_CURRENT_SOURCE_DIR}/devtools/performance_analysis.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
  )

  # Configure scripts for the install stage.
  install( CODE "
  message(\"-- Setting up spheral python virtual environment.\")

  execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target python_runtime_env)

  message(\"-- Setting up spheral executable scripts.\")

  set(SPHERAL_EXEC_DIR ${CMAKE_INSTALL_PREFIX})
  set(SPHERAL_SCRIPT_DIR ${CMAKE_INSTALL_PREFIX})

  configure_file(
    \"${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in\"
    \"${CMAKE_INSTALL_PREFIX}/bin/spheral\"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  configure_file(
    \"${CMAKE_CURRENT_SOURCE_DIR}/atstest.in\"
    \"${CMAKE_INSTALL_PREFIX}/bin/spheral-ats\"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )
  ")

  # byte compile all of the venv site-packages.
  install(CODE "
  message(\"-- Byte compiling the Spheral Virtual Envionment ...\")
  execute_process(
    COMMAND ${CMAKE_INSTALL_PREFIX}/bin/spheral -m compileall -q ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}
  )
  message(\"-- Byte compilation complete.\")
  ")

  message("-- Creating local Silo module path file.")
  file(WRITE "${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/silo_lib.pth"
    "@CONFIG_SILO_DIR@/lib"
  )
endif()
