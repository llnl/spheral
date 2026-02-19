
if (SPHERAL_ENABLE_PYTHON)

  # Our runtime-requirements file needs to be configured with the ATS
  # submodule path.
  configure_file(
    "${SPHERAL_SCRIPT_DIR}/runtime-requirements.txt.in"
    "${SPHERAL_BINARY_DIR}/scripts/runtime-requirements.txt"
  )

  set(SPHERAL_PIP_CACHE_DIR ~/.cache/spheral_pip)
  if (DEFINED ENV{SYS_TYPE})
    set(SPHERAL_ENV_SYS_TYPE $ENV{SYS_TYPE})
    set(SPHERAL_PIP_CACHE_DIR ${SPHERAL_PIP_CACHE_DIR}/$ENV{SYS_TYPE})
  endif()
  set(SPHERAL_PIP_CACHE_DIR ${SPHERAL_PIP_CACHE_DIR}/${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION})

  # Assume we have network connectivity.
  if(NOT DEFINED SPHERAL_NETWORK_CONNECTED)
    set(SPHERAL_NETWORK_CONNECTED True)
  endif()

  add_custom_target(clean_pip_cache
    COMMAND rm -rf ${SPHERAL_PIP_CACHE_DIR}
  )

  set(VENV_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/spheral/SpheralPRT.cmake")

  # ***********
  # Build Stage
  # ***********
  # Install Spheral Python Build Dependencies to a python virtual env in the build tree.
  # Need to set up the build env here so the python library targets can depend on
  # python_build_env.
  set(ENV_NAME "python_build_env")
  set(BUILD_REQ_LIST ${SPHERAL_SCRIPT_DIR}/build-requirements.txt)
  list(APPEND BUILD_REQ_LIST ${SPHERAL_BINARY_DIR}/scripts/runtime-requirements.txt)
  if(SPHERAL_ENABLE_DOCS)
    list(APPEND BUILD_REQ_LIST ${SPHERAL_SCRIPT_DIR}/docs-requirements.txt)
  endif()

  list(JOIN BUILD_REQ_LIST ";" BUILD_REQ_LIST_STR)
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/.venv/${ENV_NAME}_stamp"

    # Run the SpheralPRT.cmake script
    COMMAND "${CMAKE_COMMAND}"
    "-DPython3_EXECUTABLE=${Python3_EXECUTABLE}"
    "-DENV_REQUIREMENTS=${BUILD_REQ_LIST_STR}"
    "-DENV_PREFIX=${CMAKE_BINARY_DIR}"
    "-DPIP_CACHE_DIR=${SPHERAL_PIP_CACHE_DIR}"
    "-DNETWORK_CONNECTED=${SPHERAL_NETWORK_CONNECTED}"
    "-DENV_NAME=${ENV_NAME}"
    -P ${VENV_SCRIPT}
    # Changed to the input requirements files will trigger re-installation.
    DEPENDS Python3::Python ${BUILD_REQ_LIST}
  )

  add_custom_target(${ENV_NAME}
    DEPENDS ${CMAKE_BINARY_DIR}/.venv/${ENV_NAME}_stamp
  )

  add_custom_target(clean_${ENV_NAME}
    COMMAND rm -rf ${CMAKE_BINARY_DIR}/.venv
  )

  set_property(TARGET ${ENV_NAME} PROPERTY EXECUTABLE python)
  set_property(TARGET ${ENV_NAME} PROPERTY ACTIVATE_VENV . ${CMAKE_BINARY_DIR}/.venv/bin/activate)

  # Copy the spheral-env to run a build time environment.
  file(COPY "${SPHERAL_SCRIPT_DIR}/spheral-env.sh"
    DESTINATION "${CMAKE_BINARY_DIR}/bin/spheral"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  # *************
  # Install Stage
  # *************

  # Set up the target for building the install python virtual environment.
  set(ENV_NAME "python_runtime_env")
  set(RUNTIME_REQ_LIST ${SPHERAL_BINARY_DIR}/scripts/runtime-requirements.txt)
  list(JOIN RUNTIME_REQ_LIST ";" RUNTIME_REQ_LIST_STR)

  install(CODE "
  set(env_prefix \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX})
  execute_process(
    COMMAND ${CMAKE_COMMAND}
    -DPython3_EXECUTABLE=${Python3_EXECUTABLE}
    \"-DENV_REQUIREMENTS=${RUNTIME_REQ_LIST_STR}\"
    -DENV_PREFIX=\${env_prefix}
    -DPIP_CACHE_DIR=${SPHERAL_PIP_CACHE_DIR}
    -DNETWORK_CONNECTED=${SPHERAL_NETWORK_CONNECTED}
    -DENV_NAME=${ENV_NAME}
    -P${VENV_SCRIPT}
  )
  ")

  # Ensure devtool scripts are available on install.
  install(FILES
      "${SPHERAL_SCRIPT_DIR}/spheral_ats.py"
      "${SPHERAL_SCRIPT_DIR}/devtools/performance_analysis.py"
    DESTINATION "scripts"
  )

  install(FILES
    "${SPHERAL_SCRIPT_DIR}/spheral-env.sh"
    RENAME "spheral"
    DESTINATION "bin"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  install(FILES
    "${SPHERAL_SCRIPT_DIR}/atstest.sh"
    RENAME "spheral-ats"
    DESTINATION "bin"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  )

  # byte compile all of the venv site-packages.
  install(CODE "
  message(\"-- Byte compiling the Spheral Virtual Envionment ...\")
  set(env_prefix \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX})
  execute_process(
    COMMAND \${env_prefix}/bin/spheral -m compileall -q \${env_prefix}/${SPHERAL_SITE_PACKAGES_PATH}
  )
  message(\"-- Byte compilation complete.\")
  ")

  # This is unused so no need to do it at the moment
  # message("-- Creating local Silo module path file.")
  # file(WRITE "${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/silo_lib.pth"
  #   "@CONFIG_SILO_DIR@/lib"
  # )
endif()
