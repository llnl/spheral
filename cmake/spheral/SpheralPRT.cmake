#----------------------------------------------------------------------------------------
#                                   Spheral_Python_Runtime_Env
#----------------------------------------------------------------------------------------


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


function(Spheral_Python_Env target_name)

  # Define our arguments
  set(options )
  set(oneValueArgs PREFIX)
  set(multiValueArgs REQUIREMENTS)
  cmake_parse_arguments(${target_name} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(REQUIREMENTS_ARGS)
  foreach(_req ${${target_name}_REQUIREMENTS})
    list(APPEND REQUIREMENTS_ARGS -r)
    list(APPEND REQUIREMENTS_ARGS ${_req})
  endforeach()

  set(PIP_DOWNLOAD_CMD python -m pip download 
                      --disable-pip-version-check 
                      --exists-action i 
                      -d ${SPHERAL_PIP_CACHE_DIR})

  set(PIP_INSTALL_CMD env MPICC=${MPI_C_COMPILER} MPICXX=${MPI_CXX_COMPILER} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} 
                      python -m pip install 
                      --disable-pip-version-check 
                      --no-build-isolation 
                      --no-index 
                      --cache-dir ${SPHERAL_PIP_CACHE_DIR}
                      -f ${SPHERAL_PIP_CACHE_DIR})

  if(SPHERAL_NETWORK_CONNECTED)
    add_custom_command(
      OUTPUT ${${target_name}_PREFIX}/.venv/${target_name}_stamp

      # Create the virtual env and activate it.
      COMMAND ${Python3_EXECUTABLE} -m venv ${${target_name}_PREFIX}/.venv &&
      . ${${target_name}_PREFIX}/.venv/bin/activate &&

      # pip @ 24.1 is the first version that supports local repo paths in requirements
      # files. ATS will fail to install otherwise.
      ${PIP_DOWNLOAD_CMD} pip==24.1 &&
      ${PIP_INSTALL_CMD} pip==24.1 &&

      # Initial packages neede before any of our requirements files.
      ${PIP_DOWNLOAD_CMD} setuptools wheel cython poetry-core &&
      ${PIP_INSTALL_CMD} setuptools wheel cython poetry-core &&

      # Install reuiqrements to virtual env.
      ${PIP_DOWNLOAD_CMD} ${REQUIREMENTS_ARGS} &&
      ${PIP_INSTALL_CMD} ${REQUIREMENTS_ARGS}

      # Generate a stamp file inside the virtual environment.
      # Existence of this stamp file confirms the environment installed correctly.
      # Deletion of the venv dir will then trigger a re-installation.
      COMMAND touch ${${target_name}_PREFIX}/.venv/${target_name}_stamp

      # Changed to the input requirements files will trigger re-installation.
      DEPENDS Python3::Python ${${target_name}_REQUIREMENTS}
    )
  else()
    add_custom_command(
      OUTPUT ${${target_name}_PREFIX}/.venv/${target_name}_stamp
      COMMAND ${Python3_EXECUTABLE} -m venv ${${target_name}_PREFIX}/.venv;
      COMMAND . ${${target_name}_PREFIX}/.venv/bin/activate &&

      ${PIP_INSTALL_CMD} pip==24.1 &&

      ${PIP_INSTALL_CMD} setuptools wheel cython poetry-core &&

      ${PIP_INSTALL_CMD} ${REQUIREMENTS_ARGS}

      COMMAND touch ${${target_name}_PREFIX}/.venv/${target_name}_stamp

      DEPENDS Python3::Python ${${target_name}_REQUIREMENTS}
    )
  endif()

  add_custom_target(${target_name}
    DEPENDS ${${target_name}_PREFIX}/.venv/${target_name}_stamp
  )

  add_custom_target(clean_${target_name}
    COMMAND rm -rf ${${target_name}_PREFIX}/.venv
  )

  set_property(TARGET ${target_name} PROPERTY EXECUTABLE python)
  set_property(TARGET ${target_name} PROPERTY ACTIVATE_VENV . ${${target_name}_PREFIX}/.venv/bin/activate)
endfunction()

