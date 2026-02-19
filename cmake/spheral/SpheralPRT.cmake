#----------------------------------------------------------------------------------------
#                                   Spheral_Python_Env
#----------------------------------------------------------------------------------------

# Inputs:
# Python3_EXECUTABLE
# ENV_REQS
# ENV_PREFIX
# PIP_CACHE_DIR
# NETWORK_CONNECTED
# ENV_NAME

set(stamp_file "${ENV_PREFIX}/.venv/${ENV_NAME}_stamp")
set(need_update FALSE)

set(REQUIREMENTS_ARGS)
foreach(_req ${ENV_REQS})
  list(APPEND REQUIREMENTS_ARGS -r)
  list(APPEND REQUIREMENTS_ARGS ${_req})
  if(${_req} IS_NEWER_THAN ${stamp_file})
    set(need_update TRUE)
  endif()
endforeach()

if(NOT need_update)
  message(STATUS "-- Python env ${ENV_NAME} is up to date")
  return()
endif()
message("-- Updating Python env ${ENV_NAME}")

set(PIP_DOWNLOAD_CMD python -m pip download
  --disable-pip-version-check
  --exists-action i
  -d ${PIP_CACHE_DIR})

set(PIP_INSTALL_CMD env MPICC=${MPI_C_COMPILER} MPICXX=${MPI_CXX_COMPILER} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
  python -m pip install
  --disable-pip-version-check
  --no-build-isolation
  --no-index
  --cache-dir ${PIP_CACHE_DIR}
  -f ${PIP_CACHE_DIR})

if(NETWORK_CONNECTED)
  execute_process(
    # Create the virtual env and activate it.
    COMMAND ${Python3_EXECUTABLE} -m venv ${ENV_PREFIX}/.venv &&
    . ${ENV_PREFIX}/.venv/bin/activate &&

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
  )
else()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -m venv ${ENV_PREFIX}/.venv;
    COMMAND . ${ENV_PREFIX}/.venv/bin/activate &&

    ${PIP_INSTALL_CMD} pip==24.1 &&

    ${PIP_INSTALL_CMD} setuptools wheel cython poetry-core &&

    ${PIP_INSTALL_CMD} ${REQUIREMENTS_ARGS}
  )
endif()

# Generate a stamp file inside the virtual environment.
# Existence of this stamp file confirms the environment installed correctly.
# Deletion of the venv dir will then trigger a re-installation.
file(TOUCH "${stamp_file}")
