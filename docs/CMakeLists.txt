if(SPHERAL_ENABLE_DOCS)

  # configured documentation tools and intermediate build results
  set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  get_target_property(DOCS_ACTIVATE_VENV_CMD python_build_env ACTIVATE_VENV)

  add_custom_target(docs ALL
    COMMAND ${DOCS_ACTIVATE_VENV_CMD} && python -m sphinx -q -b html
    -c "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${SPHINX_HTML_DIR}"
    DEPENDS ${DOC_DEPENDS} python_build_env
    COMMENT "Building HTML documentation with Sphinx")
endif()
