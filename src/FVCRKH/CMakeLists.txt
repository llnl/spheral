set(FVCRKH_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND FVCRKH_sources
       FVCRKHydroBaseInst.cc
       )
endif()

set(FVCRKH_headers
  FVCRKHydroBase.hh
  )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND FVCRKH_headers
     FVCRKHydroBase.cc
     )

set(FVCRKH_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/FVCRKH
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_FVCRKH
                HEADERS     ${FVCRKH_headers}
                SOURCES     ${FVCRKH_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${FVCRKH_includes}
                )

install(TARGETS             Spheral_FVCRKH
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FVCRKH_headers}
        DESTINATION include/FVCRKH
        )

