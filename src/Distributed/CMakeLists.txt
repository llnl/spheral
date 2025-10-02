include_directories(.)

set(Distributed_sources
    Communicator.cc
    Process.cc
    RegisterMPIDataTypes.cc
    )

set(Distributed_headers
    allReduce.hh
    Communicator.hh
    Process.hh
    RegisterMPIDataTypes.hh
    RegisterMPIDataTypesInline.hh
    )

if (ENABLE_MPI)
  #----------------------------------------------------------------------------
  # MPI parallel on
  #----------------------------------------------------------------------------
  list(APPEND Distributed_sources waitAllWithDeadlockDetection.cc)

  if (ENABLE_1D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes1d.cc)
  endif()

  if (ENABLE_2D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes2d.cc)
  endif()

  if (ENABLE_3D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes3d.cc)
  endif()

  if (ENABLE_INSTANTIATIONS)
    list(APPEND Distributed_sources
         DistributedBoundaryInst.cc
         NestedGridDistributedBoundaryInst.cc
         TreeDistributedBoundaryInst.cc
         BoundingVolumeDistributedBoundaryInst.cc
         RedistributeNodesInst.cc
         DistributeByXPositionInst.cc
         SortAndDivideRedistributeNodesInst.cc
         SpaceFillingCurveRedistributeNodesInst.cc
         MortonOrderRedistributeNodesInst.cc
         PeanoHilbertOrderRedistributeNodesInst.cc
         VoronoiRedistributeNodesInst.cc
         )
  endif()

  list(APPEND Distributed_headers
       BoundingVolumeDistributedBoundary.hh
       BoundingVolumeDistributedBoundaryInline.hh
       CompareDomainNodesByPosition.hh
       DistributeByXPosition.hh
       DistributedBoundary.hh
       DistributedBoundaryInline.hh
       MortonOrderRedistributeNodes.hh
       NestedGridDistributedBoundary.hh
       NestedGridDistributedBoundaryInline.hh
       NestedGridRedistributeNodes.hh
       NestedGridRedistributeNodesInline.hh
       NestedGridUtilities.hh
       ParmetisRedistributeNodes.hh
       ParmetisRedistributeNodesInline.hh
       PeanoHilbertOrderRedistributeNodes.hh
       RedistributeNodes.hh
       RedistributeNodesInline.hh
       SortAndDivideRedistributeNodes.hh
       SortAndDivideRedistributeNodes1d.hh
       SortAndDivideRedistributeNodes2d.hh
       SortAndDivideRedistributeNodes3d.hh
       SortAndDivideRedistributeNodesInline.hh
       SpaceFillingCurveRedistributeNodes.hh
       TreeDistributedBoundary.hh
       TreeDistributedBoundaryInline.hh
       VoronoiRedistributeNodes.hh
       waitAllWithDeadlockDetection.hh
       )

  list(APPEND Distributed_sources
       DistributedBoundary.cc
       NestedGridDistributedBoundary.cc
       TreeDistributedBoundary.cc
       BoundingVolumeDistributedBoundary.cc
       RedistributeNodes.cc
       DistributeByXPosition.cc
       SortAndDivideRedistributeNodes.cc
       SpaceFillingCurveRedistributeNodes.cc
       MortonOrderRedistributeNodes.cc
       PeanoHilbertOrderRedistributeNodes.cc
       VoronoiRedistributeNodes.cc
       )

  set(MPIPY_FILE_NAME "mpi_mpi4py.py")

else()
  #----------------------------------------------------------------------------
  # MPI parallel off
  #----------------------------------------------------------------------------
  set(MPIPY_FILE_NAME "fakempi.py")

endif()

if (NOT ENABLE_CXXONLY)
  configure_file(${MPIPY_FILE_NAME}
    ${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/Spheral/mpi.py
    COPYONLY)

  install(FILES       ${MPIPY_FILE_NAME}
          DESTINATION ${SPHERAL_SITE_PACKAGES_PATH}/Spheral
          RENAME      mpi.py
          )
endif()

spheral_add_obj_library(Distributed SPHERAL_OBJ_LIBS)
