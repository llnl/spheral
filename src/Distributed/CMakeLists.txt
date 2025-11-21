include_directories(.)

set(Distributed_sources
    Communicator.cc
    Process.cc
    RegisterMPIDataTypes.cc
    )

set(Distributed_inst )

set(Distributed_headers
    allReduce.hh
    Communicator.hh
    Process.hh
    RegisterMPIDataTypes.hh
    )

if (ENABLE_MPI)
  #----------------------------------------------------------------------------
  # MPI parallel on
  #----------------------------------------------------------------------------
  list(APPEND Distributed_sources waitAllWithDeadlockDetection.cc)

  if (SPHERAL_ENABLE_1D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes1d.cc)
  endif()

  if (SPHERAL_ENABLE_2D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes2d.cc)
  endif()

  if (SPHERAL_ENABLE_3D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes3d.cc)
  endif()

  list(APPEND Distributed_inst
       DistributedBoundary
       NestedGridDistributedBoundary
       TreeDistributedBoundary
       BoundingVolumeDistributedBoundary
       RedistributeNodes
       DistributeByXPosition
       SortAndDivideRedistributeNodes
       SpaceFillingCurveRedistributeNodes
       MortonOrderRedistributeNodes
       PeanoHilbertOrderRedistributeNodes
       VoronoiRedistributeNodes
       )

 list(APPEND Distributed_headers
      BoundingVolumeDistributedBoundary.hh
      CompareDomainNodesByPosition.hh
      DistributeByXPosition.hh
      DistributedBoundary.hh
      DistributedBoundaryInline.hh
      MortonOrderRedistributeNodes.hh
      NestedGridDistributedBoundary.hh
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
      VoronoiRedistributeNodes.hh
      waitAllWithDeadlockDetection.hh
      )

  set(MPIPY_FILE_NAME "mpi_mpi4py.py")

else()
  #----------------------------------------------------------------------------
  # MPI parallel off
  #----------------------------------------------------------------------------
  set(MPIPY_FILE_NAME "fakempi.py")

endif()

if (SPHERAL_ENABLE_PYTHON)
  configure_file(${MPIPY_FILE_NAME}
    ${CMAKE_BINARY_DIR}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/Spheral/mpi.py
    COPYONLY)
endif()
instantiate(Distributed_inst Distributed_sources)

spheral_add_obj_library(Distributed SPHERAL_OBJ_LIBS)
