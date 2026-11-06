set(Neighbor_inst
  Neighbor
  TreeNeighborView
  TreeNeighbor
  ConnectivityMap
)

if (NOT SPHERAL_ENABLE_GPU)
  list(APPEND Neighbor_inst PythonNeighbor)
endif()

set(Neighbor_sources 
  NodePairList.cc
)

instantiate(Neighbor_inst Neighbor_sources)

set(Neighbor_headers
  NeighborVariant.hh
  NeighborViewVariant.hh
  NeighborSearchType.hh
  HExtent.hh
  NeighborViewContainer.hh
  ConnectivityMap.hh
  ConnectivityMapInline.hh
  NodePairIdxType.hh
  NodePairList.hh
  NodePairListView.hh
  PairwiseField.hh
  PairwiseFieldElementAccessor.hh
  PairwiseFieldInline.hh
  PairwiseFieldView.hh
  Neighbor.hh
  NeighborView.hh
  NeighborViewInline.hh
  PythonNeighbor.hh
  PythonNeighborCallView.hh
  TreeNeighborView.hh
  TreeNeighborViewInline.hh
  TreeNeighbor.hh
)

spheral_add_package(Neighbor)
