set(Neighbor_inst
  Neighbor
  PythonNeighbor
  TreeNeighborView
  TreeNeighbor
  ConnectivityMap
)

set(Neighbor_sources 
  NodePairList.cc
)

instantiate(Neighbor_inst Neighbor_sources)

set(Neighbor_headers
  NeighborVariants.hh
  NeighborSearchType.hh
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
  TreeNeighborView.hh
  TreeNeighbor.hh
)

spheral_add_package(Neighbor)

