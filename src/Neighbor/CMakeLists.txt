include_directories(.)
set(Neighbor_inst
    Neighbor
    PythonNeighbor
    PythonNeighborView
    TreeNeighbor
    ConnectivityMap
   )

set(Neighbor_sources 
    NodePairList.cc)

instantiate(Neighbor_inst Neighbor_sources)

set(Neighbor_headers
    ConnectivityMap.hh
    ConnectivityMapInline.hh
    NeighborBase.hh
    NeighborBaseInline.hh
    Neighbor.hh
    NeighborView.hh
    NeighborViewInline.hh
    NodePairIdxType.hh
    NodePairList.hh
    NodePairListView.hh
    TreeNeighbor.hh
    PythonNeighbor.hh
    PairwiseField.hh
    PairwiseFieldElementAccessor.hh
    PairwiseFieldInline.hh
    PairwiseFieldView.hh
    )

spheral_add_obj_library(Neighbor SPHERAL_OBJ_LIBS)

