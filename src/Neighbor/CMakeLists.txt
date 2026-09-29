include_directories(.)
set(Neighbor_inst
    PythonNeighbor
    TreeNeighbor
    ConnectivityMap
   )

set(Neighbor_sources 
    NodePairList.cc)

instantiate(Neighbor_inst Neighbor_sources)

set(Neighbor_headers
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
    NeighborInline.hh
    PythonNeighbor.hh
    TreeNeighbor.hh
    )

spheral_add_obj_library(Neighbor SPHERAL_OBJ_LIBS)

