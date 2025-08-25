include_directories(.)
set(Neighbor_inst
    ConnectivityMap
    GridCellIndex
    GridCellPlane
    Neighbor
    NestedGridNeighbor
    TreeNeighbor
   )


set(Neighbor_sources 
    NodePairList.cc)

instantiate(Neighbor_inst Neighbor_sources)

set(Neighbor_headers
    ConnectivityMap.hh
    ConnectivityMapInline.hh
    GridCellIndex.hh
    GridCellIndexBase.hh
    GridCellIndexInline.hh
    GridCellPlane.hh
    GridCellPlaneInline.hh
    Neighbor.hh
    NeighborInline.hh
    NestedGridNeighbor.hh
    NestedGridNeighborInline.hh
    NodePairIdxType.hh
    NodePairList.hh
    TreeNeighbor.hh
    PairwiseField.hh
    PairwiseFieldElementAccessor.hh
    PairwiseFieldInline.hh
    PairwiseFieldElementAccessor.hh
    )

spheral_add_obj_library(Neighbor SPHERAL_OBJ_LIBS)

