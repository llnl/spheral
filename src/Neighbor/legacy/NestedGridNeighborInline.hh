#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/intpow2.hh"
#include "GridCellPlane.hh"

#include <numeric>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Determine the appropriate gridlevel for a node by index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
gridLevel(const int nodeID) const {
  REQUIRE(nodeID >=0 and nodeID < (int)this->nodeList().numNodes());
  return gridLevel(this->nodeList().Hfield()(nodeID));
}

//------------------------------------------------------------------------------
// Translate a given grid cell range from one grid level to another.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NestedGridNeighbor<Dimension>::
translateGridCellRange(const GridCellIndex<Dimension>& gridCellMin,
                       const GridCellIndex<Dimension>& gridCellMax,
                       const int gridLevel,
                       const int targetGridLevel,
                       GridCellIndex<Dimension>& targetMin,
                       GridCellIndex<Dimension>& targetMax) const {
  REQUIRE(gridLevel <= mMaxGridLevels);
  REQUIRE(targetGridLevel <= mMaxGridLevels);
  if (targetGridLevel > gridLevel) {
    const int glfactor = intpow2(targetGridLevel - gridLevel);
    targetMin = gridCellMin * glfactor;
    targetMax = (gridCellMax + 1) * glfactor - 1;
  } else {
    const int glfactor = intpow2(gridLevel - targetGridLevel);
    targetMin = gridCellMin / glfactor;
    targetMax = gridCellMax / glfactor;
  }
}

//------------------------------------------------------------------------------
// Determine the appropriate gridlevel for a given H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename OtherHType>
inline
int
NestedGridNeighbor<Dimension>::
gridLevel(const OtherHType& H) const {
  REQUIRE(this->kernelExtent() > 0.0);
  double h = this->HExtent(H, this->kernelExtent()).maxElement();
  return determineGridLevel(h);
}

// Private method that does the work of determining the gridLevel for a given
// smoothing scale.
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
determineGridLevel(const double& h) const {
  REQUIRE(h > 0.0);
  REQUIRE(fuzzyLessThanOrEqual(h, topGridSize()));
  const int result =  max(0, 
                          min(mMaxGridLevels - 1,
                              int(mGridLevelConst0 - log(h)*ln2inverse)));
  ENSURE(result >= 0 and result < mMaxGridLevels);
  return result;
}

//------------------------------------------------------------------------------
// Determine the gridcell that contains a given node on the given gridlevel.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
NestedGridNeighbor<Dimension>::
gridCellIndex(const int nodeID,
              const int gridLevel) const {
  REQUIRE(nodeID >= 0 and nodeID < (int)this->nodeList().numNodes());
  return gridCellIndex(this->nodeList().positions()(nodeID), gridLevel);
}

//------------------------------------------------------------------------------
// Determine the gridcell that a given position corresponds to on a grid level.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// GridCellIndex<Dimension>
// NestedGridNeighbor<Dimension>::
// gridCellIndex(const typename Dimension::Vector& position, int gridLevel) const {
//   REQUIRE(0);
//   return GridCellIndex<Dimension>();
// }

template<>
inline
GridCellIndex<Dim<1> >
NestedGridNeighbor< Dim<1> >::
gridCellIndex(const Dim<1>::Vector& position,
              const int gridLevel) const {
  REQUIRE(gridLevel >= 0 and gridLevel < mMaxGridLevels);
  REQUIRE(mGridCellSizeInv[gridLevel] > 0.0);
  return GridCellIndex<Dim<1> >(int((position.x() - mGridOrigin.x())*mGridCellSizeInv[gridLevel]) - (position.x() < mGridOrigin.x() ? 1 : 0));
}

template<>
inline
GridCellIndex<Dim<2> >
NestedGridNeighbor< Dim<2> >::
gridCellIndex(const Dim<2>::Vector& position,
              const int gridLevel) const {
  REQUIRE(gridLevel >= 0 and gridLevel < mMaxGridLevels);
  REQUIRE(mGridCellSizeInv[gridLevel] > 0.0);
  return GridCellIndex<Dim<2> >(int((position.x() - mGridOrigin.x())*mGridCellSizeInv[gridLevel]) - (position.x() < mGridOrigin.x() ? 1 : 0),
                                int((position.y() - mGridOrigin.y())*mGridCellSizeInv[gridLevel]) - (position.y() < mGridOrigin.y() ? 1 : 0));
}

template<>
inline
GridCellIndex<Dim<3> >
NestedGridNeighbor< Dim<3> >::
gridCellIndex(const Dim<3>::Vector& position,
              const int gridLevel) const {
  REQUIRE(gridLevel >= 0 and gridLevel < mMaxGridLevels);
  REQUIRE(mGridCellSizeInv[gridLevel] > 0.0);
  return GridCellIndex<Dim<3> >(int((position.x() - mGridOrigin.x())*mGridCellSizeInv[gridLevel]) - (position.x() < mGridOrigin.x() ? 1 : 0),
                                int((position.y() - mGridOrigin.y())*mGridCellSizeInv[gridLevel]) - (position.y() < mGridOrigin.y() ? 1 : 0),
                                int((position.z() - mGridOrigin.z())*mGridCellSizeInv[gridLevel]) - (position.z() < mGridOrigin.z() ? 1 : 0));
}

//------------------------------------------------------------------------------
// Get the number of grid levels.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
numGridLevels() const {
  return mMaxGridLevels;
}

// //------------------------------------------------------------------------------
// // The first grid level with daughters.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// int
// NestedGridNeighbor<Dimension>::
// firstParentGridLevel() const {
//   return mFirstParentGridLevel;
// }

// //------------------------------------------------------------------------------
// // Test if the given cell is in the tree.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// bool
// NestedGridNeighbor<Dimension>::
// cellInTree(const GridCellIndex<Dimension>& gridCell,
//            const int gridLevel) const {
//   REQUIRE(gridLevel >= 0 and gridLevel < numGridLevels());
//   REQUIRE(mOccupiedGridCells.size() == numGridLevels());
//   return mDaughterCells[gridLevel].find(gridCell) != mDaughterCells[gridLevel].end();
// }

//------------------------------------------------------------------------------
// Test if the given cell is occupied.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NestedGridNeighbor<Dimension>::
cellOccupied(const GridCellIndex<Dimension>& gridCell,
             const int gridLevel) const {
  REQUIRE(gridLevel >= 0 and gridLevel < numGridLevels());
  REQUIRE((int)mOccupiedGridCells.size() == numGridLevels());
  return headOfGridCell(gridCell, gridLevel) != mEndOfLinkList;
}

// //------------------------------------------------------------------------------
// // Access the internal oct-tree structure.
// // Note this method intercepts calls for grid cells that are not in the map
// // and returns an empty vector for them.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// const std::vector< GridCellIndex<Dimension> >&
// NestedGridNeighbor<Dimension>::
// daughterCells(const GridCellIndex<Dimension>& gridCell,
//               const int gridLevel) const {
//   try {
//     return mDaughterCells[gridLevel][gridCell];
//   } catch (SafeIndexMapKeyError) {
//     return mEmptyNest;
//   }
// }

//------------------------------------------------------------------------------
// Access the occupied grid cells.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector< GridCellIndex<Dimension> > >&
NestedGridNeighbor<Dimension>::
occupiedGridCells() const {
  return mOccupiedGridCells;
}

template<typename Dimension>
inline
const std::vector< GridCellIndex<Dimension> >&
NestedGridNeighbor<Dimension>::
occupiedGridCells(const int gridLevel) const {
  REQUIRE(gridLevel >= 0 and gridLevel < mMaxGridLevels);
  return mOccupiedGridCells[gridLevel];
}

//------------------------------------------------------------------------------
// Get the origin of the coordinate system.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
NestedGridNeighbor<Dimension>::origin() const {
  return mGridOrigin;
}

//------------------------------------------------------------------------------
// The flag indicating the end of a linked list of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
endOfLinkList() const {
  return mEndOfLinkList;
}

//------------------------------------------------------------------------------
// The radius in gridcells a node can influence on it's home grid level.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
gridCellInfluenceRadius() const {
  return mGridCellInfluenceRadius;
}

template<typename Dimension>
inline
void
NestedGridNeighbor<Dimension>::
gridCellInfluenceRadius(const int x) {
  mGridCellInfluenceRadius = x;

  // We have to recalculate the grid level constant.
  mGridLevelConst0 = log(topGridSize()*mGridCellInfluenceRadius)*ln2inverse;
}

//------------------------------------------------------------------------------
// Return the list of inverse grid cell sizes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<double>&
NestedGridNeighbor<Dimension>::gridCellSizeInv() const {
  return mGridCellSizeInv;
}

//------------------------------------------------------------------------------
// Return the list of cells each node occupies.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<GridCellIndex<Dimension> > >&
NestedGridNeighbor<Dimension>::nodeInCell() const {
  return mNodeInCell;
}

//------------------------------------------------------------------------------
// Return the number of occupied grid levels.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
numOccupiedGridLevels() const {
  return std::accumulate(mGridLevelOccupied.begin(), mGridLevelOccupied.end(), 0);
}

//------------------------------------------------------------------------------
// Return the set indicies for the occupied grid levels.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>
NestedGridNeighbor<Dimension>::
occupiedGridLevels() const {
  REQUIRE((int)mGridLevelOccupied.size() == mMaxGridLevels);
  std::vector<int> result;
  result.reserve(mMaxGridLevels);
  for (int gridLevelID = 0; gridLevelID != mMaxGridLevels; ++gridLevelID) {
    if (mGridLevelOccupied[gridLevelID] == 1) result.push_back(gridLevelID);
  }
  return result;
}

//------------------------------------------------------------------------------
// Set the origin of the coordinate system.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NestedGridNeighbor<Dimension>::
origin(const typename Dimension::Vector& /*origin*/) {
  std::cerr << "Warning: we have determined that there is currently a bug with defining the" << std::endl
            << "         origin to be anythign other than zero." << std::endl;
//   mGridOrigin = origin;

//   // If there is enough information, recalculate the node per gridcell info.
//   if (readyToAssignNodes()) updateNodes();
}

//------------------------------------------------------------------------------
// Get the top grid cell size.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NestedGridNeighbor<Dimension>::topGridSize() const {
  if (mMaxGridLevels > 0) {
    return 1.0/(mGridCellSizeInv[0] + FLT_MIN);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the index for the head of cell node for a given grid cell.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
headOfGridCell(const GridCellIndex<Dimension>& gridCellID, int gridLevelID) const {
  REQUIRE(gridLevelID >= 0 and gridLevelID < mMaxGridLevels);
  typename std::map<GridCellIndex<Dimension>, int>::const_iterator
    headOfCellItr = mGridCellHead[gridLevelID].find(gridCellID);
  if (headOfCellItr != mGridCellHead[gridLevelID].end()) {
    CHECK((*headOfCellItr).first == gridCellID);
    return (*headOfCellItr).second;
  } else {
    return mEndOfLinkList;
  }
}

//------------------------------------------------------------------------------
// Return the index for the next node in the linked list for a cell.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NestedGridNeighbor<Dimension>::
nextNodeInCell(int nodeID) const {
  REQUIRE2(nodeID >= 0 and nodeID < (int)this->nodeList().numNodes(),
           "nodeID out of valid range:  " << nodeID << " not in [0, " << this->nodeList().numNodes() << "]");
  REQUIRE2(mNextNodeInCell.size() == this->nodeList().numNodes(),
           "mNodexNodeInCell wrong size:  " << mNextNodeInCell.size() << " != " << this->nodeList().numNodes());
  return mNextNodeInCell[nodeID];
}

//------------------------------------------------------------------------------
// Return the internal nodes in the given cell.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>
NestedGridNeighbor<Dimension>::
internalNodesInCell(const GridCellIndex<Dimension>& gridCell,
                    const int gridLevel) const {
  std::vector<int> result;
  const int firstGhost = this->nodeList().firstGhostNode();
  int nodeID = headOfGridCell(gridCell, gridLevel);
  while (nodeID != mEndOfLinkList) {
    CHECK((nodeID >= 0 and nodeID < (int)this->nodeList().numNodes()) or
          this->nodeList().numNodes() == 0);
    if (nodeID < firstGhost) result.push_back(nodeID);
    nodeID = nextNodeInCell(nodeID);
    CHECK((nodeID >= 0 and nodeID < (int)this->nodeList().numNodes()) or
          this->nodeList().numNodes() == 0 or
          nodeID == mEndOfLinkList);
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the nodes in the given cell.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>
NestedGridNeighbor<Dimension>::
nodesInCell(const GridCellIndex<Dimension>& gridCell,
            const int gridLevel) const {
  std::vector<int> result;
  int nodeID = headOfGridCell(gridCell, gridLevel);
  while (nodeID != mEndOfLinkList) {
    CHECK((nodeID >= 0 and nodeID < (int)this->nodeList().numNodes()) or
          this->nodeList().numNodes() == 0);
    result.push_back(nodeID);
    nodeID = nextNodeInCell(nodeID);
    CHECK((nodeID >= 0 and nodeID < (int)this->nodeList().numNodes()) or
          this->nodeList().numNodes() == 0 or
          nodeID == mEndOfLinkList);
  }
  return result;
}

//------------------------------------------------------------------------------
// Read out the nodes in the given cell and append them to the given vector.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NestedGridNeighbor<Dimension>::
appendNodesInCell(const GridCellIndex<Dimension>& gridCell,
                  const int gridLevel,
                  std::vector<int>& nodes) const {
  int nodeID = headOfGridCell(gridCell, gridLevel);
  while (nodeID != mEndOfLinkList) {
    CHECK((nodeID >= 0 and nodeID < (int)this->nodeList().numNodes()) or
          this->nodeList().numNodes() == 0);
    nodes.push_back(nodeID);
    nodeID = nextNodeInCell(nodeID);
    CHECK((nodeID >= 0 and nodeID < (int)this->nodeList().numNodes()) or
          this->nodeList().numNodes() == 0 or
          nodeID == mEndOfLinkList);
  }
}

//------------------------------------------------------------------------------
// Convert the normal of a plane to an equivalent integer GridCellIndex
// representation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
NestedGridNeighbor<Dimension>::
gridNormal(const typename Dimension::Vector& normal) const {
  GridCellIndex<Dimension> result;
  for (int i = 0; i != Dimension::nDim; ++i) {
    result(i) = (int) (mGridNormalMagnitude*normal(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Link a node into a grid cell linked list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NestedGridNeighbor<Dimension>::
linkNode(int nodeID, int gridLevelID, 
         const GridCellIndex<Dimension>& gridCellID) {

  REQUIRE(nodeID >= 0 and nodeID < (int)this->nodeList().numNodes());
  REQUIRE(gridLevelID >= 0 and gridLevelID < mMaxGridLevels);

  // Make the new node the head of this grid cell.
  CHECK(mNextNodeInCell.size() == this->nodeList().numNodes());
  mNextNodeInCell[nodeID] = headOfGridCell(gridCellID, gridLevelID);
  mGridCellHead[gridLevelID].unsafeIndex(gridCellID) = nodeID;
}

//------------------------------------------------------------------------------
// Unlink a node from a grid cell linked list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NestedGridNeighbor<Dimension>::
unlinkNode(int nodeID, int gridLevelID, 
           const GridCellIndex<Dimension>& gridCellID) {

  REQUIRE(nodeID >= 0 and nodeID < (int)this->nodeList().numNodes());
  REQUIRE(gridLevelID >= 0 and gridLevelID < mMaxGridLevels);

  // Find the node in the grid cell's linked list, and unlink it.
  int linkNodeID = headOfGridCell(gridCellID, gridLevelID);

  // Is this node the head of the grid cell?
  if (linkNodeID == nodeID) {
    if (nextNodeInCell(linkNodeID) == mEndOfLinkList) {
      mGridCellHead[gridLevelID].erase(mGridCellHead[gridLevelID].find(gridCellID));
    } else {
      mGridCellHead[gridLevelID][gridCellID] = nextNodeInCell(linkNodeID);
    }

  } else {
    while (linkNodeID != mEndOfLinkList and 
           nextNodeInCell(linkNodeID) != nodeID) {
      CHECK(linkNodeID >= 0 and linkNodeID < (int)this->nodeList().numNodes());
      linkNodeID = nextNodeInCell(linkNodeID);
    }

    if (linkNodeID != mEndOfLinkList) {
      CHECK(nextNodeInCell(linkNodeID) == nodeID);
      mNextNodeInCell[linkNodeID] = mNextNodeInCell[nodeID];
    }
  }
}

//------------------------------------------------------------------------------
// Produce the refined list of potential neighbors for a single node.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename OtherHType>
inline
void
NestedGridNeighbor<Dimension>::
setNestedRefineNeighborList(const typename Dimension::Vector& position,
                            const OtherHType& H,
                            const std::vector<int>& coarseNeighbors,
                            std::vector<int>& refineNeighbors) const {

  // Bizarrely on realistic test problems we seem to do best by not 
  // culling, but rather simply using the coarse set.
  // Since we have changed the refine neighbor list to just be a pointer at 
  // the coarse list, we're done.

  // const std::vector<int>& coarseList = this->coarseNeighborList();
  // std::vector<int>& refineList = this->accessRefineNeighborList();
  // refineList = coarseList;

  // Determine the maximum extent of this H tensor in each dimension.
  const Vector extent = this->HExtent(H, this->kernelExtent());
  const Vector minExtent = position - extent;
  const Vector maxExtent = position + extent;

  // Use precull to set the refined neighbor list.
  refineNeighbors = this->precullList(position, position, minExtent, maxExtent, coarseNeighbors);

//   // Set the per field refine data caches for this NodeList.
//   this->nodeList().notifyFieldsCacheRefineValues();
}

}
