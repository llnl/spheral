#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Get the set of neighbors for the given node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::vector<size_t>>&
ConnectivityMap<Dimension>::
connectivityForNode(const NodeList<Dimension>* nodeListPtr,
                    const size_t nodeID) const {
  if (mConnectivity.empty()) const_cast<ConnectivityMap<Dimension>*>(this)->buildPerPointConnectivity();
  const bool ghostValid = (mBuildGhostConnectivity or
                           mBuildOverlapConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  CONTRACT_VAR(ghostValid);
  REQUIRE(nodeID >= 0 and 
          ((nodeID < nodeListPtr->numInternalNodes()) or
          (ghostValid and nodeID < nodeListPtr->numNodes())));
  const auto nodeListID = mDataBase.nodeListIndex(nodeListPtr);
  REQUIRE(nodeListID < mDataBase.numNodeLists() and nodeListID < mOffsets.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < mConnectivity.size());
  return mConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Same as above, indicating the NodeList as an integer index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::vector<size_t>>&
ConnectivityMap<Dimension>::
connectivityForNode(const size_t nodeListID,
                    const size_t nodeID) const {
  if (mConnectivity.empty()) const_cast<ConnectivityMap<Dimension>*>(this)->buildPerPointConnectivity();
  const bool ghostValid = (mBuildGhostConnectivity or
                           mBuildOverlapConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  const auto& nodeLists = mDataBase.nodeListPtrs();
  CONTRACT_VAR(ghostValid);
  CONTRACT_VAR(nodeLists);
  REQUIRE(nodeListID < nodeLists.size() and nodeListID < mOffsets.size());
  REQUIRE((nodeID < nodeLists[nodeListID]->numInternalNodes()) or
          (ghostValid and nodeID < nodeLists[nodeListID]->numNodes()));
  REQUIRE(mOffsets[nodeListID] + nodeID < mConnectivity.size());
  return mConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Get the set of neighbors we have some overlap with (common neighbors) for the
// given node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::vector<size_t>>&
ConnectivityMap<Dimension>::
overlapConnectivityForNode(const NodeList<Dimension>* nodeListPtr,
                           const size_t nodeID) const {
  REQUIRE(mBuildOverlapConnectivity);
  const bool ghostValid = (mBuildGhostConnectivity or
                           mBuildOverlapConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  CONTRACT_VAR(ghostValid);
  REQUIRE((nodeID < nodeListPtr->numInternalNodes()) or
          (ghostValid and nodeID < nodeListPtr->numNodes()));
  const auto nodeListID = mDataBase.nodeListIndex(nodeListPtr);
  REQUIRE(nodeListID < mDataBase.numNodeLists() and nodeListID < mOffsets.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < mConnectivity.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < mOverlapConnectivity.size());
  return mOverlapConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Same as above, indicating the NodeList as an integer index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::vector<size_t>>&
ConnectivityMap<Dimension>::
overlapConnectivityForNode(const size_t nodeListID,
                           const size_t nodeID) const {
  REQUIRE(mBuildOverlapConnectivity);
  const bool ghostValid = (mBuildGhostConnectivity or
                           mBuildOverlapConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  const auto& nodeLists = mDataBase.nodeListPtrs();
  CONTRACT_VAR(ghostValid);
  CONTRACT_VAR(nodeLists);
  REQUIRE(nodeListID < nodeLists.size() and nodeListID < mOffsets.size());
  REQUIRE((nodeID < nodeLists[nodeListID]->numInternalNodes()) or
          (ghostValid and nodeID < nodeLists[nodeListID]->numNodes()));
  REQUIRE(mOffsets[nodeListID] + nodeID < mConnectivity.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < mOverlapConnectivity.size());
  return mOverlapConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Compute the number of neighbors for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                    const size_t nodeID) const {
  return mNumNeighbors(mDataBase.nodeListIndex(nodeListPtr), nodeID);
}

template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numNeighborsForNode(const size_t nodeListID,
                    const size_t nodeID) const {
  return mNumNeighbors(nodeListID, nodeID);
}

//------------------------------------------------------------------------------
// Compute the number of overlap neighbors for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numOverlapNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                           const size_t nodeID) const {
  const auto& allneighbors = overlapConnectivityForNode(nodeListPtr, nodeID);
  size_t result = 0u;
  for (const auto& neighbors: allneighbors) result += neighbors.size();
  return result;
}

template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numOverlapNeighborsForNode(const size_t nodeListID,
                           const size_t nodeID) const {
  const auto& nodeLists = mDataBase.nodeListPtrs();
  REQUIRE(nodeListID < nodeLists.size());
  return this->numOverlapNeighborsForNode(nodeLists[nodeListID], nodeID);
}

//------------------------------------------------------------------------------
// A single point to determine if in looping over nodes and neighbors the given
// pair should be calculated or not when we are doing pairs simultaneously.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
calculatePairInteraction(const size_t nodeListi, const size_t i,
                         const size_t nodeListj, const size_t j,
                         const size_t firstGhostNodej) const {
  return ((nodeListj > nodeListi) or
          (nodeListj == nodeListi and j > i) or
          (nodeListj < nodeListi and j >= firstGhostNodej));
}

//------------------------------------------------------------------------------
// Iterators for walking nodes in a prescribed order (domain decomposition 
// independent when needed).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename ConnectivityMap<Dimension>::const_iterator
ConnectivityMap<Dimension>::
begin(const size_t nodeList) const {
  REQUIRE(nodeList < mNodeTraversalIndices.size());
  return mNodeTraversalIndices[nodeList].begin();
}

template<typename Dimension>
inline
typename ConnectivityMap<Dimension>::const_iterator
ConnectivityMap<Dimension>::
end(const size_t nodeList) const {
  REQUIRE(nodeList < mNodeTraversalIndices.size());
  return mNodeTraversalIndices[nodeList].end();
}

//------------------------------------------------------------------------------
// The number of ordered nodes in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numNodes(const size_t nodeList) const {
  REQUIRE(nodeList >= 0 and nodeList < mDataBase.numNodeLists());
  return mNodeTraversalIndices[nodeList].size();
}

//------------------------------------------------------------------------------
// The ith ordered node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
ithNode(const size_t nodeList, const size_t index) const {
  REQUIRE(nodeList < mDataBase.numNodeLists());
  REQUIRE(index < mNodeTraversalIndices[nodeList].size());
  return mNodeTraversalIndices[nodeList][index];
}

//------------------------------------------------------------------------------
// Get the ith NodeList/FluidNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
ConnectivityMap<Dimension>::
nodeList(const size_t index) const {
  const auto& nodeLists = mDataBase.nodeListPtrs();
  REQUIRE(index < nodeLists.size());
  return *(nodeLists[index]);
}

//------------------------------------------------------------------------------
// intersection connectivity, if available
//------------------------------------------------------------------------------
template<typename Dimension>
const std::vector<std::vector<size_t>>&
ConnectivityMap<Dimension>::
intersectionConnectivity(const NodePairIdxType& pair) const {
  const auto itr = mIntersectionConnectivity.find(pair);
  if (itr == mIntersectionConnectivity.end()) VERIFY2(false, "ERROR: attempt to lookup missing intersection connectivity for node pair " << pair);
  return itr->second;
}

}
