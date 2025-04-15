#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
// The input iterators must dereference to const NodeList*, and must be 
// const FluidNodeLists underneath.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIterator>
inline
ConnectivityMap<Dimension>::
ConnectivityMap(const NodeListIterator& begin,
                const NodeListIterator& end,
                const bool buildGhostConnectivity,
                const bool buildOverlapConnectivity,
                const bool buildIntersectionConnectivity):
  mNodeLists(),
  mBuildGhostConnectivity(buildGhostConnectivity or buildIntersectionConnectivity),
  mBuildOverlapConnectivity(buildOverlapConnectivity),
  mBuildIntersectionConnectivity(buildIntersectionConnectivity),
  mOffsets(),
  mConnectivity(),
  mNodeTraversalIndices(),
  mKeys(FieldStorageType::CopyFields),
  mCouplingPtr(std::make_shared<NodeCoupling>()) {

  // The private method does the grunt work of filling in the connectivity once we have
  // established the set of NodeLists.
  this->rebuild(begin, end, mBuildGhostConnectivity, mBuildOverlapConnectivity, mBuildIntersectionConnectivity);

  // We'd better be valid after the constructor is finished!
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Rebuild for a given set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIterator>
inline
void
ConnectivityMap<Dimension>::
rebuild(const NodeListIterator& begin,
        const NodeListIterator& end, 
        const bool buildGhostConnectivity,
        const bool buildOverlapConnectivity,
        const bool buildIntersectionConnectivity) {
  mBuildGhostConnectivity = buildGhostConnectivity or buildIntersectionConnectivity;
  mBuildOverlapConnectivity = buildOverlapConnectivity;
  mBuildIntersectionConnectivity = buildIntersectionConnectivity;

  // Copy the set of NodeLists in the order prescribed by the NodeListRegistrar.
  NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
  const bool domainDecompIndependent = registrar.domainDecompositionIndependent();
  const unsigned numNodeLists = std::distance(begin, end);
  mNodeLists.clear();
  mOffsets.resize(numNodeLists);
  std::vector<unsigned> numNodes(numNodeLists);
  for (NodeListIterator itr = begin; itr != end; ++itr) {
    typename std::vector<const NodeList<Dimension>*>::iterator posItr = 
      registrar.findInsertionPoint(*itr, mNodeLists.begin(), mNodeLists.end());
    const unsigned i = std::distance(mNodeLists.begin(), posItr);
    CHECK(i < numNodeLists);
    mNodeLists.insert(posItr, *itr);
    numNodes[i] = ((domainDecompIndependent or mBuildGhostConnectivity or mBuildOverlapConnectivity) ?
                   (*itr)->numNodes() :
                   (*itr)->numInternalNodes());
  }

  // Construct the offsets.
  mOffsets[0] = 0;
  for (unsigned i = 1; i != numNodeLists; ++i) mOffsets[i] = mOffsets[i - 1] + numNodes[i - 1];
  CHECK(mOffsets.size() == mNodeLists.size());

  this->computeConnectivity();
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Are we computing ghost connectivity?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
buildGhostConnectivity() const {
  return mBuildGhostConnectivity;
}

//------------------------------------------------------------------------------
// Are we computing overlap connectivity?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
buildOverlapConnectivity() const {
  return mBuildOverlapConnectivity;
}

//------------------------------------------------------------------------------
// Are we computing intersection connectivity?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
buildIntersectionConnectivity() const {
  return mBuildIntersectionConnectivity;
}

//------------------------------------------------------------------------------
// Get the set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<const NodeList<Dimension>*>&
ConnectivityMap<Dimension>::
nodeLists() const {
  return mNodeLists;
}


//------------------------------------------------------------------------------
// Get the set of NodePairs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodePairList&
ConnectivityMap<Dimension>::
nodePairList() const {
  REQUIRE(mNodePairListPtr);
  return *mNodePairListPtr;
}

template<typename Dimension>
inline
std::shared_ptr<NodePairList>
ConnectivityMap<Dimension>::
nodePairListPtr() const {
  return mNodePairListPtr;
}

//------------------------------------------------------------------------------
// Get the set of neighbors for the given node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<int> >&
ConnectivityMap<Dimension>::
connectivityForNode(const NodeList<Dimension>* nodeListPtr,
                    const int nodeID) const {
  const bool ghostValid = (mBuildGhostConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  CONTRACT_VAR(ghostValid);
  REQUIRE(nodeID >= 0 and 
          ((nodeID < (int)nodeListPtr->numInternalNodes()) or
          (ghostValid and nodeID < (int)nodeListPtr->numNodes())));
  const int nodeListID = std::distance(mNodeLists.begin(),
                                       std::find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
  REQUIRE(nodeListID < (int)mConnectivity.size() and nodeListID < (int)mOffsets.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < (int)mConnectivity.size());
  return mConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Same as above, indicating the NodeList as an integer index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<int> >&
ConnectivityMap<Dimension>::
connectivityForNode(const int nodeListID,
                    const int nodeID) const {
  const bool ghostValid = (mBuildGhostConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  CONTRACT_VAR(ghostValid);
  REQUIRE(nodeListID >= 0 and nodeListID < (int)mConnectivity.size());
  REQUIRE(nodeID >= 0 and 
          ((nodeID < (int)mNodeLists[nodeListID]->numInternalNodes()) or
          (ghostValid and nodeID < (int)mNodeLists[nodeListID]->numNodes())));
  REQUIRE(mOffsets[nodeListID] + nodeID < (int)mConnectivity.size());
  return mConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Get the set of neighbors we have some overlap with (common neighbors) for the
// given node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<int> >&
ConnectivityMap<Dimension>::
overlapConnectivityForNode(const NodeList<Dimension>* nodeListPtr,
                           const int nodeID) const {
  const bool ghostValid = (mBuildGhostConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  CONTRACT_VAR(ghostValid);
  REQUIRE(nodeID >= 0 and 
          ((nodeID < (int)nodeListPtr->numInternalNodes()) or
          (ghostValid and nodeID < (int)nodeListPtr->numNodes())));
  const int nodeListID = std::distance(mNodeLists.begin(),
                                       std::find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
  REQUIRE(nodeListID < (int)mConnectivity.size() and nodeListID < (int)mOffsets.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < (int)mConnectivity.size());
  return mOverlapConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Same as above, indicating the NodeList as an integer index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<int> >&
ConnectivityMap<Dimension>::
overlapConnectivityForNode(const int nodeListID,
                           const int nodeID) const {
  const bool ghostValid = (mBuildGhostConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  CONTRACT_VAR(ghostValid);
  REQUIRE(nodeListID >= 0 and nodeListID < (int)mConnectivity.size());
  REQUIRE(nodeID >= 0 and 
          ((nodeID < (int)mNodeLists[nodeListID]->numInternalNodes()) or
          (ghostValid and nodeID < (int)mNodeLists[nodeListID]->numNodes())));
  REQUIRE(mOffsets[nodeListID] + nodeID < (int)mConnectivity.size());
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
                    const int nodeID) const {
  const std::vector< std::vector<int> >& neighbors = connectivityForNode(nodeListPtr, nodeID);
  size_t result = 0u;
  for (std::vector< std::vector<int> >::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) result += itr->size();
  return result;
}

template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numNeighborsForNode(const int nodeListID,
                    const int nodeID) const {
  REQUIRE(nodeListID < (int)mNodeLists.size());
  return this->numNeighborsForNode(mNodeLists[nodeListID], nodeID);
}

//------------------------------------------------------------------------------
// Compute the number of overlap neighbors for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numOverlapNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                           const int nodeID) const {
  const std::vector< std::vector<int> >& neighbors = overlapConnectivityForNode(nodeListPtr, nodeID);
  size_t result = 0u;
  for (std::vector< std::vector<int> >::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) result += itr->size();
  return result;
}

template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numOverlapNeighborsForNode(const int nodeListID,
                           const int nodeID) const {
  REQUIRE(nodeListID < (int)mNodeLists.size());
  return this->numOverlapNeighborsForNode(mNodeLists[nodeListID], nodeID);
}

//------------------------------------------------------------------------------
// A single point to determine if in looping over nodes and neighbors the given
// pair should be calculated or not when we are doing pairs simultaneously.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
calculatePairInteraction(const int nodeListi, const int i,
                         const int nodeListj, const int j,
                         const int firstGhostNodej) const {
  // const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  // if (domainDecompIndependent) {
  //   return ((nodeListj > nodeListi) or
  //           (nodeListj == nodeListi and 
  //            (mKeys(nodeListj, j) == mKeys(nodeListi, i) ? j > 1 : mKeys(nodeListj, j) > mKeys(nodeListi, i))));
  // } else {
    return ((nodeListj > nodeListi) or
            (nodeListj == nodeListi and j > i) or
            (nodeListj < nodeListi and j >= firstGhostNodej));
  // }
}

//------------------------------------------------------------------------------
// Iterators for walking nodes in a prescribed order (domain decomposition 
// independent when needed).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename ConnectivityMap<Dimension>::const_iterator
ConnectivityMap<Dimension>::
begin(const int nodeList) const {
  REQUIRE(nodeList >= 0 and nodeList < (int)mNodeTraversalIndices.size());
  return mNodeTraversalIndices[nodeList].begin();
}

template<typename Dimension>
inline
typename ConnectivityMap<Dimension>::const_iterator
ConnectivityMap<Dimension>::
end(const int nodeList) const {
  REQUIRE(nodeList >= 0 and nodeList < (int)mNodeTraversalIndices.size());
  return mNodeTraversalIndices[nodeList].end();
}

//------------------------------------------------------------------------------
// The number of ordered nodes in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConnectivityMap<Dimension>::
numNodes(const int nodeList) const {
  REQUIRE(nodeList >= 0 and nodeList < (int)mNodeLists.size());
  return mNodeTraversalIndices[nodeList].size();
}

//------------------------------------------------------------------------------
// The ith ordered node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
ConnectivityMap<Dimension>::
ithNode(const int nodeList, const int index) const {
  REQUIRE(nodeList >= 0 and nodeList < (int)mNodeLists.size());
  REQUIRE(index >= 0 and index < (int)mNodeTraversalIndices[nodeList].size());
  return mNodeTraversalIndices[nodeList][index];
}

//------------------------------------------------------------------------------
// Get the ith NodeList/FluidNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
ConnectivityMap<Dimension>::
nodeList(const int index) const {
  REQUIRE(index >= 0 and index < (int)mNodeLists.size());
  return *(mNodeLists[index]);
}

//------------------------------------------------------------------------------
// The NodeCoupling functor
//------------------------------------------------------------------------------
template<typename Dimension>
NodeCoupling&
ConnectivityMap<Dimension>::
coupling() {
  return *mCouplingPtr;
}

template<typename Dimension>
const NodeCoupling&
ConnectivityMap<Dimension>::
coupling() const {
  return *mCouplingPtr;
}

template<typename Dimension>
void
ConnectivityMap<Dimension>::
coupling(typename ConnectivityMap<Dimension>::NodeCouplingPtr couplingPtr) {
  mCouplingPtr = couplingPtr;
}

//------------------------------------------------------------------------------
// intersection connectivity, if available
//------------------------------------------------------------------------------
template<typename Dimension>
const std::vector<std::vector<int>>&
ConnectivityMap<Dimension>::
intersectionConnectivity(const NodePairIdxType& pair) const {
  const auto itr = mIntersectionConnectivity.find(pair);
  if (itr == mIntersectionConnectivity.end()) VERIFY2(false, "ERROR: attempt to lookup missing intersection connectivity for node pair " << pair);
  return itr->second;
}

}
