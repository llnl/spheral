#include "Utilities/DBC.hh"

#include <algorithm>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE
inline
NeighborView<Dimension, Descendant>::
NeighborView(NodeListView<Dimension>& nodeList,
             FieldView<Dimension, Vector>& nodeExtents,
             const NeighborSearchType searchType,
             const Scalar kernelExtent):
  NeighborBase<Dimension>(searchType, kernelExtent),
  mNodeListView(nodeList),
  mNodeExtentsView(nodeExtents) {
}

//------------------------------------------------------------------------------
// Force the node extent field to be computed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setNodeExtents() {
  const auto& H = mNodeListView.Hfield();
  for (auto i = 0u; i < mNodeListView.numNodes(); ++i) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

template<typename Dimension, typename Descendant>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setNodeExtents(const Container& nodeIDs) {
  const auto& H = mNodeListView.Hfield();
  for (auto i: nodeIDs) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setInternalNodeExtents() {
  const auto& H = mNodeListView.Hfield();
  const auto n = mNodeListView.numInternalNodes();
  for (auto i = 0u; i < n; ++i) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setGhostNodeExtents() {
  const auto& H = mNodeListView.Hfield();
  const auto ifirst = mNodeListView.firstGhostNode();
  const auto n = mNodeListView.numNodes();
  for (auto i = ifirst; i < n; ++i) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

//------------------------------------------------------------------------------
// Set the master list of nodes based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
template<typename Container>
SPHERAL_HOST_DEVICE 
std::tuple<typename NeighborView<Dimension, Descendant>::size_type,
           typename NeighborView<Dimension, Descendant>::size_type>
NeighborView<Dimension, Descendant>::
setMasterList(const int nodeID,
              Container& masterList,
              Container& coarseNeighbors,
              const bool ghostConnectivity) const {
  REQUIRE(nodeID >= 0 and nodeID < int(mNodeListView.numInternalNodes()));
  const auto& positions = mNodeListView.positions();
  const auto& Hfield = mNodeListView.Hfield();
  return this->setMasterList(positions(nodeID), Hfield(nodeID), masterList, coarseNeighbors, ghostConnectivity);
}

//------------------------------------------------------------------------------
// Set the refined list of potential neighbors based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
template<typename Container>
SPHERAL_HOST_DEVICE 
typename NeighborView<Dimension, Descendant>::size_type
NeighborView<Dimension, Descendant>::
setRefineNeighborList(const int nodeID,
                      const Container& coarseNeighbors,
                      Container& refineNeighbors) const {
  REQUIRE(nodeID >= 0 and nodeID < int(mNodeListView.numInternalNodes()));
  const auto& positions = mNodeListView.positions();
  const auto& Hfield = mNodeListView.Hfield();
  return this->setRefineNeighborList(positions(nodeID), Hfield(nodeID), coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// General precull routine.  Select nodes based on a range of positions and
// extents.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
template<typename Container>
SPHERAL_HOST_DEVICE 
typename NeighborView<Dimension, Descendant>::size_type
NeighborView<Dimension, Descendant>::
precullList(const Vector& minMasterPosition, const Vector& maxMasterPosition,
            const Vector& minMasterExtent, const Vector& maxMasterExtent,
            const typename NeighborView<Dimension, Descendant>::size_type numCoarse,
            Container& coarseList) const {
  REQUIRE(size_t(numCoarse) <= coarseList.size());

  // Grab data we need
  const auto& positions = nodeListView().positions();
  size_type result = 0;

  // What kind of preculling we're doing determines the applied test.
  switch(mSearchType) {

  case NeighborSearchType::GatherScatter:
    for (auto k = 0u; k < numCoarse; ++k) {
      const auto  j = coarseList[k];
      const auto& nodePosition = positions(j);
      const auto  minNodeExtent = nodePosition - mNodeExtentsView(j);
      const auto  maxNodeExtent = nodePosition + mNodeExtentsView(j);
      if (testPointInBox(nodePosition,
                         minMasterExtent,
                         maxMasterExtent) or
          testBoxIntersection(minMasterPosition,
                              maxMasterPosition,
                              minNodeExtent,
                              maxNodeExtent)) {
        coarseList[result] = j;
        ++result;
      }
    }
    break;

  case NeighborSearchType::Gather:
    for (auto k = 0u; k < numCoarse; ++k) {
      const auto  j = coarseList[k];
      const auto& nodePosition = positions(j);
      if (testPointInBox(nodePosition,
                         minMasterExtent,
                         maxMasterExtent)) {
        coarseList[result] = j;
        ++result;
      }
    }
    break;

  case NeighborSearchType::Scatter:
    for (auto k = 0u; k < numCoarse; ++k) {
      const auto  j = coarseList[k];
      const auto& nodePosition = positions(j);
      const auto  minNodeExtent = nodePosition - mNodeExtentsView(j);
      const auto  maxNodeExtent = nodePosition + mNodeExtentsView(j);
      if (testBoxIntersection(minMasterPosition,
                              maxMasterPosition,
                              minNodeExtent,
                              maxNodeExtent)) {
        coarseList[result] = j;
        ++result;
      }
    }
    break;

  default:
    printf("NeighborView::precullList Unhandeled NeighborSearchType\n");
    abort();
  }

  ENSURE(result <= coarseList.size());
  return result;
}

//------------------------------------------------------------------------------
// Static method for setting master/coarse on multiple NodeLists/Neighbors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
template<typename NodeListIterator, typename ContainerOfContainers, typename Container>
SPHERAL_HOST_DEVICE 
typename NeighborView<Dimension, Descendant>::size_type                        // num NodeLists
NeighborView<Dimension, Descendant>::
setMasterNeighborGroup(const typename Dimension::Vector& position,
                       const typename Dimension::SymTensor& H,
                       const NodeListIterator& nodeListBegin,
                       const NodeListIterator& nodeListEnd,
                       const Scalar kernelExtent,
                       ContainerOfContainers& masterLists,
                       ContainerOfContainers& coarseNeighbors,
                       Container& numMasters,
                       Container& numCoarseNeighbors,
                       const bool ghostConnectivity) {

  // Pre-conditions
  const auto maxNumNodeLists = masterLists.size();
  REQUIRE(coarseNeighbors.size() == maxNumNodeLists and
          numMasters.size() == maxNumNodeLists and
          numCoarseNeighbors.size() == maxNumNodeLists);

  // Initializations
  Vector minMasterPosition = Vector::one() * std::numeric_limits<Scalar>::max();
  Vector maxMasterPosition = Vector::one() * std::numeric_limits<Scalar>::lowest();
  Vector minMasterExtent, maxMasterExtent;

  // Make sure the return values are sized appropriately.
  const size_type numNodeLists = std::distance(nodeListBegin, nodeListEnd);
  if (numNodeLists > maxNumNodeLists) return ~numNodeLists;

  // Loop over all the NodeLists.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = **(nodeListBegin + nodeListi);
    const auto& neighbor = nodeList.neighborView();

    // Begin by setting the basic master/coarse neighbor info for each individual NodeList.
    std::tie(numMasters[nodeListi], numCoarseNeighbors[nodeListi]) = neighbor.setMasterList(position, H, masterLists[nodeListi], coarseNeighbors[nodeListi], ghostConnectivity);
    if (numMasters[nodeListi] < 0 or numCoarseNeighbors[nodeListi] < 0) return ~numNodeLists;

    // Accumulate the min/max master node positions and extents.
    const auto& positions = nodeList.positions();
    const auto& extents = neighbor.nodeExtentField();
    for (auto imaster = 0; imaster < numMasters[nodeListi]; ++imaster) {
      const auto  i = masterLists[nodeListi][imaster];
      const auto& ri = positions(i);
      const auto  minExtenti = ri - extents(i);
      const auto  maxExtenti = ri + extents(i);
      for (auto k = 0u; k < Dimension::nDim; ++k) {
        minMasterPosition(k) = std::min(minMasterPosition(k), ri(k));
        maxMasterPosition(k) = std::max(maxMasterPosition(k), ri(k));
        minMasterExtent(k) = std::min(minMasterExtent(k), minExtenti(k));
        maxMasterExtent(k) = std::max(maxMasterExtent(k), maxExtenti(k));
      }
    }
  }

  // Don't forget to compare these positions and extents with the given position/extent. :-P
  const auto extent = HExtent(H, kernelExtent);
  const auto minExtenti = position - extent;
  const auto maxExtenti = position + extent;
  for (auto k = 0u; k < Dimension::nDim; ++k) {
    minMasterPosition(k) = std::min(minMasterPosition(k), position(k));
    maxMasterPosition(k) = std::max(maxMasterPosition(k), position(k));
    minMasterExtent(k) = std::min(minMasterExtent(k), minExtenti(k));
    maxMasterExtent(k) = std::max(maxMasterExtent(k), maxExtenti(k));
  }

  // Loop over the nodes again, and cull the coarse neighbor lists according
  // to the overall min/max master distribution we just calcuated.
  // WARNING!  After this step the neighbor information in the NodeLists is only
  // guaranteed complete for the set of NodeLists passed to this method!
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = **(nodeListBegin + nodeListi);
    const auto& neighbor = nodeList.neighbor();
    numCoarseNeighbors[nodeListi] = neighbor.precullList(minMasterPosition, maxMasterPosition,
                                                         minMasterExtent, maxMasterExtent,
                                                         coarseNeighbors[nodeListi]);
  }
}

}
