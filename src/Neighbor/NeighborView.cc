//---------------------------------Spheral++----------------------------------//
// NeighborView -- An abstract base class for the Neighbor objects.
//----------------------------------------------------------------------------//
#include "Neighbor/NeighborView.hh"
#include "Utilities/testBoxIntersection.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
NeighborView<Dimension>::
NeighborView(NodeListView<Dimension>& nodeList,
             FieldView<Dimension, Vector>& nodeExtents,
             const NeighborSearchType searchType,
             const double kernelExtent):
  NeighborBase<Dimension>(searchType, kernelExtent),
  mNodeListView(nodeList),
  mNodeExtentsView(nodeExtents) {
}

//------------------------------------------------------------------------------
// Force the node extent field to be computed.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
void
NeighborView<Dimension>::
setNodeExtents() {
  const auto& H = mNodeListView.Hfield();
  for (auto i = 0u; i < mNodeListView.numNodes(); ++i) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
void
NeighborView<Dimension>::
setNodeExtents(const std::span<int>& nodeIDs) {
  const auto& H = mNodeListView.Hfield();
  for (auto i: nodeIDs) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
void
NeighborView<Dimension>::
setInternalNodeExtents() {
  const auto& H = mNodeListView.Hfield();
  const auto n = mNodeListView.numInternalNodes();
  for (auto i = 0u; i < n; ++i) {
    mNodeExtentsView(i) = Hextent(H(i), mKernelExtent);
  }
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
void
NeighborView<Dimension>::
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
template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<typename NeighborView<Dimension>::size_type, typename NeighborView<Dimension>::size_type>
NeighborView<Dimension>::
setMasterList(const int nodeID,
              std::span<int>& masterList,
              std::span<int>& coarseNeighbors,
              const bool ghostConnectivity) const {
  REQUIRE(nodeID >= 0 and nodeID < int(mNodeListView.numInternalNodes()));
  const auto& positions = mNodeListView.positions();
  const auto& Hfield = mNodeListView.Hfield();
  return this->setMasterList(positions(nodeID), Hfield(nodeID), masterList, coarseNeighbors, ghostConnectivity);
}

//------------------------------------------------------------------------------
// Set the refined list of potential neighbors based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
typename NeighborView<Dimension>::size_type
NeighborView<Dimension>::
setRefineNeighborList(const int nodeID,
                      const std::span<int>& coarseNeighbors,
                      std::span<int>& refineNeighbors) const {
  REQUIRE(nodeID >= 0 and nodeID < int(mNodeListView.numInternalNodes()));
  const auto& positions = mNodeListView.positions();
  const auto& Hfield = mNodeListView.Hfield();
  return this->setRefineNeighborList(positions(nodeID), Hfield(nodeID), coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// General precull routine.  Select nodes based on a range of positions and
// extents.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
typename NeighborView<Dimension>::size_type
NeighborView<Dimension>::
precullList(const Vector& minMasterPosition, const Vector& maxMasterPosition,
            const Vector& minMasterExtent, const Vector& maxMasterExtent,
            const typename NeighborView<Dimension>::size_type numCoarse,
            std::span<int>& coarseList) const {
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

}


