#include "Utilities/DBC.hh"

#include <algorithm>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Static method for setting master/coarse on multiple NodeLists/Neighbors.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIteratorType>
SPHERAL_HOST_DEVICE 
typename NeighborView<Dimension>::size_type                        // num NodeLists
NeighborView<Dimension>::
setMasterNeighborGroup(const typename Dimension::Vector& position,
                       const typename Dimension::SymTensor& H,
                       const NodeListIteratorType& nodeListBegin,
                       const NodeListIteratorType& nodeListEnd,
                       const Scalar kernelExtent,
                       std::span<std::span<int>>& masterLists,
                       std::span<std::span<int>>& coarseNeighbors,
                       std::span<typename NeighborView<Dimension>::size_type>& numMasters,
                       std::span<typename NeighborView<Dimension>::size_type>& numCoarseNeighbors,
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
