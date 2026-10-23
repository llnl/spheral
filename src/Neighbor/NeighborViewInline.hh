#include "DataBase/DataBaseView.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
#include <limits>

namespace Spheral {

namespace { // anonymous

//------------------------------------------------------------------------------
// Specializations for symmetric tensor smoothing transformation (ASPH).
// The ASPH tensor has units of inverse length.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE
inline
Dim<1>::Vector
HextentImpl(const Dim<1>::SymTensor& H,
            const Dim<1>::Scalar& kernelExtent) {
  REQUIRE(H.Determinant() > 0.0);
  return Dim<1>::Vector(kernelExtent/H.xx());
}

SPHERAL_HOST_DEVICE
inline
Dim<2>::Vector
HextentImpl(const Dim<2>::SymTensor& H,
            const Dim<2>::Scalar& kernelExtent) {
  const auto Hdet = H.Determinant();
  const auto M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Dim<2>::Vector(std::sqrt(M.yy()), std::sqrt(M.xx()));
}
  
SPHERAL_HOST_DEVICE
inline
Dim<3>::Vector
HextentImpl(const Dim<3>::SymTensor& H,
            const Dim<3>::Scalar& kernelExtent) {
  const auto Hdet = H.Determinant();
  const auto M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Dim<3>::Vector(std::sqrt(M.yy()*M.zz() - M.yz()*M.zy()),
                                          std::sqrt(M.xx()*M.zz() - M.xz()*M.zx()),
                                          std::sqrt(M.xx()*M.yy() - M.xy()*M.yx()));
}

}             // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE
inline
NeighborView<Dimension, Descendant>::
NeighborView(DataBaseView<Dimension> dataBase,
             const NeighborSearchType searchType,
             const Scalar kernelExtent):
  mDataBaseView(dataBase),
  mNodeExtentView(),
  mSearchType(searchType),
  mKernelExtent(kernelExtent) {
}

//------------------------------------------------------------------------------
// Force the node extent field to be computed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setNodeExtents(const DataBaseView<Dimension>& dataBase) {
  REQUIRE(mNodeExtentView.numFields() == dataBase.numNodeLists());
  for (auto [k, nv]: enumerate(dataBase.nodeListViews())) {
    const auto& H = nv.Hfield();
    for (auto i = 0u; i < nv.numNodes(); ++i) {
      mNodeExtentView(k,i) = Hextent(H(i), mKernelExtent);
    }
  }
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setInternalNodeExtents(const DataBaseView<Dimension>& dataBase) {
  REQUIRE(mNodeExtentView.numFields() == dataBase.numNodeLists());
  for (auto [k, nv]: enumerate(dataBase.nodeListViews())) {
    const auto& H = nv.Hfield();
    for (auto i = 0u; i < nv.numInternalNodes(); ++i) {
      mNodeExtentView(k,i) = Hextent(H(i), mKernelExtent);
    }
  }
}

template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setGhostNodeExtents(const DataBaseView<Dimension>& dataBase) {
  REQUIRE(mNodeExtentView.numFields() == dataBase.numNodeLists());
  for (auto [k, nv]: enumerate(dataBase.nodeListViews())) {
    const auto& H = nv.Hfield();
    for (auto i = nv.firstGhostNode(); i < nv.numNodes(); ++i) {
      mNodeExtentView(k,i) = Hextent(H(i), mKernelExtent);
    }
  }
}

template<typename Dimension, typename Descendant>
template<template<typename ...> typename Container>
SPHERAL_HOST_DEVICE 
inline
void
NeighborView<Dimension, Descendant>::
setNodeExtents(const DataBaseView<Dimension>& dataBase,
               const Container<size_t>& nodeIDs) {
  REQUIRE(mNodeExtentView.numFields() == dataBase.numNodeLists());
  const auto& nodeListViews = dataBase.nodeListViews();
  for (auto [k, i]: nodeIDs) {
    const auto& Hi = nodeListViews[k].Hfield()(i);
    mNodeExtentView(k,i) = Hextent(Hi, mKernelExtent);
  }
}

//------------------------------------------------------------------------------
// General precull routine.  Select nodes based on a range of positions and
// extents.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
template<template<typename ...> typename Container>
SPHERAL_HOST_DEVICE 
size_t
NeighborView<Dimension, Descendant>::
precullList(const FieldListView<Dimension, Vector>& positions,
            const Vector& minMasterPosition, const Vector& maxMasterPosition,
            const Vector& minMasterExtent, const Vector& maxMasterExtent,
            const size_t numCoarse,
            Container<NodeID>& coarseList) const {
  REQUIRE(numCoarse <= coarseList.size());

  // What kind of preculling we're doing determines the applied test.
  size_t result = 0u;
  switch(mSearchType) {

  case NeighborSearchType::GatherScatter:
    for (auto k = 0u; k < numCoarse; ++k) {
      const auto& [j, i] = coarseList[k];
      const auto& nodePosition = positions(j,i);
      const auto  minNodeExtent = nodePosition - mNodeExtentView(j,i);
      const auto  maxNodeExtent = nodePosition + mNodeExtentView(j,i);
      if (testPointInBox(nodePosition,
                         minMasterExtent,
                         maxMasterExtent) or
          testBoxIntersection(minMasterPosition,
                              maxMasterPosition,
                              minNodeExtent,
                              maxNodeExtent)) {
        coarseList[result] = std::make_pair(j,i);
        ++result;
      }
    }
    break;

  case NeighborSearchType::Gather:
    for (auto k = 0u; k < numCoarse; ++k) {
      const auto& [j, i] = coarseList[k];
      const auto& nodePosition = positions(j,i);
      if (testPointInBox(nodePosition,
                         minMasterExtent,
                         maxMasterExtent)) {
        coarseList[result] = std::make_pair(j,i);
        ++result;
      }
    }
    break;

  case NeighborSearchType::Scatter:
    for (auto k = 0u; k < numCoarse; ++k) {
      const auto& [j, i] = coarseList[k];
      const auto& nodePosition = positions(j,i);
      const auto  minNodeExtent = nodePosition - mNodeExtentView(j,i);
      const auto  maxNodeExtent = nodePosition + mNodeExtentView(j,i);
      if (testBoxIntersection(minMasterPosition,
                              maxMasterPosition,
                              minNodeExtent,
                              maxNodeExtent)) {
        coarseList[result] = std::make_pair(j,i);
        ++result;
      }
    }
    break;

  default:
    VERIFY2(false, "Neighbor::precullList Unhandeled NeighborSearchType");
  }

  ENSURE(result <= numCoarse);
  return result;
}

// //------------------------------------------------------------------------------
// // Static method for setting master/coarse on multiple NodeLists/Neighbors.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename Descendant>
// template<typename NodeListIterator, typename ContainerOfContainers, typename Container>
// SPHERAL_HOST_DEVICE 
// bool                       // success/fail
// Neighbor<Dimension, Descendant>::
// setMasterNeighborGroup(const typename Dimension::Vector& position,
//                        const typename Dimension::SymTensor& H,
//                        const NodeListIterator& nodeListBegin,
//                        const NodeListIterator& nodeListEnd,
//                        const Scalar kernelExtent,
//                        ContainerOfContainers& masterLists,
//                        ContainerOfContainers& coarseNeighbors,
//                        Container& numMasters,
//                        Container& numCoarseNeighbors,
//                        const bool ghostConnectivity) {

//   // Pre-conditions
//   const auto numNodeLists = std::distance(nodeListBegin, nodeListEnd);
//   REQUIRE(masterLists.size() == numNodeLists and
//           coarseNeighbors.size() == numNodeLists and
//           numMasters.size() == numNodeLists and
//           numCoarseNeighbors.size() == numNodeLists);

//   // Initializations
//   Vector minMasterPosition = Vector::one() * std::numeric_limits<Scalar>::max();
//   Vector maxMasterPosition = Vector::one() * std::numeric_limits<Scalar>::lowest();
//   Vector minMasterExtent, maxMasterExtent;

//   // Loop over all the NodeLists.
//   bool result = true;
//   for (auto [nodeListi, nodeListPtr]: enumerate(nodeListBegin, nodeListEnd)) { 
//     const auto& neighbor = nodeListPtr->neighborView();

//     // Begin by setting the basic master/coarse neighbor info for each individual NodeList.
//     bool status;
//     std::tie(status, numMasters[nodeListi], numCoarseNeighbors[nodeListi]) = neighbor.setMasterList(position, H, masterLists[nodeListi], coarseNeighbors[nodeListi], ghostConnectivity);
//     result &= status;

//     // Accumulate the min/max master node positions and extents.
//     if (result) {
//       const auto& positions = nodeListPtr->positions();
//       const auto& extents = neighbor.nodeExtentField();
//       for (auto imaster = 0; imaster < numMasters[nodeListi]; ++imaster) {
//         const auto  i = masterLists[nodeListi][imaster];
//         const auto& ri = positions(i);
//         const auto  minExtenti = ri - extents(i);
//         const auto  maxExtenti = ri + extents(i);
//         for (auto k = 0u; k < Dimension::nDim; ++k) {
//           minMasterPosition(k) = std::min(minMasterPosition(k), ri(k));
//           maxMasterPosition(k) = std::max(maxMasterPosition(k), ri(k));
//           minMasterExtent(k) = std::min(minMasterExtent(k), minExtenti(k));
//           maxMasterExtent(k) = std::max(maxMasterExtent(k), maxExtenti(k));
//         }
//       }
//     }
//   }

//   // Don't forget to compare these positions and extents with the given position/extent. :-P
//   if (result) {
//     const auto extent = HExtent(H, kernelExtent);
//     const auto minExtenti = position - extent;
//     const auto maxExtenti = position + extent;
//     for (auto k = 0u; k < Dimension::nDim; ++k) {
//       minMasterPosition(k) = std::min(minMasterPosition(k), position(k));
//       maxMasterPosition(k) = std::max(maxMasterPosition(k), position(k));
//       minMasterExtent(k) = std::min(minMasterExtent(k), minExtenti(k));
//       maxMasterExtent(k) = std::max(maxMasterExtent(k), maxExtenti(k));
//     }

//     // Loop over the nodes again, and cull the coarse neighbor lists according
//     // to the overall min/max master distribution we just calcuated.
//     // WARNING!  After this step the neighbor information in the NodeLists is only
//     // guaranteed complete for the set of NodeLists passed to this method!
//     for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
//       const auto& nodeList = **(nodeListBegin + nodeListi);
//       const auto& neighbor = nodeList.neighbor();
//       numCoarseNeighbors[nodeListi] = neighbor.precullList(minMasterPosition, maxMasterPosition,
//                                                            minMasterExtent, maxMasterExtent,
//                                                            coarseNeighbors[nodeListi]);
//     }
//   }

//   return result;
// }

//------------------------------------------------------------------------------
// Calculate the maximum radial extent of a given smoothing tensor.
//------------------------------------------------------------------------------
// Make the general case the SPH approximation: assume the argument is a 
// scalar representing the inverse smoothing scale.
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE 
inline
typename Dimension::Vector
NeighborView<Dimension, Descendant>::
Hextent(const typename Dimension::Scalar& H, 
        const typename Dimension::Scalar& kernelExtent) {
  REQUIRE(H > 0.0);
  return Vector(kernelExtent/H,
                kernelExtent/H,
                kernelExtent/H);
}

//------------------------------------------------------------------------------
// Specializations for symmetric tensor smoothing transformation (ASPH).
// The ASPH tensor has units of inverse length.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendant>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Vector
NeighborView<Dimension, Descendant>::
Hextent(const SymTensor& H,
        const Scalar& kernelExtent) {
  return HextentImpl(H, kernelExtent);
}

}
