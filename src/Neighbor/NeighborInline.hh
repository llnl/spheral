#include "Utilities/DBC.hh"
#include <math.h>

namespace Spheral {

using std::abs;
using std::sqrt;

//------------------------------------------------------------------------------
// Calculate the maximum radial extent of a given smoothing tensor.
//------------------------------------------------------------------------------
// Make the general case the SPH approximation: assume the argument is a 
// scalar representing the inverse smoothing scale.
template<typename Dimension>
inline
typename Dimension::Vector
Neighbor<Dimension>::
HExtent(const typename Dimension::Scalar& H, 
        const double kernelExtent) {
  CHECK(H > 0.0);
  const double r = kernelExtent/H;
  return Vector(r);
}

//------------------------------------------------------------------------------
// Specializations for symmetric tensor smoothing transformation (ASPH).
// The ASPH tensor has units of inverse length.
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Vector
Neighbor< Dim<1> >::
HExtent(const Dim<1>::SymTensor& H,
        const double kernelExtent) {
  CHECK(H.Determinant() > 0.0);
  const double r = kernelExtent/H.xx();
  return Vector(r);
}

template<>
inline
Dim<2>::Vector
Neighbor< Dim<2> >::
HExtent(const Dim<2>::SymTensor& H,
        const double kernelExtent) {
  const double Hdet = H.Determinant();
  const SymTensor M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Vector(std::sqrt(M.yy()), sqrt(M.xx()));
}
  
template<>
inline
Dim<3>::Vector
Neighbor< Dim<3> >::
HExtent(const Dim<3>::SymTensor& H,
        const double kernelExtent) {
  const double Hdet = H.Determinant();
  const SymTensor M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Vector(sqrt(M.yy()*M.zz() - M.yz()*M.zy()),
                        sqrt(M.xx()*M.zz() - M.xz()*M.zx()),
                        sqrt(M.xx()*M.yy() - M.xy()*M.yx()));
}

//------------------------------------------------------------------------------
// Access the sampling kernel extent.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Neighbor<Dimension>::kernelExtent() const {
  return mKernelExtent;
}

template<typename Dimension>
inline
void
Neighbor<Dimension>::kernelExtent(double kernelExtent) {
  CHECK(kernelExtent > 0.0);
  mKernelExtent = kernelExtent;
  this->reinitialize();
}

//------------------------------------------------------------------------------
// Static method for setting master/coarse on multiple NodeLists/Neighbors.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIteratorType>
void
Neighbor<Dimension>::
setMasterNeighborGroup(const typename Dimension::Vector& position,
                       const typename Dimension::SymTensor& H,
                       const NodeListIteratorType& nodeListBegin,
                       const NodeListIteratorType& nodeListEnd,
                       double kernelExtent,
                       std::vector<std::vector<int>>& masterLists,
                       std::vector<std::vector<int>>& coarseNeighbors,
                       const bool ghostConnectivity) {

  typedef typename Dimension::Vector Vector;
  Vector minMasterPosition = Vector::one() * FLT_MAX;
  Vector maxMasterPosition = Vector::one() * -FLT_MAX;
  Vector minMasterExtent, maxMasterExtent;

  // Make sure the return values are sized appropriately.
  const auto numNodeLists = std::distance(nodeListBegin, nodeListEnd);
  masterLists = std::vector<std::vector<int>>(numNodeLists);
  coarseNeighbors = std::vector<std::vector<int>>(numNodeLists);

  // Loop over all the NodeLists.
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = **(nodeListBegin + nodeListi);
    const auto& neighbor = nodeList.neighbor();

    // Begin by setting the basic master/coarse neighbor info for each individual
    // NodeList.
    neighbor.setMasterList(position, H, masterLists[nodeListi], coarseNeighbors[nodeListi], ghostConnectivity);

    // Accumulate the min/max master node positions and extents.
    const auto& positions = nodeList.positions();
    const auto& extents = neighbor.nodeExtentField();
    const auto nmaster = masterLists[nodeListi].size();
    for (auto imaster = 0u; imaster < nmaster; ++imaster) {
      const auto  i = masterLists[nodeListi][imaster];
      const auto& ri = positions(i);
      const auto  minExtenti = ri - extents(i);
      const auto  maxExtenti = ri + extents(i);
      for (int k = 0; k < Dimension::nDim; ++k) {
        minMasterPosition(k) = std::min(minMasterPosition(k), ri(k));
        maxMasterPosition(k) = std::max(maxMasterPosition(k), ri(k));
        minMasterExtent(k) = std::min(minMasterExtent(k), minExtenti(k));
        maxMasterExtent(k) = std::max(maxMasterExtent(k), maxExtenti(k));
      }
    }
  }

  // Don't forget to compare these positions and extents with the given 
  // position/extent. :-P
  const auto extent = HExtent(H, kernelExtent);
  const auto minExtenti = position - extent;
  const auto maxExtenti = position + extent;
  for (int k = 0; k < Dimension::nDim; ++k) {
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
    coarseNeighbors[nodeListi] = neighbor.precullList(minMasterPosition, maxMasterPosition,
                                                      minMasterExtent, maxMasterExtent,
                                                      coarseNeighbors[nodeListi]);

//     // Set the per field coarse data caches for this NodeList.
//     (*nodeListItr)->notifyFieldsCacheCoarseValues();
  }
}

}
