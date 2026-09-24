//---------------------------------Spheral++----------------------------------//
// NeighborView -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborView__
#define __Spheral_NeighborView__

#include "Neighbor/NeighborBase.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"

#include <span>
#include <tuple>
#include <cstdint>

namespace Spheral {

template<typename Dimension> class NodeListView;
template<typename Dimension, typename DataType> class FieldView;

template<typename Dimension>
class NeighborView: public NeighborBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using size_type = int64_t;

  using NeighborBase<Dimension>::Hextent;

  // Constructors and destructors
  SPHERAL_HOST_DEVICE NeighborView(NodeListView<Dimension>& nodeList, 
                                   FieldView<Dimension, Vector>& nodeExtents,
                                   const NeighborSearchType searchType,
                                   const double kernelExtent);
  SPHERAL_HOST_DEVICE virtual ~NeighborView() = default;

  // Allow access to the field of node extents.
  SPHERAL_HOST_DEVICE const FieldView<Dimension, Vector>& nodeExtentField() const { return mNodeExtentsView; }

  // Access the node list.
  SPHERAL_HOST_DEVICE const NodeListView<Dimension>& nodeListView() const         { return mNodeListView; }

  // Force the node extent field to be computed.
  SPHERAL_HOST_DEVICE void setNodeExtents();
  SPHERAL_HOST_DEVICE void setNodeExtents(const std::span<int>& nodeIDs);
  SPHERAL_HOST_DEVICE void setInternalNodeExtents();
  SPHERAL_HOST_DEVICE void setGhostNodeExtents();

  // Set or refine the neighbor lists for a given node ID.
  // The std::spans should be sized with enough space before calling, and the final sizes are returned as the
  // size_type results.
  // If the std::spans are not allocated with sufficient space, then the required size is returned as the
  // ones complement (~size), which will be negative.
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setMasterList(const int nodeID,
                                                                             std::span<int>& masterList,
                                                                             std::span<int>& coarseNeighbors,
                                                                             const bool ghostConnectivity = false) const;
  SPHERAL_HOST_DEVICE virtual size_type setRefineNeighborList(const int nodeID,
                                                              const std::span<int>& coarseNeighbors,
                                                              std::span<int>& refineNeighbors) const;

  // ********** Descendent NeighborView types must provide these methods. **********
  // Set or refine the neighbor lists for the given position and smoothing scale.
  // These methods return the sizes of the resulting arrays in tuples
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setMasterList(const Vector& position,
                                                                             const Scalar& H,
                                                                             std::span<int>& masterList,
                                                                             std::span<int>& coarseNeighbors,
                                                                             const bool ghostConnectivity = false) const = 0;
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setMasterList(const Vector& position,
                                                                             const SymTensor& H,
                                                                             std::span<int>& masterList,
                                                                             std::span<int>& coarseNeighbors,
                                                                             const bool ghostConnectivity = false) const = 0;

  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setRefineNeighborList(const Vector& position,
                                                                                     const Scalar& H,
                                                                                     const std::span<int>& coarseNeighbors,
                                                                                     std::span<int>& refineNeighbors) const = 0;
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setRefineNeighborList(const Vector& position,
                                                                                     const SymTensor& H,
                                                                                     const std::span<int>& coarseNeighbors,
                                                                                     std::span<int>& refineNeighbors) const = 0;

  // Set Neighbors for the given position.
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setMasterList(const Vector& position,
                                                                             std::span<int>& masterList,
                                                                             std::span<int>& coarseNeighbors,
                                                                             const bool ghostConnectivity = false) const = 0;
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setRefineNeighborList(const Vector& position,
                                                                                     const std::span<int>& coarseNeighbors,
                                                                                     std::span<int>& refineNeighbors) const = 0;

  // Set the neighbor lists based on proximity to planes.
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                                             const GeomPlane<Dimension>& exitPlane,
                                                                             std::span<int>& masterList,
                                                                             std::span<int>& coarseNeighbors) const = 0;

  // Force the update of internal data for the NodeList.
  SPHERAL_HOST_DEVICE virtual void updateNodes() = 0;
  SPHERAL_HOST_DEVICE virtual std::tuple<size_type, size_type> updateNodes(const std::span<int>& nodeIDs) = 0;
  //****************************************************************************

  // Helper method to cull lists of neighbors based on min/max positions and min/max extents.
  // Note: this method modifies the passed set in-place, and returns the size
  SPHERAL_HOST_DEVICE size_type precullList(const Vector& minMasterPosition, const Vector& maxMasterPostion,
                                            const Vector& minMasterExtent, const Vector& maxMasterExtent,
                                            const size_type numCoarse,
                                            std::span<int>& coarseList) const;

  // A static method for set master/coarse neighbor sets on multiple NodeLists
  // with culling of the coarse set.  This method is static because it is 
  // intended to deal with multiple Neighbor objects (not just this particular
  // one), and must access private data of the Neighbors.
  template<typename NodeListIteratorType>
  SPHERAL_HOST_DEVICE static size_type setMasterNeighborGroup(const Vector& position,
                                                              const SymTensor& H,
                                                              const NodeListIteratorType& nodeListBegin,
                                                              const NodeListIteratorType& nodeListEnd,
                                                              const Scalar kernelExtent,
                                                              std::span<std::span<int>>& masterLists,
                                                              std::span<std::span<int>>& coarseNeighbors,
                                                              std::span<size_type>& numMasters,
                                                              std::span<size_type>& numCoarseNeighbors,
                                                              const bool ghostConnectivity = false);

private:
  //--------------------------- Private Interface ---------------------------//
  using NeighborBase<Dimension>::mSearchType;
  using NeighborBase<Dimension>::mKernelExtent;
  NodeListView<Dimension> mNodeListView;
  FieldView<Dimension, Vector> mNodeExtentsView;
};

}

#include "NeighborViewInline.hh"

#endif
