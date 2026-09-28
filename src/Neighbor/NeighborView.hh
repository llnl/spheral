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

#include <tuple>
#include <cstdint>

namespace Spheral {

template<typename Dimension> class NodeListView;
template<typename Dimension, typename DataType> class FieldView;

template<typename Dimension, typename Descendant>
class NeighborView: public NeighborBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using DescendantType = Descendant;
  using size_type = int64_t;

  using NeighborBase<Dimension>::Hextent;

  // Constructors and destructors
  SPHERAL_HOST_DEVICE NeighborView(NodeListView<Dimension>& nodeList, 
                                   FieldView<Dimension, Vector>& nodeExtents,
                                   const NeighborSearchType searchType,
                                   const Scalar kernelExtent);
  SPHERAL_HOST_DEVICE virtual ~NeighborView() = default;

  // Allow access to the field of node extents.
  SPHERAL_HOST_DEVICE const FieldView<Dimension, Vector>& nodeExtentField() const { return mNodeExtentsView; }

  // Access the node list.
  SPHERAL_HOST_DEVICE const NodeListView<Dimension>& nodeListView() const         { return mNodeListView; }

  // Force the node extent field to be computed.
  SPHERAL_HOST_DEVICE void setNodeExtents();
  template<typename Container>
  SPHERAL_HOST_DEVICE void setNodeExtents(const Container& nodeIDs);
  SPHERAL_HOST_DEVICE void setInternalNodeExtents();
  SPHERAL_HOST_DEVICE void setGhostNodeExtents();

  // Set or refine the neighbor lists for a given node ID.
  // The std::spans should be sized with enough space before calling, and the final sizes are returned as the
  // size_type results.
  // If the std::spans are not allocated with sufficient space, then the required size is returned as the
  // ones complement (~size), which will be negative.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setMasterList(const int nodeID,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const;
  template<typename Container>
  SPHERAL_HOST_DEVICE size_type setRefineNeighborList(const int nodeID,
                                                      const Container& coarseNeighbors,
                                                      Container& refineNeighbors) const;

  // ********** Descendant NeighborView types must provide these methods. **********
  // Set or refine the neighbor lists for the given position and smoothing scale.
  // These methods return the sizes of the resulting arrays in tuples
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setMasterList(const Vector& position,
                                                                     const Scalar& H,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return static_cast<Descendant*>(this)->setMasterList(position, H, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setMasterList(const Vector& position,
                                                                     const SymTensor& H,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return static_cast<Descendant*>(this)->setMasterList(position, H, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setRefineNeighborList(const Vector& position,
                                                                             const Scalar& H,
                                                                             const Container& coarseNeighbors,
                                                                             Container& refineNeighbors) const {
    return static_cast<Descendant*>(this)->setRefineNeighborList(position, H, coarseNeighbors, refineNeighbors);
  }

  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setRefineNeighborList(const Vector& position,
                                                                             const SymTensor& H,
                                                                             const Container& coarseNeighbors,
                                                                             Container& refineNeighbors) const {
    return static_cast<Descendant*>(this)->setRefineNeighborList(position, H, coarseNeighbors, refineNeighbors);
  }

  // Set Neighbors for the given position.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setMasterList(const Vector& position,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return static_cast<Descendant*>(this)->setMasterList(position, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setRefineNeighborList(const Vector& position,
                                                                             const Container& coarseNeighbors,
                                                                             Container& refineNeighbors) const {
    return static_cast<Descendant*>(this)->setRefineNeighborList(position, coarseNeighbors, refineNeighbors);
  }

  // Set the neighbor lists based on proximity to planes.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                                     const GeomPlane<Dimension>& exitPlane,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors) const {
    return static_cast<Descendant*>(this)->setMasterList(enterPlane, exitPlane, masterList, coarseNeighbors);
  }

  // Force the update of internal data for the NodeList.
  SPHERAL_HOST_DEVICE virtual void updateNodes() = 0;
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<size_type, size_type> updateNodes(const Container& nodeIDs) { return static_cast<Descendant*>(this)->updateNodes(nodeIDs); }
  //****************************************************************************

  // Helper method to cull lists of neighbors based on min/max positions and min/max extents.
  // Note: this method modifies the passed set in-place, and returns the size
  template<typename Container>
  SPHERAL_HOST_DEVICE size_type precullList(const Vector& minMasterPosition, const Vector& maxMasterPostion,
                                            const Vector& minMasterExtent, const Vector& maxMasterExtent,
                                            const size_type numCoarse,
                                            Container& coarseList) const;

  // A static method for set master/coarse neighbor sets on multiple NodeLists
  // with culling of the coarse set.  This method is static because it is 
  // intended to deal with multiple Neighbor objects (not just this particular
  // one), and must access private data of the Neighbors.
  template<typename NodeListIterator, typename ContainerOfContainers, typename Container>
  SPHERAL_HOST_DEVICE static size_type setMasterNeighborGroup(const Vector& position,
                                                              const SymTensor& H,
                                                              const NodeListIterator& nodeListBegin,
                                                              const NodeListIterator& nodeListEnd,
                                                              const Scalar kernelExtent,
                                                              ContainerOfContainers& masterLists,
                                                              ContainerOfContainers& coarseNeighbors,
                                                              Container& numMasters,
                                                              Container& numCoarseNeighbors,
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
