//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
// Split into Value/View pattern Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_Neighbor__
#define __Spheral_Neighbor__

#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

template<typename Dimension> class NodeListView;
template<typename Dimension, typename DataType> class FieldView;

enum class NeighborSearchType {
  Gather = 1,
  Scatter = 2,
  GatherScatter = 3
};

template<typename Dimension, typename Descendant>
class Neighbor {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors and destructors
  SPHERAL_HOST_DEVICE Neighbor(const NeighborSearchType searchType,
                               const Scalar kernelExtent);
  SPHERAL_HOST_DEVICE Neighbor(const Neighbor& rhs) = default;
  SPHERAL_HOST_DEVICE Neighbor(Neighbor&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~Neighbor() = default;

  // Access our data
  SPHERAL_HOST_DEVICE const FieldView<Dimension, Vector>& nodeExtentView()   const { return mNodeExtentView; }
  SPHERAL_HOST_DEVICE NeighborSearchType neighborSearchType()                const { return mSearchType; }
  SPHERAL_HOST_DEVICE Scalar kernelExtent()                                  const { return mKernelExtent; }

  // Set internal data
  void neighborSearchType(const NeighborSearchType& x)                             { mSearchType = x; }
  void kernelExtent(const Scalar x)                                                { VERIFY(x > 0.0); mKernelExtent = x; }

  // Force the node extent field to be computed.
  SPHERAL_HOST_DEVICE void setNodeExtents();
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE void setNodeExtents(const Container<size_t>& nodeIDs);
  SPHERAL_HOST_DEVICE void setInternalNodeExtents();
  SPHERAL_HOST_DEVICE void setGhostNodeExtents();

  // Helper method to cull lists of neighbors based on min/max positions and min/max extents.
  // Note: this method modifies the passed coarseList in-place, and returns the resulting size
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE size_t precullList(const Vector& minMasterPosition, const Vector& maxMasterPostion,
                                         const Vector& minMasterExtent, const Vector& maxMasterExtent,
                                         const size_t numCoarse,
                                         Container<size_t>& coarseList) const;

  // A static method for set master/coarse neighbor sets on multiple NodeLists
  // with culling of the coarse set.  This method is static because it is 
  // intended to deal with multiple Neighbor objects (not just this particular
  // one), and must access private data of the Neighbors.
  // Returns a bool indicating success if there was a enough space in the return arrays
  // The numMasters/numCoarseNeighbors indicate the actual sizes of masters and neighbors required
  template<typename NodeListIterator,typename ContainerOfContainers, typename Container>
  SPHERAL_HOST_DEVICE static bool setMasterNeighborGroup(const Vector& position,
                                                         const SymTensor& H,
                                                         const NodeListIterator& nodeListBegin,
                                                         const NodeListIterator& nodeListEnd,
                                                         const Scalar kernelExtent,
                                                         ContainerOfContainers& masterLists,
                                                         ContainerOfContainers& coarseNeighbors,
                                                         Container& numMasters,
                                                         Container& numCoarseNeighbors,
                                                         const bool ghostConnectivity = false);

  // Determine the maximum extent of a given H smoothing scale along the
  // Cartesian axes.
  SPHERAL_HOST_DEVICE static Vector Hextent(const Scalar& H, const Scalar& kernelExtent);
  SPHERAL_HOST_DEVICE static Vector Hextent(const SymTensor& H, const Scalar& kernelExtent);

  // Cast to the descendant type
  SPHERAL_HOST_DEVICE Descendant* asDescendant()                                   { return static_cast<Descendant*>(this); }

  //.............................................................................
  // Set or refine the neighbor lists for a given node ID.
  // The Containers should be sized with enough space before calling, and the final sizes are returned as the
  // size_type results.
  // If the std::spans are not allocated with sufficient space, then the required size is returned as the
  // result, and the status (first bool) is false.
  // In the refine nodes nodes we assume the refineNeighbors is at least as large as coarseNeighbors, so there
  // should always be enough room and the bool status flag is not needed.
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const size_t nodeID,
                                                                     Container<size_t>& masterList,
                                                                     Container<size_t>& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const;
    
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const size_t nodeID,
                                                   const size_t numCoarse,
                                                   const Container<size_t>& coarseNeighbors,
                                                   Container<size_t>& refineNeighbors) const;

  //.............................................................................
  // The required interface for Neighbor objects (device or CPU)
  // Set neighbors given (position, H) or (posiion, h)
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                                     const SymTensor& H,
                                                                     Container<size_t>& masterList,
                                                                     Container<size_t>& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return asDescendant()->setMasterList(position, H, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<template<typename ...> typename Container>
        SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const Vector& position,
                                                         const SymTensor& H,
                                                         const size_t numCoarse,
                                                         const Container<size_t>& coarseNeighbors,
                                                         Container<size_t>& refineNeighbors) const {
    return asDescendant()->setRefineNeighborList(position, H, numCoarse, coarseNeighbors, refineNeighbors);
  }

  // Set Neighbors for the given position
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                                     Container<size_t>& masterList,
                                                                     Container<size_t>& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return asDescendant()->setMasterList(position, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const Vector& position,
                                                   const size_t numCoarse,
                                                   const Container<size_t>& coarseNeighbors,
                                                   Container<size_t>& refineNeighbors) const {
    return asDescendant()->setRefineNeighborList(position, numCoarse, coarseNeighbors, refineNeighbors);
  }

  // Set the neighbor lists based on proximity to planes.
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                                     const GeomPlane<Dimension>& exitPlane,
                                                                     Container<size_t>& masterList,
                                                                     Container<size_t>& coarseNeighbors) const {
    return asDescendant()->setMasterList(enterPlane, exitPlane, masterList, coarseNeighbors);
  }

  // //.............................................................................
  // // Requirements for rebuilding the internal data of a Neighbor object.
  // // Expected to be CPU Value object tasks for now.
  // // Force the update of internal data for the NodeList.
  // virtual void updateNodes() = 0;
  // void updateNodes(const std::span<size_t>& nodeIDs) = 0;

  // // Reinitialize based on the desired target h and box dimensions.
  // virtual void reinitialize() = 0;
  // virtual void reinitialize(const Vector& xmin,
  //                           const Vector& xmax,
  //                           const Scalar htarget) = 0;


protected:
  //-------------------------- Protected Interface --------------------------//
  FieldView<Dimension, Vector>& mNodeExtentView;
  NeighborSearchType mSearchType;
  Scalar mKernelExtent;
};

}

#include "Neighbor/NeighborInline.hh"

#endif
