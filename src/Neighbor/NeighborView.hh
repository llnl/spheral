//---------------------------------Spheral++----------------------------------//
// NeighborView -- The View class for Neighbor
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
// Split into Value/View pattern Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborView__
#define __Spheral_NeighborView__

#include "Neighbor/NeighborSearchType.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"

#include <utility>    // std::pair

namespace Spheral {

template<typename Dimension> class NodeListView;
template<typename Dimension, typename DataType> class FieldView;
template<typename Dimension, typename DataType> class FieldListView;
template<typename Dimension> class DataBaseView;

template<typename Dimension, typename Descendant>
class NeighborView {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using NodeID = std::pair<size_t, size_t>;  // (NodeList ID, node ID)

  // Constructors and destructors
  SPHERAL_HOST_DEVICE NeighborView(DataBaseView<Dimension> dataBase,
                                   const NeighborSearchType searchType,
                                   const Scalar kernelExtent);
  SPHERAL_HOST_DEVICE NeighborView(const NeighborView& rhs) = default;
  SPHERAL_HOST_DEVICE NeighborView(NeighborView&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~NeighborView() = default;

  // Access our data
  SPHERAL_HOST_DEVICE const DataBaseView<Dimension>& dataBaseView()              const { return mDataBaseView; }
  SPHERAL_HOST_DEVICE const FieldListView<Dimension, Vector>& nodeExtentView()   const { return mNodeExtentView; }
  SPHERAL_HOST_DEVICE NeighborSearchType neighborSearchType()                    const { return mSearchType; }
  SPHERAL_HOST_DEVICE Scalar kernelExtent()                                      const { return mKernelExtent; }

  // Force the node extent field to be computed.
  SPHERAL_HOST_DEVICE void setNodeExtents(const DataBaseView<Dimension>& dataBase);
  SPHERAL_HOST_DEVICE void setInternalNodeExtents(const DataBaseView<Dimension>& dataBase);
  SPHERAL_HOST_DEVICE void setGhostNodeExtents(const DataBaseView<Dimension>& dataBase);
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE void setNodeExtents(const DataBaseView<Dimension>& dataBase,
                                          const Container<size_t>& nodeIDs);

  // Helper method to cull lists of neighbors based on min/max positions and min/max extents.
  // Note: this method modifies the passed coarseList in-place, and returns the resulting size
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE size_t precullList(const FieldListView<Dimension, Vector>& positions,
                                         const Vector& minMasterPosition, const Vector& maxMasterPostion,
                                         const Vector& minMasterExtent, const Vector& maxMasterExtent,
                                         const size_t numCoarse,
                                         Container<NodeID>& coarseList) const;

  // Determine the maximum extent of a given H smoothing scale along the
  // Cartesian axes.
  SPHERAL_HOST_DEVICE static Vector Hextent(const Scalar& H, const Scalar& kernelExtent);
  SPHERAL_HOST_DEVICE static Vector Hextent(const SymTensor& H, const Scalar& kernelExtent);

  // Cast to the descendant type
  SPHERAL_HOST_DEVICE Descendant* asDescendant()                                       { return static_cast<Descendant*>(*this); }

  //.............................................................................
  // Set or refine the neighbor lists
  // The Containers should be sized with enough space before calling, and the final sizes are returned as the
  // size_type results.
  // If the std::spans are not allocated with sufficient space, then the required size is returned as the
  // result, and the status (first bool) is false.
  // In the refine nodes nodes we assume the refineNeighborViews is at least as large as coarseNeighborViews, so there
  // should always be enough room and the bool status flag is not needed.

  // Set neighbors given (position, H)
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                                     const SymTensor& H,
                                                                     Container<NodeID>& masterList,
                                                                     Container<NodeID>& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return asDescendant()->setMasterList(position, H, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const Vector& position,
                                                   const SymTensor& H,
                                                   const size_t numCoarse,
                                                   const Container<NodeID>& coarseNeighbors,
                                                   Container<NodeID>& refineNeighbors) const {
    return asDescendant()->setRefineNeighborList(position, H, numCoarse, coarseNeighbors, refineNeighbors);
  }

  // Set Neighbors for the given position
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                                     Container<NodeID>& masterList,
                                                                     Container<NodeID>& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const {
    return asDescendant()->setMasterList(position, masterList, coarseNeighbors, ghostConnectivity);
  }

  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const Vector& position,
                                                   const size_t numCoarse,
                                                   const Container<NodeID>& coarseNeighbors,
                                                   Container<NodeID>& refineNeighbors) const {
    return asDescendant()->setRefineNeighborList(position, numCoarse, coarseNeighbors, refineNeighbors);
  }

  // Set the neighbor lists based on proximity to planes.
  template<template<typename ...> typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                                     const GeomPlane<Dimension>& exitPlane,
                                                                     Container<NodeID>& masterList,
                                                                     Container<NodeID>& coarseNeighbors) const {
    return asDescendant()->setMasterList(enterPlane, exitPlane, masterList, coarseNeighbors);
  }

  //.............................................................................
  // Requirements for rebuilding the internal data of a NeighborView object.
  // Expected to be CPU Value object tasks for now.
  // Force the update of internal data for the NodeList.
  SPHERAL_HOST_DEVICE void updateNodes(const DataBaseView<Dimension>& dataBase) { asDescendant()->updateNodes(dataBase); }
  SPHERAL_HOST_DEVICE void updateNodes(const DataBaseView<Dimension>& dataBase,
                                       const std::span<NodeID>& nodeIDs)        { asDescendant()->updateNodes(dataBase, nodeIDs); }

protected:
  //-------------------------- Protected Interface --------------------------//
  DataBaseView<Dimension> mDataBaseView;
  FieldListView<Dimension, Vector> mNodeExtentView;
  NeighborSearchType mSearchType;
  Scalar mKernelExtent;
};

}

#include "Neighbor/NeighborViewInline.hh"

#endif
