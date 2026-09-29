//---------------------------------Spheral++----------------------------------//
// PythonNeighborView
//
// PythonNeighborView just provides the virtual function hooks to support the
// NeighborView interface.
//
// NOTE: PythonNeighbor and its descendants are CPU-only (no GPU capability).
//
// Created by J. Michael Owen, Thu Sep 24 15:41:42 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_PythonNeighborView__
#define __Spheral_PythonNeighborView__

#include "Neighbor/Neighbor.hh"

namespace Spheral {

template<typename Dimension>
class PythonNeighborView:
    public Neighbor<Dimension, PythonNeighborView<Dimension>> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using Container = std::vector<int>;
  using ContainerOfContainers = std::vector<std::vector<size_t>>;

  PythonNeighborView(NodeListView<Dimension>& nodeListView,
                     FieldView<Dimension, Vector>& nodeExtentView,
                     const NeighborSearchType searchType,
                     const Scalar kernelExtent);
  virtual ~PythonNeighborView() = default;

  //.............................................................................
  // The required interface for Neighbor objects
  virtual std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                         const Scalar& H,
                                                         Container& masterList,
                                                         Container& coarseNeighbors,
                                                         const bool ghostConnectivity = false) const = 0;

  virtual std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                         const SymTensor& H,
                                                         Container& masterList,
                                                         Container& coarseNeighbors,
                                                         const bool ghostConnectivity = false) const = 0;

  virtual std::tuple<bool, size_t, size_t> setRefineNeighborList(const Vector& position,
                                                                 const Scalar& H,
                                                                 const Container& coarseNeighbors,
                                                                 Container& refineNeighbors) const = 0;

  virtual std::tuple<bool, size_t, size_t> setRefineNeighborList(const Vector& position,
                                                                 const SymTensor& H,
                                                                 const Container& coarseNeighbors,
                                                                 Container& refineNeighbors) const = 0;

  // Set Neighbors for the given position.
  virtual std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                         Container& masterList,
                                                         Container& coarseNeighbors,
                                                         const bool ghostConnectivity = false) const = 0;

  virtual std::tuple<bool, size_t, size_t> setRefineNeighborList(const Vector& position,
                                                                 const Container& coarseNeighbors,
                                                                 Container& refineNeighbors) const = 0;

  // Set the neighbor lists based on proximity to planes.
  virtual std::tuple<bool, size_t, size_t> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                         const GeomPlane<Dimension>& exitPlane,
                                                         Container& masterList,
                                                         Container& coarseNeighbors) const = 0;

  // Force the update of internal data for the NodeList.
  virtual void updateNodes() = 0;
  virtual void updateNodes(const std::vector<int>& nodeIDs) = 0;

  // Reinitialize based on the desired target h and box dimensions.
  virtual void reinitialize() = 0;
  virtual void reinitialize(const Vector& xmin, const Vector& xmax, const Scalar htarget) = 0;
};

}

#endif
