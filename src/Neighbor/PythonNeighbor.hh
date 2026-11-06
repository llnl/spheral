//---------------------------------Spheral++----------------------------------//
// PythonNeighbor
//
// PythonNeighbor just provides the virtual function hooks to support the
// NeighborView interface.  The virutal Neighbor interface is passed along
// for descendants to provide.
//
// NOTE: PythonNeighbor and its descendants are CPU-only (no GPU capability).
//
// Created by J. Michael Owen, Thu Sep 24 15:41:42 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_PythonNeighbor__
#define __Spheral_PythonNeighbor__

#include "Neighbor/Neighbor.hh"
#include "Neighbor/NeighborView.hh"
#include "Neighbor/NeighborViewVariant.hh"
#include "Neighbor/PythonNeighborCallView.hh"

#include <utility>
#include <vector>

namespace Spheral {

// Forward declarations
template<typename Dimension> class GeomPlane;
template<typename Dimension> class DataBase;

template<typename Dimension>
class PythonNeighbor:
    public Neighbor<Dimension>,
    public NeighborView<Dimension, PythonNeighbor<Dimension>> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using BaseValueType = Neighbor<Dimension>;
  using BaseViewType  = NeighborView<Dimension, PythonNeighbor<Dimension>>;
  using ViewType = PythonNeighborCallView<Dimension>;

  using NodeID = typename BaseValueType::NodeID;
  template<typename T> using Container = typename BaseViewType::template Container<T>;

  PythonNeighbor(DataBase<Dimension>& dataBase,
                 const NeighborSearchType searchType,
                 const Scalar kernelExtent);
  virtual ~PythonNeighbor() = default;

  //.............................................................................
  // The required Neighbor interface
  virtual void resizeNodeExtents() override;

  // Clear the current Neighbor state
  virtual void clear() override;

  // Return a lightweight callback adapter.  The adapter does not own this
  // PythonNeighbor, so this object must outlive all copies of the view.
  ViewType view() { return ViewType(this); }
  NeighborViewVariant<Dimension> variantView() { return this->view(); }

  //.............................................................................
  // The required interface for NeighborView objects
  virtual std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                         const SymTensor& H,
                                                         Container<NodeID>& masterList,
                                                         Container<NodeID>& coarseNeighbors,
                                                         const bool ghostConnectivity = false) const = 0;

  virtual size_t setRefineNeighborList(const Vector& position,
                                       const SymTensor& H,
                                       const size_t numCoarse,
                                       const Container<NodeID>& coarseNeighbors,
                                       Container<NodeID>& refineNeighbors) const = 0;

  // Set Neighbors for the given position.
  virtual std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                         Container<NodeID>& masterList,
                                                         Container<NodeID>& coarseNeighbors,
                                                         const bool ghostConnectivity = false) const = 0;

  virtual size_t setRefineNeighborList(const Vector& position,
                                       const size_t numCoarse,
                                       const Container<NodeID>& coarseNeighbors,
                                       Container<NodeID>& refineNeighbors) const = 0;

  // Set the neighbor lists based on proximity to planes.
  virtual std::tuple<bool, size_t, size_t> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                         const GeomPlane<Dimension>& exitPlane,
                                                         Container<NodeID>& masterList,
                                                         Container<NodeID>& coarseNeighbors) const = 0;

private:
  //-------------------------- Private Interface --------------------------//
  using BaseValueType::mDataBase;
  using BaseValueType::mNodeExtent;
  using BaseViewType::mNodeExtentView;
};

}

#endif
