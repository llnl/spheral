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

namespace Spheral {

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

  using NodeID = typename Neighbor<Dimension>::NodeID;
  using Container = std::vector<NodeID>;

  PythonNeighbor(DataBase<Dimension>& dataBase,
                 const NeighborSearchType searchType,
                 const Scalar kernelExtent);
  virtual ~PythonNeighbor() = default;

  //.............................................................................
  // The required interface for NeighborView objects
  virtual std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                         const SymTensor& H,
                                                         Container& masterList,
                                                         Container& coarseNeighbors,
                                                         const bool ghostConnectivity = false) const = 0;

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

private:
  //-------------------------- Private Interface --------------------------//
  using Neighbor<Dimension>::mDataBase;
  using Neighbor<Dimension>::mNodeExtent;
  using Neighbor<Dimension>::mSearchType;
  using Neighbor<Dimension>::mKernelExtent;
};

}

#endif
