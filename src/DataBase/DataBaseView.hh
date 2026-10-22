//---------------------------------Spheral++----------------------------------//
// DataBaseView -- A device capable view of DataBase
//
// Created by JMO, Tue Oct 20 13:45:31 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_DataBaseView__
#define __Spheral_DataBaseView__

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

#include "Geometry/Dimension.hh"
// #include "NodeList/NodeListView.hh"
// #include "NodeList/FluidNodeListView.hh"
// #include "NodeList/SolidNodeListView.hh"
// #include "NodeList/DEMNodeListView.hh"
#include "Utilities/span.hh"

namespace Spheral {

// Forward declarations
template<typename Dimension, typename DataType> class FieldView;
template<typename Dimension, typename DataType> class FieldListView;
template<typename Dimension> class NodeListView;
template<typename Dimension> class FluidNodeListView;
template<typename Dimension> class SolidNodeListView;
template<typename Dimension> class DEMNodeListView;

template<typename Dimension>
class DataBaseView: public chai::CHAICopyable {
public:
  //--------------------------- Public Interface ---------------------------//
  using DataBaseDimension = Dimension;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

#ifdef SPHERAL_UNIFIED_MEMORY
  template<typename T> using Container = SPHERAL_SPAN_TYPE<T>;
#else
  template<typename T> using Container = typename chai::ManagedArray<T>;
#endif

  // It is convenient to be able to query the DataBase for the problem
  // dimensionality for Python.
  static int nDim;

  // Constructors.
  SPHERAL_HOST_DEVICE DataBaseView() = default;
  SPHERAL_HOST_DEVICE DataBaseView(const DataBaseView& rhs) = default;
  SPHERAL_HOST_DEVICE DataBaseView(DataBaseView&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~DataBaseView() = default;

  // Assignment.
  SPHERAL_HOST_DEVICE DataBaseView& operator=(DataBaseView& rhs) = default;
  SPHERAL_HOST_DEVICE DataBaseView& operator=(DataBaseView&& rhs) = default;

  // Number of NodeLists we have in the DataBase.
  SPHERAL_HOST_DEVICE size_t numNodeLists()                                               const { return mNodeListViews.size(); }
  SPHERAL_HOST_DEVICE size_t numFluidNodeLists()                                          const { return mFluidNodeListViews.size(); }
  SPHERAL_HOST_DEVICE size_t numSolidNodeLists()                                          const { return mSolidNodeListViews.size(); }
  SPHERAL_HOST_DEVICE size_t numDEMNodeLists()                                            const { return mDEMNodeListViews.size(); }

  // Numbers of nodes.
  SPHERAL_HOST_DEVICE size_t numInternalNodes() const;
  SPHERAL_HOST_DEVICE size_t numGhostNodes() const;
  SPHERAL_HOST_DEVICE size_t numNodes() const;

  // Numbers of fluid nodes.
  SPHERAL_HOST_DEVICE size_t numFluidInternalNodes() const;
  SPHERAL_HOST_DEVICE size_t numFluidGhostNodes() const;
  SPHERAL_HOST_DEVICE size_t numFluidNodes() const;

  // Numbers of solid nodes.
  SPHERAL_HOST_DEVICE size_t numSolidInternalNodes() const;
  SPHERAL_HOST_DEVICE size_t numSolidGhostNodes() const;
  SPHERAL_HOST_DEVICE size_t numSolidNodes() const;

  // Numbers of DEM nodes.
  SPHERAL_HOST_DEVICE size_t numDEMInternalNodes() const;
  SPHERAL_HOST_DEVICE size_t numDEMGhostNodes() const;
  SPHERAL_HOST_DEVICE size_t numDEMNodes() const;

  // Access the various NodeLists
  SPHERAL_HOST_DEVICE const Container<NodeListView<Dimension>>& nodeListViews()           const { return mNodeListViews; }
  SPHERAL_HOST_DEVICE const Container<FluidNodeListView<Dimension>>& fluidNodeListViews() const { return mFluidNodeListViews; }
  SPHERAL_HOST_DEVICE const Container<SolidNodeListView<Dimension>>& solidNodeListViews() const { return mSolidNodeListViews; }
  SPHERAL_HOST_DEVICE const Container<DEMNodeListView<Dimension>>& DEMNodeListViews()     const { return mDEMNodeListViews; }

  // Query about NodeLists
  SPHERAL_HOST_DEVICE bool haveNodeList(const NodeListView<Dimension>& x) const;
  SPHERAL_HOST_DEVICE bool haveNodeList(const FluidNodeListView<Dimension>& x) const;
  SPHERAL_HOST_DEVICE bool haveNodeList(const SolidNodeListView<Dimension>& x) const;
  SPHERAL_HOST_DEVICE bool haveNodeList(const DEMNodeListView<Dimension>& x) const;

  SPHERAL_HOST_DEVICE size_t nodeListIndex(const NodeListView<Dimension>& x) const;
  SPHERAL_HOST_DEVICE size_t nodeListIndex(const FluidNodeListView<Dimension>& x) const;
  SPHERAL_HOST_DEVICE size_t nodeListIndex(const SolidNodeListView<Dimension>& x) const;
  SPHERAL_HOST_DEVICE size_t nodeListIndex(const DEMNodeListView<Dimension>& x) const;

  // Special DEM methods
  SPHERAL_HOST_DEVICE void setDEMHfieldFromParticleRadius(const size_t startUniqueIndex,
                                                          const Scalar kernelExtent);
  SPHERAL_HOST_DEVICE Scalar maxDEMNeighborSearchBuffer() const;

  // Get the maximum kernel extent from the Neighbor
  SPHERAL_HOST_DEVICE Scalar maxKernelExtent()                                            const { VERIFY(false); return 0.0; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  Container<NodeListView<Dimension>> mNodeListViews;
  Container<FluidNodeListView<Dimension>> mFluidNodeListViews;
  Container<SolidNodeListView<Dimension>> mSolidNodeListViews;
  Container<DEMNodeListView<Dimension>> mDEMNodeListViews;
};

//------------------------------------------------------------------------------
// Static varaible initialization
//------------------------------------------------------------------------------
template<typename Dimension> int DataBaseView<Dimension>::nDim = Dimension::nDim;

}

#endif
