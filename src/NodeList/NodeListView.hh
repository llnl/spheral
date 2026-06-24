//---------------------------------Spheral++----------------------------------//
// NodeListView -- CPU/GPU view of NodeList data needed by ConnectivityMap.
//
// This intentionally exposes only NodeList-owned data used in neighbor finding.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeListView_hh__
#define __Spheral_NodeListView_hh__

#include "config.hh"

#include "Field/FieldView.hh"

namespace Spheral {

template<typename Dimension>
class NodeListView: public chai::CHAICopyable {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;

  using PositionView = FieldView<Dimension, Vector>;
  using HfieldView = FieldView<Dimension, SymTensor>;

  // Constructors, destructor.
  SPHERAL_HOST_DEVICE NodeListView() = default;
  SPHERAL_HOST_DEVICE NodeListView(const size_t numNodes,
                                   const size_t firstGhostNode);
  SPHERAL_HOST_DEVICE NodeListView(const size_t numNodes,
                                   const size_t firstGhostNode,
                                   const PositionView& positions,
                                   const HfieldView& Hfield);
  SPHERAL_HOST_DEVICE NodeListView(const NodeListView& rhs) = default;
  SPHERAL_HOST_DEVICE NodeListView(NodeListView&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~NodeListView() = default;

  // Assignment.
  SPHERAL_HOST_DEVICE NodeListView& operator=(NodeListView& rhs) = default;
  SPHERAL_HOST_DEVICE NodeListView& operator=(NodeListView&& rhs) = default;

  // Node count accessors.
  SPHERAL_HOST_DEVICE size_t numNodes() const;
  SPHERAL_HOST_DEVICE size_t numInternalNodes() const;
  SPHERAL_HOST_DEVICE size_t firstGhostNode() const;

  // ConnectivityMap-facing Field views.
  SPHERAL_HOST_DEVICE PositionView& positions();
  SPHERAL_HOST_DEVICE const PositionView& positions() const;
  SPHERAL_HOST_DEVICE HfieldView& Hfield();
  SPHERAL_HOST_DEVICE const HfieldView& Hfield() const;

  SPHERAL_HOST_DEVICE Vector& position(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE SymTensor& H(const size_t nodeID) const;

  // CHAI/ManagedArray specific operations for the contained FieldViews.
  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  //--------------------------- Protected Interface ---------------------------//
  size_t mNumNodes = 0u;
  size_t mFirstGhostNode = 0u;
  PositionView mPositionsView;
  HfieldView mHfieldView;
};

}

#include "NodeListViewInline.hh"

#endif
