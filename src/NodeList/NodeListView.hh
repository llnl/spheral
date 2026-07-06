//---------------------------------Spheral++----------------------------------//
// NodeListView -- allocation-free view/base state for NodeList.
//
// This intentionally avoids owning containers such as std::string and std::vector.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeListView_hh__
#define __Spheral_NodeListView_hh__

#include "config.hh"

namespace Spheral {

enum class NodeType {
  InternalNode = 0,
  GhostNode = 1
};

}

#include "Field/FieldView.hh"

namespace Spheral {

template<typename Dimension>
class NodeListView: public chai::CHAICopyable {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;

  using MassView = FieldView<Dimension, Scalar>;
  using PositionView = FieldView<Dimension, Vector>;
  using VelocityView = FieldView<Dimension, Vector>;
  using HfieldView = FieldView<Dimension, SymTensor>;
  using WorkView = FieldView<Dimension, Scalar>;

  // Constructors, destructor.
  SPHERAL_HOST_DEVICE NodeListView() = default;
  SPHERAL_HOST_DEVICE NodeListView(const size_t numNodes,
                                   const size_t firstGhostNode,
                                   const Scalar hmin = 1.0e-20,
                                   const Scalar hmax = 1.0e20,
                                   const Scalar hminratio = 0.1,
                                   const Scalar nPerh = 2.01,
                                   const size_t maxNumNeighbors = 500);
  SPHERAL_HOST_DEVICE NodeListView(const size_t numNodes,
                                   const size_t firstGhostNode,
                                   const PositionView& positions,
                                   const HfieldView& Hfield,
                                   const Scalar hmin = 1.0e-20,
                                   const Scalar hmax = 1.0e20,
                                   const Scalar hminratio = 0.1,
                                   const Scalar nPerh = 2.01,
                                   const size_t maxNumNeighbors = 500);
  SPHERAL_HOST_DEVICE NodeListView(const size_t numNodes,
                                   const size_t firstGhostNode,
                                   const MassView& mass,
                                   const PositionView& positions,
                                   const VelocityView& velocity,
                                   const HfieldView& Hfield,
                                   const WorkView& work,
                                   const Scalar hmin = 1.0e-20,
                                   const Scalar hmax = 1.0e20,
                                   const Scalar hminratio = 0.1,
                                   const Scalar nPerh = 2.01,
                                   const size_t maxNumNeighbors = 500);
  SPHERAL_HOST_DEVICE NodeListView(const NodeListView& rhs) = default;
  SPHERAL_HOST_DEVICE NodeListView(NodeListView&& rhs) = default;
  SPHERAL_HOST_DEVICE virtual ~NodeListView() = default;

  // Assignment.
  SPHERAL_HOST_DEVICE NodeListView& operator=(NodeListView& rhs) = default;
  SPHERAL_HOST_DEVICE NodeListView& operator=(NodeListView&& rhs) = default;

  // Node count accessors.
  SPHERAL_HOST_DEVICE size_t numNodes() const;
  SPHERAL_HOST_DEVICE size_t numInternalNodes() const;
  SPHERAL_HOST_DEVICE size_t numGhostNodes() const;
  SPHERAL_HOST_DEVICE size_t firstGhostNode() const;
  SPHERAL_HOST_DEVICE NodeType nodeType(size_t i) const;

  // Smoothing scale controls.
  SPHERAL_HOST_DEVICE Scalar nodesPerSmoothingScale() const;
  SPHERAL_HOST_DEVICE void nodesPerSmoothingScale(Scalar val);
  SPHERAL_HOST_DEVICE size_t maxNumNeighbors() const;
  SPHERAL_HOST_DEVICE void maxNumNeighbors(size_t val);
  SPHERAL_HOST_DEVICE Scalar hmin() const;
  SPHERAL_HOST_DEVICE void hmin(Scalar val);
  SPHERAL_HOST_DEVICE Scalar hmax() const;
  SPHERAL_HOST_DEVICE void hmax(Scalar val);
  SPHERAL_HOST_DEVICE Scalar hminratio() const;
  SPHERAL_HOST_DEVICE void hminratio(Scalar val);

  // Field views.
  SPHERAL_HOST_DEVICE MassView& mass();
  SPHERAL_HOST_DEVICE const MassView& mass() const;
  SPHERAL_HOST_DEVICE PositionView& positions();
  SPHERAL_HOST_DEVICE const PositionView& positions() const;
  SPHERAL_HOST_DEVICE VelocityView& velocity();
  SPHERAL_HOST_DEVICE const VelocityView& velocity() const;
  SPHERAL_HOST_DEVICE HfieldView& Hfield();
  SPHERAL_HOST_DEVICE const HfieldView& Hfield() const;
  SPHERAL_HOST_DEVICE WorkView& work();
  SPHERAL_HOST_DEVICE const WorkView& work() const;

  SPHERAL_HOST_DEVICE Scalar& mass(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Vector& position(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Vector& velocity(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE SymTensor& H(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& work(const size_t nodeID) const;

  // CHAI/ManagedArray specific operations for the contained FieldViews.
  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  //--------------------------- Protected Interface ---------------------------//
  size_t mNumNodes = 0u;
  size_t mFirstGhostNode = 0u;
  Scalar mhmin = 1.0e-20;
  Scalar mhmax = 1.0e20;
  Scalar mhminratio = 0.1;
  Scalar mNodesPerSmoothingScale = 2.01;
  size_t mMaxNumNeighbors = 500u;
  MassView mMassView;
  PositionView mPositionsView;
  VelocityView mVelocityView;
  HfieldView mHfieldView;
  WorkView mWorkView;
};

}

#include "NodeListViewInline.hh"

#endif
