//---------------------------------Spheral++----------------------------------//
// NodeListView -- allocation-free view/base state for NodeList.
//
// This intentionally avoids owning containers such as std::string and std::vector.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeListView_hh__
#define __Spheral_NodeListView_hh__

#include "config.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

enum class NodeType {
  InternalNode = 0,
  GhostNode = 1
};

}

#include "Field/FieldView.hh"

namespace Spheral {

template<typename Dimension>
class NodeListView {// : public chai::CHAICopyable {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;

  using ScalarView = FieldView<Dimension, Scalar>;
  using VectorView = FieldView<Dimension, Vector>;
  using SymTensorView = FieldView<Dimension, SymTensor>;

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
                                   const VectorView& positions,
                                   const SymTensorView& Hfield,
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

  // Comparison
  SPHERAL_HOST_DEVICE bool operator==(const NodeListView&) const = default;

  // Node count accessors.
  SPHERAL_HOST_DEVICE size_t numNodes()                         const { return mNumNodes; }
  SPHERAL_HOST_DEVICE size_t numInternalNodes()                 const { return mFirstGhostNode; }
  SPHERAL_HOST_DEVICE size_t numGhostNodes()                    const { CHECK(mFirstGhostNode <= mNumNodes); return mNumNodes - mFirstGhostNode; }
  SPHERAL_HOST_DEVICE size_t firstGhostNode()                   const { return mFirstGhostNode; }
  SPHERAL_HOST_DEVICE NodeType nodeType(size_t i)               const { REQUIRE(i < mNumNodes); return (i < mFirstGhostNode ? NodeType::InternalNode : NodeType::GhostNode); }

  // Smoothing scale controls.
  SPHERAL_HOST_DEVICE Scalar nodesPerSmoothingScale()           const { return mNodesPerSmoothingScale; }
  SPHERAL_HOST_DEVICE void nodesPerSmoothingScale(Scalar val)         { mNodesPerSmoothingScale = val; }
  SPHERAL_HOST_DEVICE size_t maxNumNeighbors()                  const { return mMaxNumNeighbors; }
  SPHERAL_HOST_DEVICE void maxNumNeighbors(size_t val)                { mMaxNumNeighbors = val; }
  SPHERAL_HOST_DEVICE Scalar hmin()                             const { return mhmin; }
  SPHERAL_HOST_DEVICE void hmin(Scalar val)                           { mhmin = val; }
  SPHERAL_HOST_DEVICE Scalar hmax()                             const { return mhmax; }
  SPHERAL_HOST_DEVICE void hmax(Scalar val)                           { mhmax = val; }
  SPHERAL_HOST_DEVICE Scalar hminratio()                        const { return mhminratio; }
  SPHERAL_HOST_DEVICE void hminratio(Scalar val)                      { mhminratio = val; }

  // Field views.
  SPHERAL_HOST_DEVICE ScalarView& mass()                              { return mMassView; }
  SPHERAL_HOST_DEVICE const ScalarView& mass()                  const { return mMassView; }
  SPHERAL_HOST_DEVICE ScalarView& massRZ()                            { return mMassRZView; }
  SPHERAL_HOST_DEVICE const ScalarView& massRZ()                const { return mMassRZView; }
  SPHERAL_HOST_DEVICE VectorView& positions()                         { return mPositionsView; }
  SPHERAL_HOST_DEVICE const VectorView& positions()             const { return mPositionsView; }
  SPHERAL_HOST_DEVICE VectorView& velocity()                          { return mVelocityView; }
  SPHERAL_HOST_DEVICE const VectorView& velocity()              const { return mVelocityView; }
  SPHERAL_HOST_DEVICE SymTensorView& Hfield()                         { return mHfieldView; }
  SPHERAL_HOST_DEVICE const SymTensorView& Hfield()             const { return mHfieldView; }
  SPHERAL_HOST_DEVICE ScalarView& work()                              { return mWorkView; }
  SPHERAL_HOST_DEVICE const ScalarView& work()                  const { return mWorkView; }

  SPHERAL_HOST_DEVICE Scalar& mass(const size_t nodeID)         const { REQUIRE(nodeID < mMassView.numElements()); return mMassView[nodeID]; }
  SPHERAL_HOST_DEVICE Scalar& massRZ(const size_t nodeID)       const { REQUIRE(nodeID < mMassRZView.numElements()); return mMassRZView[nodeID]; }
  SPHERAL_HOST_DEVICE Vector& position(const size_t nodeID)     const { REQUIRE(nodeID < mPositionsView.numElements()); return mPositionsView[nodeID]; }
  SPHERAL_HOST_DEVICE Vector& velocity(const size_t nodeID)     const { REQUIRE(nodeID < mVelocityView.numElements()); return mVelocityView[nodeID]; }
  SPHERAL_HOST_DEVICE SymTensor& H(const size_t nodeID)         const { REQUIRE(nodeID < mHfieldView.numElements()); return mHfieldView[nodeID]; }
  SPHERAL_HOST_DEVICE Scalar& work(const size_t nodeID)         const { REQUIRE(nodeID < mWorkView.numElements()); return mWorkView[nodeID]; }

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
  ScalarView mMassView;
  ScalarView mMassRZView;
  VectorView mPositionsView;
  VectorView mVelocityView;
  SymTensorView mHfieldView;
  ScalarView mWorkView;
};

}

#include "NodeListViewInline.hh"

#endif
