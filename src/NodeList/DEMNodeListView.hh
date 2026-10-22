//---------------------------------Spheral++----------------------------------//
// DEMNodeListView -- allocation-free view/base state for DEMNodeList.
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMNodeListView_hh__
#define __Spheral_DEMNodeListView_hh__

#include "NodeList/NodeListView.hh"

namespace Spheral {

template<typename Dimension>
class DEMNodeListView: public virtual NodeListView<Dimension> {
public:
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;

  using NodeListViewType = NodeListView<Dimension>;

  using ScalarView = FieldView<Dimension, Scalar>;
  using IntView = FieldView<Dimension, int>;
  using IndexView = FieldView<Dimension, size_t>;

  // Constructors, destructor.
  SPHERAL_HOST_DEVICE DEMNodeListView() = default;
  SPHERAL_HOST_DEVICE DEMNodeListView(const size_t numNodes,
                                      const size_t firstGhostNode,
                                      const Scalar hmin = 1.0e-20,
                                      const Scalar hmax = 1.0e20,
                                      const Scalar hminratio = 0.1,
                                      const Scalar nPerh = 2.01,
                                      const Scalar neighborSearchBuffer = 0.1,
                                      const size_t maxNumNeighbors = 500);
  SPHERAL_HOST_DEVICE DEMNodeListView(const size_t numNodes,
                                      const size_t firstGhostNode,
                                      const ScalarView& particleRadius,
                                      const IntView& compositeParticleIndex,
                                      const IndexView& uniqueIndex,
                                      const Scalar neighborSearchBuffer = 0.1);
  SPHERAL_HOST_DEVICE DEMNodeListView(const DEMNodeListView&) = default;
  SPHERAL_HOST_DEVICE DEMNodeListView(DEMNodeListView&&) = default;
  SPHERAL_HOST_DEVICE virtual ~DEMNodeListView() = default;

  // Assignment.
  SPHERAL_HOST_DEVICE DEMNodeListView& operator=(DEMNodeListView&) = default;
  SPHERAL_HOST_DEVICE DEMNodeListView& operator=(DEMNodeListView&&) = default;

  // Comparison
  SPHERAL_HOST_DEVICE bool operator==(const DEMNodeListView&) const = default;

  // Field views.
  SPHERAL_HOST_DEVICE ScalarView& particleRadius()                          { return mParticleRadiusView; }
  SPHERAL_HOST_DEVICE const ScalarView& particleRadius()              const { return mParticleRadiusView; }
  SPHERAL_HOST_DEVICE IntView& compositeParticleIndex()                     { return mCompositeParticleIndexView; }
  SPHERAL_HOST_DEVICE const IntView& compositeParticleIndex()         const { return mCompositeParticleIndexView; }
  SPHERAL_HOST_DEVICE IndexView& uniqueIndex()                              { return mUniqueIndexView; }
  SPHERAL_HOST_DEVICE const IndexView& uniqueIndex()                  const { return mUniqueIndexView; }

  // DEM particle properties
  SPHERAL_HOST_DEVICE Scalar& particleRadius(size_t nodeID) const;
  SPHERAL_HOST_DEVICE int& compositeParticleIndex(size_t nodeID) const;
  SPHERAL_HOST_DEVICE size_t& uniqueIndex(size_t nodeID) const;

  SPHERAL_HOST_DEVICE Scalar neighborSearchBuffer()                   const { return mNeighborSearchBuffer; }
  SPHERAL_HOST_DEVICE void neighborSearchBuffer(const Scalar x)             { mNeighborSearchBuffer = x; }

  SPHERAL_HOST_DEVICE void setHfieldFromParticleRadius(const size_t startUniqueIndex,
                                                       const Scalar kernelExtent);

  // CHAI/ManagedArray specific operations for the contained FieldViews.
  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  Scalar mNeighborSearchBuffer = 0.1;
  ScalarView mParticleRadiusView;
  IntView mCompositeParticleIndexView;
  IndexView mUniqueIndexView;
};

}

#include "DEMNodeListViewInline.hh"

#endif
