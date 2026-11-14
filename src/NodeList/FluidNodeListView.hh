//---------------------------------Spheral++----------------------------------//
// FluidNodeListView -- allocation-free view/base state for FluidNodeList.
//----------------------------------------------------------------------------//
#ifndef __Spheral_FluidNodeListView_hh__
#define __Spheral_FluidNodeListView_hh__

#include "NodeList/NodeListView.hh"

namespace Spheral {

template<typename Dimension>
class FluidNodeListView: public NodeListView<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using NodeListViewType = NodeListView<Dimension>;
  using ScalarFieldView = FieldView<Dimension, Scalar>;

  // Constructors, destructor.
  SPHERAL_HOST_DEVICE FluidNodeListView() = default;
  SPHERAL_HOST_DEVICE FluidNodeListView(const size_t numNodes,
                                        const size_t firstGhostNode,
                                        const Scalar hmin = 1.0e-20,
                                        const Scalar hmax = 1.0e20,
                                        const Scalar hminratio = 0.1,
                                        const Scalar nPerh = 2.01,
                                        const size_t maxNumNeighbors = 500,
                                        const Scalar rhoMin = 1.0e-10,
                                        const Scalar rhoMax = 1.0e100);
  SPHERAL_HOST_DEVICE FluidNodeListView(const size_t numNodes,
                                        const size_t firstGhostNode,
                                        const ScalarFieldView& massDensity,
                                        const ScalarFieldView& specificThermalEnergy,
                                        const Scalar rhoMin = 1.0e-10,
                                        const Scalar rhoMax = 1.0e100);
  SPHERAL_HOST_DEVICE FluidNodeListView(const NodeListViewType& nodeListView,
                                        const ScalarFieldView& massDensity,
                                        const ScalarFieldView& massDensityRZ,
                                        const ScalarFieldView& specificThermalEnergy,
                                        const Scalar rhoMin = 1.0e-10,
                                        const Scalar rhoMax = 1.0e100);
  SPHERAL_HOST_DEVICE FluidNodeListView(const FluidNodeListView&) = default;
  SPHERAL_HOST_DEVICE FluidNodeListView(FluidNodeListView&&) = default;
  SPHERAL_HOST_DEVICE ~FluidNodeListView() = default;

  // Assignment.
  SPHERAL_HOST_DEVICE FluidNodeListView& operator=(FluidNodeListView&) = default;
  SPHERAL_HOST_DEVICE FluidNodeListView& operator=(FluidNodeListView&&) = default;

  // Comparison
  SPHERAL_HOST_DEVICE bool operator==(const FluidNodeListView& rhs) const {
    return NodeListViewType::operator==(rhs);
  }

  // Field views.
  SPHERAL_HOST_DEVICE ScalarFieldView& massDensity()                              { return mMassDensityView; }
  SPHERAL_HOST_DEVICE const ScalarFieldView& massDensity()                  const { return mMassDensityView; }
  SPHERAL_HOST_DEVICE ScalarFieldView& massDensityRZ()                            { return mMassDensityRZView; }
  SPHERAL_HOST_DEVICE const ScalarFieldView& massDensityRZ()                const { return mMassDensityRZView; }
  SPHERAL_HOST_DEVICE ScalarFieldView& specificThermalEnergy()                    { return mSpecificThermalEnergyView; }
  SPHERAL_HOST_DEVICE const ScalarFieldView& specificThermalEnergy()        const { return mSpecificThermalEnergyView; }

  // Index for values
  SPHERAL_HOST_DEVICE Scalar& massDensity(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& massDensityRZ(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& specificThermalEnergy(const size_t nodeID) const;

  // Properties
  SPHERAL_HOST_DEVICE Scalar rhoMin()                              const { return mRhoMin; }
  SPHERAL_HOST_DEVICE Scalar rhoMax()                              const { return mRhoMax; }
  SPHERAL_HOST_DEVICE void rhoMin(const Scalar x)                        { mRhoMin = x; }
  SPHERAL_HOST_DEVICE void rhoMax(const Scalar x)                        { mRhoMax = x; }

  // CHAI/ManagedArray specific operations for the contained FieldViews.
  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  //--------------------------- Protected Interface ---------------------------//
  Scalar mRhoMin = 1.0e-10;
  Scalar mRhoMax = 1.0e100;
  ScalarFieldView mMassDensityView;
  ScalarFieldView mMassDensityRZView;
  ScalarFieldView mSpecificThermalEnergyView;
};

}

#include "FluidNodeListViewInline.hh"

#endif
