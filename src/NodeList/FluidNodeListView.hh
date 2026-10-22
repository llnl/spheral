//---------------------------------Spheral++----------------------------------//
// FluidNodeListView -- allocation-free view/base state for FluidNodeList.
//----------------------------------------------------------------------------//
#ifndef __Spheral_FluidNodeListView_hh__
#define __Spheral_FluidNodeListView_hh__

#include "NodeList/NodeListView.hh"

namespace Spheral {

template<typename Dimension>
class FluidNodeListView: public virtual NodeListView<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using NodeListViewType = NodeListView<Dimension>;
  using ScalarView = FieldView<Dimension, Scalar>;

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
                                        const ScalarView& massDensity,
                                        const ScalarView& specificThermalEnergy,
                                        const Scalar rhoMin = 1.0e-10,
                                        const Scalar rhoMax = 1.0e100);
  SPHERAL_HOST_DEVICE FluidNodeListView(const FluidNodeListView&) = default;
  SPHERAL_HOST_DEVICE FluidNodeListView(FluidNodeListView&&) = default;
  SPHERAL_HOST_DEVICE virtual ~FluidNodeListView() = default;

  SPHERAL_HOST_DEVICE FluidNodeListView& operator=(FluidNodeListView&) = default;
  SPHERAL_HOST_DEVICE FluidNodeListView& operator=(FluidNodeListView&&) = default;

  SPHERAL_HOST_DEVICE bool operator==(const FluidNodeListView&) const = default;

  SPHERAL_HOST_DEVICE ScalarView& massDensity()                              { return mMassDensityView; }
  SPHERAL_HOST_DEVICE const ScalarView& massDensity()                  const { return mMassDensityView; }
  SPHERAL_HOST_DEVICE ScalarView& massDensityRZ()                            { return mMassDensityRZView; }
  SPHERAL_HOST_DEVICE const ScalarView& massDensityRZ()                const { return mMassDensityRZView; }
  SPHERAL_HOST_DEVICE ScalarView& specificThermalEnergy()                    { return mSpecificThermalEnergyView; }
  SPHERAL_HOST_DEVICE const ScalarView& specificThermalEnergy()        const { return mSpecificThermalEnergyView; }

  SPHERAL_HOST_DEVICE Scalar& massDensity(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& massDensityRZ(const size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& specificThermalEnergy(const size_t nodeID) const;

  SPHERAL_HOST_DEVICE Scalar rhoMin()                              const { return mRhoMin; }
  SPHERAL_HOST_DEVICE Scalar rhoMax()                              const { return mRhoMax; }
  SPHERAL_HOST_DEVICE void rhoMin(const Scalar x)                        { mRhoMin = x; }
  SPHERAL_HOST_DEVICE void rhoMax(const Scalar x)                        { mRhoMax = x; }

  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  //--------------------------- Protected Interface ---------------------------//
  Scalar mRhoMin = 1.0e-10;
  Scalar mRhoMax = 1.0e100;
  ScalarView mMassDensityView;
  ScalarView mMassDensityRZView;
  ScalarView mSpecificThermalEnergyView;
};

}

#include "FluidNodeListViewInline.hh"

#endif
