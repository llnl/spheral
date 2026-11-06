#include "Utilities/DBC.hh"

namespace Spheral {

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
FluidNodeListView<Dimension>::
FluidNodeListView(const size_t numNodes,
                  const size_t firstGhostNode,
                  const Scalar hmin,
                  const Scalar hmax,
                  const Scalar hminratio,
                  const Scalar nPerh,
                  const size_t maxNumNeighbors,
                  const Scalar rhoMin,
                  const Scalar rhoMax):
  NodeListViewType(numNodes, firstGhostNode, hmin, hmax, hminratio, nPerh, maxNumNeighbors),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax),
  mMassDensityView(),
  mMassDensityRZView(),
  mSpecificThermalEnergyView() {
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
FluidNodeListView<Dimension>::
FluidNodeListView(const size_t numNodes,
                  const size_t firstGhostNode,
                  const ScalarView& massDensity,
                  const ScalarView& specificThermalEnergy,
                  const Scalar rhoMin,
                  const Scalar rhoMax):
  NodeListViewType(numNodes, firstGhostNode),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax),
  mMassDensityView(massDensity),
  mMassDensityRZView(),
  mSpecificThermalEnergyView(specificThermalEnergy) {
  REQUIRE(mMassDensityView.numElements() == this->mNumNodes);
  REQUIRE(mSpecificThermalEnergyView.numElements() == this->mNumNodes);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename FluidNodeListView<Dimension>::Scalar&
FluidNodeListView<Dimension>::
massDensity(const size_t nodeID) const {
  REQUIRE(nodeID < mMassDensityView.numElements());
  return mMassDensityView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename FluidNodeListView<Dimension>::Scalar&
FluidNodeListView<Dimension>::
massDensityRZ(const size_t nodeID) const {
  REQUIRE(nodeID < mMassDensityRZView.numElements());
  return mMassDensityRZView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename FluidNodeListView<Dimension>::Scalar&
FluidNodeListView<Dimension>::
specificThermalEnergy(const size_t nodeID) const {
  REQUIRE(nodeID < mSpecificThermalEnergyView.numElements());
  return mSpecificThermalEnergyView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST
inline
void
FluidNodeListView<Dimension>::
move(chai::ExecutionSpace space) {
  NodeListViewType::move(space);
#ifndef SPHERAL_UNIFIED_MEMORY
  mMassDensityView.move(space);
  mMassDensityRZView.move(space);
  mSpecificThermalEnergyView.move(space);
#endif
}

template<typename Dimension>
SPHERAL_HOST
inline
void
FluidNodeListView<Dimension>::
touch(chai::ExecutionSpace space) {
  NodeListViewType::touch(space);
#ifndef SPHERAL_UNIFIED_MEMORY
  mMassDensityView.touch(space);
  mMassDensityRZView.touch(space);
  mSpecificThermalEnergyView.touch(space);
#endif
}

}
