#include "Utilities/DBC.hh"

namespace Spheral {

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
SolidNodeListView<Dimension>::
SolidNodeListView(const size_t numNodes,
                  const size_t firstGhostNode,
                  const Scalar hmin,
                  const Scalar hmax,
                  const Scalar hminratio,
                  const Scalar nPerh,
                  const size_t maxNumNeighbors,
                  const Scalar rhoMin,
                  const Scalar rhoMax):
  NodeListView<Dimension>(numNodes, firstGhostNode, hmin, hmax, hminratio, nPerh, maxNumNeighbors),
  FluidNodeListViewType(numNodes, firstGhostNode, hmin, hmax, hminratio, nPerh, maxNumNeighbors, rhoMin, rhoMax) {
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
SolidNodeListView<Dimension>::
SolidNodeListView(const size_t numNodes,
                  const size_t firstGhostNode,
                  const SymTensorView& deviatoricStress,
                  const ScalarView& plasticStrain,
                  const ScalarView& plasticStrainRate,
                  const SymTensorView& damage,
                  const IntView& fragmentIDs,
                  const IntView& particleTypes):
  NodeListView<Dimension>(numNodes, firstGhostNode),
  FluidNodeListViewType(numNodes, firstGhostNode),
  mDeviatoricStressView(deviatoricStress),
  mPlasticStrainView(plasticStrain),
  mPlasticStrainRateView(plasticStrainRate),
  mDamageView(damage),
  mFragmentIDsView(fragmentIDs),
  mParticleTypesView(particleTypes) {
  REQUIRE(mDeviatoricStressView.numElements() == this->mNumNodes);
  REQUIRE(mPlasticStrainView.numElements() == this->mNumNodes);
  REQUIRE(mPlasticStrainRateView.numElements() == this->mNumNodes);
  REQUIRE(mDamageView.numElements() == this->mNumNodes);
  REQUIRE(mFragmentIDsView.numElements() == this->mNumNodes);
  REQUIRE(mParticleTypesView.numElements() == this->mNumNodes);
}

#define SPHERAL_SOLID_VIEW_ELEMENT(METHOD, MEMBER, TYPE)                     \
  template<typename Dimension>                                               \
  SPHERAL_HOST_DEVICE inline TYPE& SolidNodeListView<Dimension>::METHOD(     \
      const size_t nodeID) const {                                            \
    REQUIRE(nodeID < MEMBER.numElements());                                   \
    return MEMBER[nodeID];                                                    \
  }

SPHERAL_SOLID_VIEW_ELEMENT(deviatoricStress, mDeviatoricStressView, typename Dimension::SymTensor)
SPHERAL_SOLID_VIEW_ELEMENT(plasticStrain, mPlasticStrainView, typename Dimension::Scalar)
SPHERAL_SOLID_VIEW_ELEMENT(plasticStrainRate, mPlasticStrainRateView, typename Dimension::Scalar)
SPHERAL_SOLID_VIEW_ELEMENT(damage, mDamageView, typename Dimension::SymTensor)
SPHERAL_SOLID_VIEW_ELEMENT(fragmentID, mFragmentIDsView, int)
SPHERAL_SOLID_VIEW_ELEMENT(particleType, mParticleTypesView, int)

#undef SPHERAL_SOLID_VIEW_ELEMENT

template<typename Dimension>
SPHERAL_HOST
inline
void
SolidNodeListView<Dimension>::move(chai::ExecutionSpace space) {
  FluidNodeListViewType::move(space);
#ifndef SPHERAL_UNIFIED_MEMORY
  mDeviatoricStressView.move(space);
  mPlasticStrainView.move(space);
  mPlasticStrainRateView.move(space);
  mDamageView.move(space);
  mFragmentIDsView.move(space);
  mParticleTypesView.move(space);
#endif
}

template<typename Dimension>
SPHERAL_HOST
inline
void
SolidNodeListView<Dimension>::touch(chai::ExecutionSpace space) {
  FluidNodeListViewType::touch(space);
#ifndef SPHERAL_UNIFIED_MEMORY
  mDeviatoricStressView.touch(space);
  mPlasticStrainView.touch(space);
  mPlasticStrainRateView.touch(space);
  mDamageView.touch(space);
  mFragmentIDsView.touch(space);
  mParticleTypesView.touch(space);
#endif
}

}
