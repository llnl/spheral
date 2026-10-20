//---------------------------------Spheral++----------------------------------//
// SolidNodeListView -- allocation-free view/base state for SolidNodeList.
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidNodeListView_hh__
#define __Spheral_SolidNodeListView_hh__

#include "NodeList/FluidNodeListView.hh"

namespace Spheral {

template<typename Dimension>
class SolidNodeListView: public virtual FluidNodeListView<Dimension> {
public:
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using FluidNodeListViewType = FluidNodeListView<Dimension>;
  using ScalarView = FieldView<Dimension, Scalar>;
  using SymTensorView = FieldView<Dimension, SymTensor>;
  using IntView = FieldView<Dimension, int>;
  using DeviatoricStressView = SymTensorView;
  using PlasticStrainView = ScalarView;
  using DamageView = SymTensorView;
  using FragmentIDsView = IntView;
  using ParticleTypesView = IntView;

  SPHERAL_HOST_DEVICE SolidNodeListView() = default;
  SPHERAL_HOST_DEVICE SolidNodeListView(const size_t numNodes,
                                        const size_t firstGhostNode,
                                        const Scalar hmin = 1.0e-20,
                                        const Scalar hmax = 1.0e20,
                                        const Scalar hminratio = 0.1,
                                        const Scalar nPerh = 2.01,
                                        const size_t maxNumNeighbors = 500,
                                        const Scalar rhoMin = 1.0e-10,
                                        const Scalar rhoMax = 1.0e100);
  SPHERAL_HOST_DEVICE SolidNodeListView(const size_t numNodes,
                                        const size_t firstGhostNode,
                                        const SymTensorView& deviatoricStress,
                                        const ScalarView& plasticStrain,
                                        const ScalarView& plasticStrainRate,
                                        const SymTensorView& damage,
                                        const IntView& fragmentIDs,
                                        const IntView& particleTypes);
  SPHERAL_HOST_DEVICE SolidNodeListView(const SolidNodeListView&) = default;
  SPHERAL_HOST_DEVICE SolidNodeListView(SolidNodeListView&&) = default;
  SPHERAL_HOST_DEVICE virtual ~SolidNodeListView() = default;

  SPHERAL_HOST_DEVICE SolidNodeListView& operator=(SolidNodeListView&) = default;
  SPHERAL_HOST_DEVICE SolidNodeListView& operator=(SolidNodeListView&&) = default;

  SPHERAL_HOST_DEVICE SymTensorView& deviatoricStress()                       { return mDeviatoricStressView; }
  SPHERAL_HOST_DEVICE const SymTensorView& deviatoricStress()           const { return mDeviatoricStressView; }
  SPHERAL_HOST_DEVICE ScalarView& plasticStrain()                            { return mPlasticStrainView; }
  SPHERAL_HOST_DEVICE const ScalarView& plasticStrain()                const { return mPlasticStrainView; }
  SPHERAL_HOST_DEVICE ScalarView& plasticStrainRate()                        { return mPlasticStrainRateView; }
  SPHERAL_HOST_DEVICE const ScalarView& plasticStrainRate()            const { return mPlasticStrainRateView; }
  SPHERAL_HOST_DEVICE SymTensorView& damage()                                { return mDamageView; }
  SPHERAL_HOST_DEVICE const SymTensorView& damage()                    const { return mDamageView; }
  SPHERAL_HOST_DEVICE IntView& fragmentIDs()                                 { return mFragmentIDsView; }
  SPHERAL_HOST_DEVICE const IntView& fragmentIDs()                     const { return mFragmentIDsView; }
  SPHERAL_HOST_DEVICE IntView& particleTypes()                               { return mParticleTypesView; }
  SPHERAL_HOST_DEVICE const IntView& particleTypes()                   const { return mParticleTypesView; }

  SPHERAL_HOST_DEVICE SymTensor& deviatoricStress(size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& plasticStrain(size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& plasticStrainRate(size_t nodeID) const;
  SPHERAL_HOST_DEVICE SymTensor& damage(size_t nodeID) const;
  SPHERAL_HOST_DEVICE int& fragmentID(size_t nodeID) const;
  SPHERAL_HOST_DEVICE int& particleType(size_t nodeID) const;

  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  SymTensorView mDeviatoricStressView;
  ScalarView mPlasticStrainView;
  ScalarView mPlasticStrainRateView;
  SymTensorView mDamageView;
  IntView mFragmentIDsView;
  IntView mParticleTypesView;
};

}

#include "SolidNodeListViewInline.hh"

#endif
