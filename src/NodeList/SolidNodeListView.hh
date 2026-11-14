//---------------------------------Spheral++----------------------------------//
// SolidNodeListView -- allocation-free view/base state for SolidNodeList.
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidNodeListView_hh__
#define __Spheral_SolidNodeListView_hh__

#include "NodeList/FluidNodeListView.hh"

namespace Spheral {

template<typename Dimension>
class SolidNodeListView: public FluidNodeListView<Dimension> {
public:
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using FluidNodeListViewType = FluidNodeListView<Dimension>;
  using ScalarFieldView = FieldView<Dimension, Scalar>;
  using SymTensorFieldView = FieldView<Dimension, SymTensor>;
  using IntFieldView = FieldView<Dimension, int>;

  // Constructors, destructor.
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
                                        const SymTensorFieldView& deviatoricStress,
                                        const ScalarFieldView& plasticStrain,
                                        const ScalarFieldView& plasticStrainRate,
                                        const SymTensorFieldView& damage,
                                        const IntFieldView& fragmentIDs,
                                        const IntFieldView& particleTypes);
  SPHERAL_HOST_DEVICE SolidNodeListView(const FluidNodeListViewType& fluidNodeListView,
                                        const SymTensorFieldView& deviatoricStress,
                                        const ScalarFieldView& plasticStrain,
                                        const ScalarFieldView& plasticStrainRate,
                                        const SymTensorFieldView& damage,
                                        const IntFieldView& fragmentIDs,
                                        const IntFieldView& particleTypes);
  SPHERAL_HOST_DEVICE SolidNodeListView(const SolidNodeListView&) = default;
  SPHERAL_HOST_DEVICE SolidNodeListView(SolidNodeListView&&) = default;
  SPHERAL_HOST_DEVICE ~SolidNodeListView() = default;

  // Assignment.
  SPHERAL_HOST_DEVICE SolidNodeListView& operator=(SolidNodeListView&) = default;
  SPHERAL_HOST_DEVICE SolidNodeListView& operator=(SolidNodeListView&&) = default;

  // Comparison
  SPHERAL_HOST_DEVICE bool operator==(const SolidNodeListView& rhs) const {
    return FluidNodeListViewType::operator==(rhs);
  }

  // Field views.
  SPHERAL_HOST_DEVICE SymTensorFieldView& deviatoricStress()                      { return mDeviatoricStressView; }
  SPHERAL_HOST_DEVICE const SymTensorFieldView& deviatoricStress()          const { return mDeviatoricStressView; }
  SPHERAL_HOST_DEVICE ScalarFieldView& plasticStrain()                            { return mPlasticStrainView; }
  SPHERAL_HOST_DEVICE const ScalarFieldView& plasticStrain()                const { return mPlasticStrainView; }
  SPHERAL_HOST_DEVICE ScalarFieldView& plasticStrainRate()                        { return mPlasticStrainRateView; }
  SPHERAL_HOST_DEVICE const ScalarFieldView& plasticStrainRate()            const { return mPlasticStrainRateView; }
  SPHERAL_HOST_DEVICE SymTensorFieldView& damage()                                { return mDamageView; }
  SPHERAL_HOST_DEVICE const SymTensorFieldView& damage()                    const { return mDamageView; }
  SPHERAL_HOST_DEVICE IntFieldView& fragmentIDs()                                 { return mFragmentIDsView; }
  SPHERAL_HOST_DEVICE const IntFieldView& fragmentIDs()                     const { return mFragmentIDsView; }
  SPHERAL_HOST_DEVICE IntFieldView& particleTypes()                               { return mParticleTypesView; }
  SPHERAL_HOST_DEVICE const IntFieldView& particleTypes()                   const { return mParticleTypesView; }

  // Index for values
  SPHERAL_HOST_DEVICE SymTensor& deviatoricStress(size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& plasticStrain(size_t nodeID) const;
  SPHERAL_HOST_DEVICE Scalar& plasticStrainRate(size_t nodeID) const;
  SPHERAL_HOST_DEVICE SymTensor& damage(size_t nodeID) const;
  SPHERAL_HOST_DEVICE int& fragmentID(size_t nodeID) const;
  SPHERAL_HOST_DEVICE int& particleType(size_t nodeID) const;

  // CHAI/ManagedArray specific operations for the contained FieldViews.
  SPHERAL_HOST void move(chai::ExecutionSpace space);
  SPHERAL_HOST void touch(chai::ExecutionSpace space);

protected:
  SymTensorFieldView mDeviatoricStressView;
  ScalarFieldView mPlasticStrainView;
  ScalarFieldView mPlasticStrainRateView;
  SymTensorFieldView mDamageView;
  IntFieldView mFragmentIDsView;
  IntFieldView mParticleTypesView;
};

}

#include "SolidNodeListViewInline.hh"

#endif
