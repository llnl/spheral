//---------------------------------Spheral++----------------------------------//
// HVolumePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the volume based on the local hull constructions.
//
// Created by JMO, Wed Aug 13 10:52:16 PDT 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_HVolumePolicy_hh__
#define __Spheral_HVolumePolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension>
class HVolumePolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  HVolumePolicy(const Scalar kernelExtent);
  virtual ~HVolumePolicy() {}
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // We'll make the updateAsIncrement a no-op.
  virtual void updateAsIncrement(const KeyType& /*key*/,
                                 State<Dimension>& /*state*/,
                                 StateDerivatives<Dimension>& /*derivs*/,
                                 const double /*multiplier*/,
                                 const double /*t*/,
                                 const double /*dt*/) override {}

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mKernelExtent;

  HVolumePolicy(const HVolumePolicy& rhs);
  HVolumePolicy& operator=(const HVolumePolicy& rhs);
};

}

#endif
