//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyVolumePolicy
//
// Created by JMO, Sun Aug 18 19:21:27 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SpecificThermalEnergyVolumePolicy_hh__
#define __Spheral_SpecificThermalEnergyVolumePolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

template<typename Dimension>
class SpecificThermalEnergyVolumePolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  SpecificThermalEnergyVolumePolicy();
  virtual ~SpecificThermalEnergyVolumePolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  SpecificThermalEnergyVolumePolicy(const SpecificThermalEnergyVolumePolicy& rhs) = delete;
  SpecificThermalEnergyVolumePolicy& operator=(const SpecificThermalEnergyVolumePolicy& rhs) = delete;
};

}

#endif
