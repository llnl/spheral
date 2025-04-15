//---------------------------------Spheral++----------------------------------//
// SpecificFromTotalThermalEnergyPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the specific thermal energy from the total
// energy.
// 
// Created by JMO, Thu Jan  7 23:08:39 PST 2016
//----------------------------------------------------------------------------//

#ifndef __Spheral_SpecificFromTotalThermalEnergyPolicy_hh__
#define __Spheral_SpecificFromTotalThermalEnergyPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class SpecificFromTotalThermalEnergyPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  SpecificFromTotalThermalEnergyPolicy();
  virtual ~SpecificFromTotalThermalEnergyPolicy() = default;
  
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
  SpecificFromTotalThermalEnergyPolicy(const SpecificFromTotalThermalEnergyPolicy& rhs) = delete;
  SpecificFromTotalThermalEnergyPolicy& operator=(const SpecificFromTotalThermalEnergyPolicy& rhs) = delete;
};

}

#endif
