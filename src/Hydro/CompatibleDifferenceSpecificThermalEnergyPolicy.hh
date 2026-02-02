//---------------------------------Spheral++----------------------------------//
// CompatibleDifferenceSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy 
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method used in FSISPH as described in
//
// Pearl, J. M., Raskin, C. D., & Michael Owen, J. (2022). FSISPH: An SPH
// formulation for impacts between dissimilar materials. Journal of
// Computational Physics, 469, 111533.
//----------------------------------------------------------------------------//

#ifndef __Spheral_CompatibleDifferenceSpecificThermalEnergyPolicy_hh__
#define __Spheral_CompatibleDifferenceSpecificThermalEnergyPolicy_hh__

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
class CompatibleDifferenceSpecificThermalEnergyPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  CompatibleDifferenceSpecificThermalEnergyPolicy(const DataBase<Dimension>& db);
  virtual ~CompatibleDifferenceSpecificThermalEnergyPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // If the derivative stored values for the pair-accelerations has not been updated,
  // we need to just time advance normally.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Don't advance this policy implicitly
  virtual bool independent() const override { return false; }

  // Forbidden methods
  CompatibleDifferenceSpecificThermalEnergyPolicy(const CompatibleDifferenceSpecificThermalEnergyPolicy& rhs) = delete;
  CompatibleDifferenceSpecificThermalEnergyPolicy& operator=(const CompatibleDifferenceSpecificThermalEnergyPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>* mDataBasePtr;
};

}

#endif
