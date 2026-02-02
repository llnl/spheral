//---------------------------------Spheral++----------------------------------//
// NonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method in the case where we do *not* assume that pairwise forces are
// equal and opposite.
//
// Created by JMO, Sat Aug 10 23:03:39 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_NonSymmetricSpecificThermalEnergyPolicy_hh__
#define __Spheral_NonSymmetricSpecificThermalEnergyPolicy_hh__

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
class NonSymmetricSpecificThermalEnergyPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  NonSymmetricSpecificThermalEnergyPolicy(const DataBase<Dimension>& db);
  virtual ~NonSymmetricSpecificThermalEnergyPolicy() = default;
  
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
  NonSymmetricSpecificThermalEnergyPolicy(const NonSymmetricSpecificThermalEnergyPolicy& rhs) = delete;
  NonSymmetricSpecificThermalEnergyPolicy& operator=(const NonSymmetricSpecificThermalEnergyPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>* mDataBasePtr;
};

}

#endif
