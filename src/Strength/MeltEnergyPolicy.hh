//---------------------------------Spheral++----------------------------------//
// MeltEnergyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent melt energy state.
//
// Created by JMO, Fri Jul 13 12:57:54 PDT 2018
//----------------------------------------------------------------------------//
#ifndef __Spheral_MeltEnergyPolicy_hh__
#define __Spheral_MeltEnergyPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class MeltEnergyPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  MeltEnergyPolicy();
  virtual ~MeltEnergyPolicy();
  
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
  MeltEnergyPolicy(const MeltEnergyPolicy& rhs) = delete;
  MeltEnergyPolicy& operator=(const MeltEnergyPolicy& rhs) = delete;
};

}

#endif
