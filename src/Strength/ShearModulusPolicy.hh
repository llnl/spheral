//---------------------------------Spheral++----------------------------------//
// ShearModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent shear modulus state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ShearModulusPolicy_hh__
#define __Spheral_ShearModulusPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class ShearModulusPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  ShearModulusPolicy();
  virtual ~ShearModulusPolicy() = default;
  
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
  ShearModulusPolicy(const ShearModulusPolicy& rhs) = delete;
  ShearModulusPolicy& operator=(const ShearModulusPolicy& rhs) = delete;
};

}

#endif
