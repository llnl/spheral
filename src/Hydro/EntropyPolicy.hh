//---------------------------------Spheral++----------------------------------//
// EntropyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the entropy.
//
// Created by JMO, Sun Mar  6 21:57:47 PST 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_EntropyPolicy_hh__
#define __Spheral_EntropyPolicy_hh__

#include <string>

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;

template<typename Dimension>
class EntropyPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  EntropyPolicy();
  virtual ~EntropyPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbiddnen
  EntropyPolicy(const EntropyPolicy& rhs) = delete;
  EntropyPolicy& operator=(const EntropyPolicy& rhs) = delete;
};

}

#endif
