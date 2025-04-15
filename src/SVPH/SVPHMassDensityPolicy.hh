//---------------------------------Spheral++----------------------------------//
// SVPHMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the mass density in the state according to the SVPH 
// formalism.
//
// Created by JMO, Sat Aug 10 18:52:20 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHMassDensityPolicy_hh__
#define __Spheral_SVPHMassDensityPolicy_hh__

#include <string>

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class SVPHMassDensityPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  SVPHMassDensityPolicy(const Scalar& rhoMin,
                        const Scalar& rhoMax);
  virtual ~SVPHMassDensityPolicy() = default;
  
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
  SVPHMassDensityPolicy() = delete;
  SVPHMassDensityPolicy(const SVPHMassDensityPolicy& rhs) = delete;
  SVPHMassDensityPolicy& operator=(const SVPHMassDensityPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mRhoMin, mRhoMax;
};

}

#endif
