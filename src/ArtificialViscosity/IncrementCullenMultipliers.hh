//---------------------------------Spheral++----------------------------------//
// IncrementCullenMultipliers
//
// Created by JMO, Mon Dec 28 22:18:44 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementCullenMultipliers_hh__
#define __Spheral_IncrementCullenMultipliers_hh__

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class IncrementCullenMultipliers: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  IncrementCullenMultipliers(const Scalar minValue, const Scalar maxValue, const bool hopkinsForm);
  virtual ~IncrementCullenMultipliers();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mHopkinsForm;
  Scalar mAlphaMin, mAlphaMax;
  IncrementCullenMultipliers(const IncrementCullenMultipliers& rhs);
  IncrementCullenMultipliers& operator=(const IncrementCullenMultipliers& rhs);
};

}

#endif
