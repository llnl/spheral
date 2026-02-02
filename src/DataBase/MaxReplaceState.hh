//---------------------------------Spheral++----------------------------------//
// MaxReplaceState -- Replaces the state with the max of the state/deriv value
//
// J.M. Pearl 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_MaxReplaceState_hh__
#define __Spheral_MaxReplaceState_hh__

#include "FieldUpdatePolicy.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class MaxReplaceState: public FieldUpdatePolicy<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  MaxReplaceState(std::initializer_list<std::string> depends = {});
  virtual ~MaxReplaceState() {}
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "new "; }

private:
  //--------------------------- Private Interface ---------------------------//
  MaxReplaceState(const MaxReplaceState& rhs);
  MaxReplaceState& operator=(const MaxReplaceState& rhs);
};

}

#include "MaxReplaceStateInline.hh"

#endif
