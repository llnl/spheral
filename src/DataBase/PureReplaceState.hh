//---------------------------------Spheral++----------------------------------//
// PureReplaceState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PureReplaceState_hh__
#define __Spheral_PureReplaceState_hh__

#include "FieldUpdatePolicy.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class PureReplaceState: public FieldUpdatePolicy<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  PureReplaceState(std::initializer_list<std::string> depends = {});
  virtual ~PureReplaceState() {}
  
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

  // Don't use as evolved state for implicit integration
  virtual bool independent() const override { return false; }

private:
  //--------------------------- Private Interface ---------------------------//
  PureReplaceState(const PureReplaceState& rhs);
  PureReplaceState& operator=(const PureReplaceState& rhs);
};

}

#include "PureReplaceStateInline.hh"

#endif
