//---------------------------------Spheral++----------------------------------//
// ReplaceState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// This version assumes there is a derivative based update available.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceState_hh__
#define __Spheral_ReplaceState_hh__

#include "PureReplaceState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class ReplaceState: public PureReplaceState<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  ReplaceState(std::initializer_list<std::string> depends = {});
  virtual ~ReplaceState() {}
  
  // An alternate method to be called when you want to specify that the "Replace" information
  // in the derivatives is invalid, and instead the value should be treated as a time advancement
  // algorithm instead.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Don't use as evolved state for implicit integration
  virtual bool independent() const override { return false; }

private:
  //--------------------------- Private Interface ---------------------------//
  ReplaceState(const ReplaceState& rhs);
  ReplaceState& operator=(const ReplaceState& rhs);
};

}

#include "ReplaceStateInline.hh"

#endif
