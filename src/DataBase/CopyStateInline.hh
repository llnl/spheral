//---------------------------------Spheral++----------------------------------//
// CopyState -- An implementation of UpdatePolicyBase appropriate for
// copying one state field to another.
// Assumes that the copied state is dependent upon the master state, *and* that
// that is the only dependency.
//
// Created by JMO, Tue Oct 5 11:08:48 2004
//----------------------------------------------------------------------------//

#include "State.hh"
#include "StateDerivatives.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
CopyState<Dimension, ValueType>::
CopyState(const std::string& masterState,
          const std::string& copyState):
  UpdatePolicyBase<Dimension>(),
  mMasterStateName(masterState),
  mCopyStateName(copyState) {
}

//------------------------------------------------------------------------------
// Update the state
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
CopyState<Dimension, ValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  REQUIRE(key == mCopyStateName);

  // The state we're updating
  ValueType& f = state.template get<ValueType>(key);

  // The master state we're copying
  const ValueType& fmaster = state.template get<ValueType>(mMasterStateName);

  // Copy the master state using the assignment operator
  f = fmaster;
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
bool
CopyState<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const CopyState<Dimension, ValueType>*>(&rhs);
  if (rhsPtr == 0) return false;
  return (mMasterStateName == rhsPtr->mMasterStateName && 
          mCopyStateName == rhsPtr->mCopyStateName);
}

}

