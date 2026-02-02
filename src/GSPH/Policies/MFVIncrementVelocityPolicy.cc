//---------------------------------Spheral++----------------------------------//
// MFVIncrementVelocityPolicy -- specialized policy for hydros that allow for mass
//                      flux between nodes. The momentum time derivative
//                      is used to update the velocity. The "hydro acceleration"
//                      is also added in to be compatible w/ phys packages
//                      that apply a pure acceleration. 
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
// TODO : HydroAcceleration needs to be added in
//----------------------------------------------------------------------------//

#include "GSPH/Policies/MFVIncrementVelocityPolicy.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

#include <limits.h>

namespace Spheral {
//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
MFVIncrementVelocityPolicy<Dimension>::
MFVIncrementVelocityPolicy(std::initializer_list<std::string> depends):
  FieldUpdatePolicy<Dimension, Vector>(depends){
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVIncrementVelocityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  const auto massKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const auto momDerivFieldKey = StateBase<Dimension>::buildFieldKey(prefix() + GSPHFieldNames::momentum, nodeListKey);
  //const auto accDerivFieldKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::hydroAcceleration, nodeListKey);
  
  const auto&  m = state.field(massKey, Scalar());
        auto&  v = state.field(key,     Vector());

  const auto&  DmDt = derivs.field(prefix() + massKey, Scalar());
  const auto&  DpDt = derivs.field(momDerivFieldKey,   Vector());
  //const auto&  DvDt = derivs.field(accDerivFieldKey,   Vector());

  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto m1 = m(i)+DmDt(i)*multiplier;
    const auto DpDti = DpDt(i);
    if (m1 > tiny) v(i) += (DpDti - DmDt(i)*v(i)) * multiplier * safeInv(m1);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MFVIncrementVelocityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const auto* rhsPtr = dynamic_cast<const MFVIncrementVelocityPolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

