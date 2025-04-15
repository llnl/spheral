//---------------------------------Spheral++----------------------------------//
// SVPHMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the mass density in the state according to the SVPH 
// formalism.
//
// Created by JMO, Sat Aug 10 18:52:20 PDT 2013
//----------------------------------------------------------------------------//
#include "SVPH/SVPHMassDensityPolicy.hh"
#include "SVPH/SVPHFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHMassDensityPolicy<Dimension>::
SVPHMassDensityPolicy(const Scalar& rhoMin,
                      const Scalar& rhoMax):
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::mass,
                                        SVPHFieldNames::A_SVPH}),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax) {
}

//------------------------------------------------------------------------------
// Update the mass density using our SVPH defintion for momentum conservation.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHMassDensityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(fieldKey == HydroFieldNames::massDensity);
  const KeyType Akey = StateBase<Dimension>::buildFieldKey(SVPHFieldNames::A_SVPH, nodeListKey);
  const KeyType massKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);

  Field<Dimension, Scalar>& rho = state.field(key, 0.0);
  const Field<Dimension, Scalar>& A = state.field(Akey, 0.0);
  const Field<Dimension, Scalar>& mass = state.field(Akey, 0.0);

  const NodeList<Dimension>& nodeList = rho.nodeList();
  const unsigned n = nodeList.numInternalNodes();
  for (unsigned i = 0; i != n; ++i) {
    CHECK(A(i) > 0.0);
    CHECK(mass(i) > 0.0);
    rho(i) = std::max(mRhoMin, std::min(mRhoMax, A(i)*mass(i)));
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SVPHMassDensityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a SVPHMassDensity object.
  const SVPHMassDensityPolicy<Dimension>* rhsPtr = dynamic_cast<const SVPHMassDensityPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

