//---------------------------------Spheral++----------------------------------//
// JohnsonCookFailureStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the JC failure strain.
//
// Created by JMO, Thu Jul 12 13:40:45 PDT 2018
//----------------------------------------------------------------------------//
#include "JohnsonCookFailureStrainPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookFailureStrainPolicy<Dimension>::
JohnsonCookFailureStrainPolicy(const Field<Dimension, Scalar>& D1,
                               const Field<Dimension, Scalar>& D2,
                               const double D3,
                               const double D4,
                               const double D5,
                               const double epsilondot0,
                               const double sigmamax,
                               const double efailmin,
                               const double Tcrit):
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::pressure,
                                        HydroFieldNames::specificThermalEnergy,
                                        SolidFieldNames::deviatoricStress,
                                        SolidFieldNames::plasticStrain,
                                        SolidFieldNames::meltSpecificEnergy}),
  mD1(D1),
  mD2(D2),
  mD3(D3),
  mD4(D4),
  mD5(D5),
  mepsilondot0(epsilondot0),
  msigmamax(sigmamax),
  mefailmin(efailmin),
  mTcrit(Tcrit) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookFailureStrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::flaws);
  auto& efail = state.field(key, 0.0);

  // Get the state fields.
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType epsKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  const KeyType stressKey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  const KeyType psrKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrainRate, nodeListKey);
  const KeyType meltKey = State<Dimension>::buildFieldKey(SolidFieldNames::meltSpecificEnergy, nodeListKey);
  CHECK(state.registered(PKey));
  CHECK(state.registered(epsKey));
  CHECK(state.registered(stressKey));
  CHECK(derivs.registered(psrKey));
  CHECK(state.registered(meltKey));

  const auto& P = state.field(PKey, 0.0);
  const auto& eps = state.field(epsKey, 0.0);
  const auto& S = state.field(stressKey, SymTensor::zero());
  const auto& psr = derivs.field(psrKey, 0.0);
  const auto& epsMelt = state.field(meltKey, 0.0);

  // Iterate over the internal nodes.
  const auto n = efail.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto sigmaVM = sqrt(1.5*S(i).doubledot(S(i)));
    CHECK(sigmaVM >= 0.0);
    const auto Pi = min(0.0, P(i));                        // Only use negative pressure
    if (sigmaVM < 1.0e-10 or abs(psr(i)) < 1.0e-10) {

      // Negligible von-Mises stress or negligible plastic strain.
      efail(i) = mD1(i) + mD2(i);

    } else if (-Pi < mTcrit/sigmaVM) {

      // Ordinary JC definition
      efail(i) = (mD1(i) + mD2(i)*exp(min(35.0, mD3*Pi/sigmaVM)))*
                 (1.0 + mD5*eps(i)*safeInv(epsMelt(i)));
      if (psr(i) > mepsilondot0) {
        efail(i) *= 1.0 + mD4*log(psr(i)*safeInv(mepsilondot0));
      }

    } else if (-Pi < -msigmamax) {

      // Linear interpolation to minimum failure
      auto efailTcrit = (mD1(i) + mD2(i)*exp(min(35.0, mD3*mTcrit)))*
                        (1.0 + mD5*eps(i)*safeInv(epsMelt(i)));
      if (psr(i) > mepsilondot0) {
        efail(i) *= 1.0 + mD4*log(psr(i)*safeInv(mepsilondot0));
      }
      const auto psi = max(0.0, min(1.0, (-Pi + mTcrit*sigmaVM)*safeInvVar(-msigmamax + mTcrit*sigmaVM)));
      CHECK(psi >= 0.0 and psi <= 1.0);
      efail(i) = (1.0 - psi)*efailTcrit + psi*mefailmin;

    } else {

      // Just the minimum
      efail(i) = mefailmin;

    }
    efail(i) = max(mefailmin, efail(i));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
JohnsonCookFailureStrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const JohnsonCookFailureStrainPolicy<Dimension>* rhsPtr = dynamic_cast<const JohnsonCookFailureStrainPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}
