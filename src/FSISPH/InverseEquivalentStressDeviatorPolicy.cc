//---------------------------------Spheral++----------------------------------//
// InverseEquivalentStressDeviatorPolicy -- updates the inverse equivalent stress
//      deviator for the FSISPH hydro package. This allows yield applied in
//      a pairwise manner during the eval derivs loop.
//
// J.M. Pearl 2024
//----------------------------------------------------------------------------//

#include "InverseEquivalentStressDeviatorPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"

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
namespace {

//------------------------------------------------------------------------------
// Helper method to compute the J2 constant from the deviatoric stress.
//------------------------------------------------------------------------------
inline
double
computeJ2(const Dim<3>::SymTensor& S) {
  return 0.5*(S.doubledot(S));
}

inline
double
computeJ2(const Dim<2>::SymTensor& S) {
  // the third diagonal component of S (acutally negative to make Tr(S)=0)
  const auto S33 = S.Trace();
  return 0.5*(S.doubledot(S) + S33*S33);
}

inline
double
computeJ2(const Dim<1>::SymTensor& S) {
  // S_22 == S_33 = -S_11/2
  return 0.75*S.xx()*S.xx();
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
InverseEquivalentStressDeviatorPolicy<Dimension>::
InverseEquivalentStressDeviatorPolicy():
  FieldUpdatePolicy<Dimension,typename Dimension::Scalar>({SolidFieldNames::plasticStrain}) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
InverseEquivalentStressDeviatorPolicy<Dimension>::
~InverseEquivalentStressDeviatorPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
// Note this updates both the deviatoric stress and the plastic strain
// via the von Mises yielding correction.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InverseEquivalentStressDeviatorPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == FSIFieldNames::inverseEquivalentDeviatoricStress);

  // Get the state we depend on.
  const auto buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  auto&       invSeff = state.field(key, 0.0);
  const auto& deviatoricStress = state.field(buildKey(SolidFieldNames::deviatoricStress), SymTensor::zero());

  // Iterate over the internal nodes.
  const auto n = invSeff.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {

    // Equivalent stress deviator.
    const auto J2 = computeJ2(deviatoricStress(i));
    CHECK(J2 >= 0.0);

    // Scale the stress deviator.
    invSeff(i) = safeInvVar(sqrt(3.0*J2));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
InverseEquivalentStressDeviatorPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const InverseEquivalentStressDeviatorPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

