//---------------------------------Spheral++----------------------------------//
// TensorStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor strain.
//
// Created by JMO, Mon Oct 17 10:56:28 PDT 2005
//----------------------------------------------------------------------------//
#include "TensorStrainPolicy.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Kernel/TableKernel.hh"

#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorStrainPolicy<Dimension>::
TensorStrainPolicy(const TensorStrainAlgorithm strainType):
  FieldUpdatePolicy<Dimension, SymTensor>({HydroFieldNames::position,
                                           HydroFieldNames::H,
                                           SolidFieldNames::YoungsModulus,
                                           SolidFieldNames::bulkModulus,
                                           SolidFieldNames::shearModulus,
                                           HydroFieldNames::pressure,
                                           SolidFieldNames::deviatoricStress}),
  mStrainType(strainType) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorStrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::effectiveStrainTensor);
  auto& stateField = state.field(key, SymTensor::zero());

  const auto tiny = 1.0e-15;

  // Alias for shorter call building State Field keys
  auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };

  // Get the state fields.
  auto&       strain = state.field(buildKey(SolidFieldNames::strainTensor), SymTensor::zero());
  const auto& E = state.field(buildKey(SolidFieldNames::YoungsModulus), 0.0);
  const auto& K = state.field(buildKey(SolidFieldNames::bulkModulus), 0.0);
  const auto& mu = state.field(buildKey(SolidFieldNames::shearModulus), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& plasticStrain = state.field(buildKey(SolidFieldNames::plasticStrain), 0.0);
  const auto& S = state.field(buildKey(SolidFieldNames::deviatoricStress), SymTensor::zero());
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());
  const auto& gradv = derivs.field(buildKey(HydroFieldNames::internalVelocityGradient), Tensor::zero());
  const auto& DSDt = derivs.field(buildKey(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress), SymTensor::zero());

  // Check if a porosity model has registered a modifier for the deviatoric stress.
  // They should have added it as a dependency of this policy if so.
  const auto porosityScaling = state.registered(buildKey(SolidFieldNames::fDSjutzi));
  const Field<Dimension, Scalar>* fDSptr = nullptr;
  const Field<Dimension, Scalar>* alphaPtr = nullptr;
  const Field<Dimension, Scalar>* DalphaDtPtr = nullptr;
  if (porosityScaling) {
    fDSptr = &state.field(buildKey(SolidFieldNames::fDSjutzi), 0.0);
    alphaPtr = &state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    DalphaDtPtr = &derivs.field(buildKey(IncrementBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityAlpha), 0.0);
  }

  // Iterate over the internal nodes.
  const auto ni = stateField.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    double fDSi = 1.0;

    // Begin the big bonanza of options!

    // PseudoPlasticStrain.
    if (mStrainType == TensorStrainAlgorithm::PseudoPlasticStrain) {

      auto DSDti = DSDt(i);
      if (porosityScaling) {
        fDSi = (*fDSptr)(i);
        const auto alphai = (*alphaPtr)(i);
        const auto DalphaDti = (*DalphaDtPtr)(i);
        CHECK(alphai >= 1.0);
        DSDti = (fDSi*DSDti - S(i)*DalphaDti/alphai)/alphai;
      }
      strain(i) += multiplier*safeInv(mu(i), 1.0e-10)*DSDti;
      stateField(i) = strain(i);

    } else {

      // First apply the rotational term to the current strain history.
      const auto gradvi = fDSi * gradv(i);
      const auto spin = gradvi.SkewSymmetric();
      strain(i) += multiplier*(spin*strain(i) + strain(i)*spin).Symmetric();

      // Update the strain history with the current instantaneous deformation.
      const auto eigenv = gradvi.Symmetric().eigenVectors();
      auto sgradvi = constructSymTensorWithMaxDiagonal(eigenv.eigenValues, 0.0);
      sgradvi.rotationalTransform(eigenv.eigenVectors);
      strain(i) += multiplier*sgradvi;

      const auto volstrain = strain(i).Trace();

      // Update the effective strain according to the specified algorithm.
      switch(mStrainType) {

      case(TensorStrainAlgorithm::BenzAsphaugStrain):
        CHECK2(E(i) >= 0.0, "Bad Youngs modulus for " << stateField.nodeList().name() << " " << i << " : " << E(i));
        stateField(i) = (S(i) - P(i)*SymTensor::one())/(E(i) + tiny);
        break;

      case(TensorStrainAlgorithm::StrainHistory):
        stateField(i) = strain(i);
        break;

      case(TensorStrainAlgorithm::MeloshRyanAsphaugStrain):
        stateField(i) = ((K(i) - 2.0*mu(i)/Dimension::nDim)*volstrain*SymTensor::one() + 2.0*mu(i)*strain(i))/(E(i) + tiny);
        break;

      case(TensorStrainAlgorithm::PlasticStrain):
        stateField(i) = plasticStrain(i)*SymTensor::one();
        break;

      default:
        VERIFY2(false, "TensorStrainPolicy ERROR:  no update for case " << static_cast<int>(mStrainType) << "!");
        break;

      }
    }

    //------------------------------------------------------------------------------------
    // NOTE: this is a temporary fix to address the difference between Spheral's damaged
    //       pressure and the standard EOS pressure originally used by Benz and Asphaug.
    //       (1-D) is squared to counter-act the 1-D rolled into the pressure. Original 
    //        code is commented out below. (This will enhance the deviatoric portion)
    //       -JMPearl 11/14/2023

    //const auto fDs = 1.0 - D(i).eigenValues().maxElement();
    //stateField(i) *= safeInvVar(max(0.0, fDs*fDs), tiny);

    // Damage enhancement of the effective strain.
    stateField(i) *= safeInvVar(max(0.0, 1.0 - D(i).Trace()/Dimension::nDim), tiny);
    //------------------------------------------------------------------------------------


    // Apply limiting to the effective strain.
    stateField(i) = max(stateField(i), 1.0e-7*max(1.0, std::abs(stateField(i).Trace())/Dimension::nDim));
    // ENSURE2(fuzzyGreaterThanOrEqual(stateField(i).eigenValues().minElement(), 0.0, 1.0e-5),
    //         "Effective strain bad eigenvalues!  " << stateField(i).eigenValues());

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TensorStrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const TensorStrainPolicy<Dimension>*>(&rhs) != nullptr;
}

}

