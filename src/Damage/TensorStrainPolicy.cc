//---------------------------------Spheral++----------------------------------//
// TensorStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor strain.
//
// Created by JMO, Mon Oct 17 10:56:28 PDT 2005
//----------------------------------------------------------------------------//
#include "Damage/TensorStrainPolicy.hh"
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
#include "Geometry/GeometryRegistrar.hh"
#include "Geometry/RZGeometryOps.hh"
#include "Kernel/TableKernel.hh"

#include <functional>

namespace Spheral {

namespace {  // anonymous

// //------------------------------------------------------------------------------
// // Velocity gradient
// //------------------------------------------------------------------------------
// // Generic definition
// template<typename VectorType, typename TensorType, typename ReturnType>
// inline
// ReturnType
// velocityGradient(const TensorType& DvDxi,
//                  const VectorType& posi,
//                  const VectorType& veli,
//                  const ReturnType& dummy) {
//   return DvDxi;
// }

// // RZ
// template<>
// inline
// Dim<3>::Tensor
// velocityGradient<Dim<2>::Vector, Dim<2>::Tensor, Dim<3>::Tensor>(const Dim<2>::Tensor& DvDxi,
//                                                                  const Dim<2>::Vector& posi,
//                                                                  const Dim<2>::Vector& veli,
//                                                                  const Dim<3>::Tensor& dummy) {
//   return Dim<3>::Tensor(DvDxi[0], DvDxi[1], 0.0,
//                         DvDxi[2], DvDxi[3], 0.0,
//                         0.0,      0.0,      veli.y()*safeInvVar(posi.y(), 1.0e-10));
// }

// //------------------------------------------------------------------------------
// // Deviatoric stress
// //------------------------------------------------------------------------------
// // Generic definition
// template<typename Dimension, typename SymTensorType>
// inline
// SymTensorType
// deviatoricStress(const typename Dimension::SymTensor& Si) {
//   return Si;
// }

// // RZ
// template<>
// inline
// Dim<3>::SymTensor
// deviatoricStress<Dim<2>, Dim<3>::SymTensor>(const Dim<2>::SymTensor& Si) {
//   return Dim<3>::SymTensor(Si[0], Si[1], 0.0,
//                            Si[1], Si[2], 0.0,
//                            0.0,   0.0,   -(Si[0] + Si[2]));
// }

// //------------------------------------------------------------------------------
// // Build a tensor of the requested type
// //------------------------------------------------------------------------------
// // Generic definition
// template<typename Dimension, typename OutTensorType>
// OutTensorType
// buildTensor(const size_t i,
//             const Field<Dimension, typename Dimension::SymTensor>& tfield,
//             const Field<Dimension, typename Dimension::Scalar>* ttfieldptr) {
//   return tfield(i);
// }

// // RZ specialization
// template<>
// Dim<3>::SymTensor
// buildTensor<Dim<2>, Dim<3>::SymTensor>(const size_t i,
//                                        const Field<Dim<2>, Dim<2>::SymTensor>& tfield,
//                                        const Field<Dim<2>, Dim<2>::Scalar>* ttfieldptr) {
//   REQUIRE(ttfieldptr != nullptr);
//   const auto& ti = tfield(i);
//   return Dim<3>::SymTensor(ti[0], ti[1], 0.0,
//                            ti[1], ti[2], 0.0,
//                            0.0,   0.0,   (*ttfieldptr)(i));
// }


// //------------------------------------------------------------------------------
// // Assign tensor components
// //------------------------------------------------------------------------------
// // Generic definition
// template<typename Dimension, typename InTensorType>
// void
// assignTensor(const size_t i,
//              const InTensorType& ti,
//              Field<Dimension, typename Dimension::SymTensor>& tfield,
//              Field<Dimension, typename Dimension::Scalar>* ttfieldptr) {
//   tfield(i) = ti;
// }

// // RZ specialization
// template<>
// void
// assignTensor<Dim<2>, Dim<3>::SymTensor>(const size_t i,
//                                         const Dim<3>::SymTensor& ti,
//                                         Field<Dim<2>, Dim<2>::SymTensor>& tfield,
//                                         Field<Dim<2>, Dim<2>::Scalar>* ttfieldptr) {
//   REQUIRE(ttfieldptr != nullptr);
//   tfield(i)[0] = ti[0];
//   tfield(i)[1] = ti[1];
//   tfield(i)[2] = ti[3];
//   (*ttfieldptr)(i) = ti[5];
// }

}            // anonymous

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
// Update the field (dispatch)
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorStrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  if (GeometryRegistrar::coords() == CoordinateType::Cartesian) {
    updateImpl<&SymTensor::eigenVectors, &SymTensor::Trace, Dimension::nDim>(key, state, derivs, multiplier, t, dt);
  } else {
    updateImpl<&SymTensor::eigenVectors3D, &SymTensor::Trace3D, 3>(key, state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// Update the field (implementation)
//------------------------------------------------------------------------------
template<typename Dimension>
template<auto EigenVectorsMethod, auto TraceMethod, int effDims>
void
TensorStrainPolicy<Dimension>::
updateImpl(const KeyType& key,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs,
           const double multiplier,
           const double t,
           const double dt) {
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
  const auto& S = state.field(buildKey(SolidFieldNames::deviatoricStress), Dimension::SymTensor::zero());
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());
  const auto& gradv = derivs.field(buildKey(HydroFieldNames::internalVelocityGradient), Tensor::zero());
  const auto& DSDt = derivs.field(buildKey(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress), Dimension::SymTensor::zero());

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
    const auto& Si = S(i);
    auto& straini = strain(i);
    auto& effStraini = stateField(i);

    // Begin the big bonanza of options!

    // PseudoPlasticStrain.
    if (mStrainType == TensorStrainAlgorithm::PseudoPlasticStrain) {

      auto DSDti = DSDt(i);
      if (porosityScaling) {
        fDSi = (*fDSptr)(i);
        const auto alphai = (*alphaPtr)(i);
        const auto DalphaDti = (*DalphaDtPtr)(i);
        CHECK(alphai >= 1.0);
        DSDti = (fDSi*DSDti - Si*DalphaDti/alphai)/alphai;
      }
      straini += multiplier*safeInv(mu(i), 1.0e-10)*DSDti;
      effStraini = straini;

    } else {

      // First apply the rotational term to the current strain history.
      const auto gradvi = fDSi * gradv(i);
      const auto spin = gradvi.SkewSymmetric();
      straini += multiplier*(straini*spin - spin*straini).Symmetric();

      // Update the strain history with the current instantaneous deformation.
      // const auto eigenv = gradvi.Symmetric().eigenVectors3D();
      const auto eigenv = std::invoke(EigenVectorsMethod, gradvi.Symmetric());
      SymTensor sgradvi(eigenv.eigenValues, std::numeric_limits<double>::lowest(), 0.0);
      sgradvi.rotationalTransform(eigenv.eigenVectors);
      straini += multiplier*sgradvi;

      const auto volstrain = straini.Trace();

      // Update the effective strain according to the specified algorithm.
      switch(mStrainType) {

      case(TensorStrainAlgorithm::BenzAsphaugStrain):
        CHECK2(E(i) >= 0.0, "Bad Youngs modulus for " << E.nodeList().name() << " " << i << " : " << E(i));
        effStraini = (Si - P(i)*SymTensor::one())/(E(i) + tiny);
        break;

      case(TensorStrainAlgorithm::StrainHistory):
        effStraini = straini;
        break;

      case(TensorStrainAlgorithm::MeloshRyanAsphaugStrain):
        effStraini = ((K(i) - 2.0*mu(i)/effDims)*volstrain*SymTensor::one() + 2.0*mu(i)*straini)/(E(i) + tiny);
        break;

      case(TensorStrainAlgorithm::PlasticStrain):
        effStraini = plasticStrain(i)*SymTensor::one();
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
    const auto Davgi = std::invoke(TraceMethod, D(i))/effDims;
    effStraini *= safeInvVar(std::max(0.0, 1.0 - Davgi), tiny);
    //------------------------------------------------------------------------------------


    // Apply limiting to the effective strain.
    const auto effStrainAvgi = std::invoke(TraceMethod, effStraini)/effDims;
    effStraini = max(effStraini, 1.0e-7*max(1.0, std::abs(effStrainAvgi)));
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

