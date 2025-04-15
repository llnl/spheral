//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#include "ProbabilisticDamagePolicy.hh"
#include "ProbabilisticDamageModel.hh"
#include "oneMinusDamage.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
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
// Force all eigen values positive.
//------------------------------------------------------------------------------
inline
void
abs_in_place(Dim<1>::Vector& vec) {
  vec[0] = std::abs(vec[0]);
}

inline
void
abs_in_place(Dim<2>::Vector& vec) {
  vec[0] = std::abs(vec[0]);
  vec[1] = std::abs(vec[1]);
}

inline
void
abs_in_place(Dim<3>::Vector& vec) {
  vec[0] = std::abs(vec[0]);
  vec[1] = std::abs(vec[1]);
  vec[2] = std::abs(vec[2]);
}

//------------------------------------------------------------------------------
// Sort the eigen values (& associated eigen vectors) in decreasing order.
//------------------------------------------------------------------------------
inline
void
sortEigen(Dim<1>::SymTensor::EigenStructType&) {
}

inline
void
sortEigen(Dim<2>::SymTensor::EigenStructType& eigeni) {
  if (eigeni.eigenValues(0) < eigeni.eigenValues(1)) {
    std::swap(eigeni.eigenValues(0), eigeni.eigenValues(1));
    eigeni.eigenVectors = Dim<2>::Tensor(eigeni.eigenVectors.xy(), eigeni.eigenVectors.xx(),
                                         eigeni.eigenVectors.yy(), eigeni.eigenVectors.yx());
  }
  ENSURE(eigeni.eigenValues(0) >= eigeni.eigenValues(1));
}

inline
void
sortEigen(Dim<3>::SymTensor::EigenStructType& eigeni) {
  int i = 0;
  int j = 1;
  int k = 2;
  bool flipped = true;
  while (flipped) {
    flipped = false;
    if (eigeni.eigenValues(i) < eigeni.eigenValues(j)) {
      flipped = true;
      std::swap(i, j);
    }
    if (eigeni.eigenValues(j) < eigeni.eigenValues(k)) {
      flipped = true;
      std::swap(j, k);
    }
  }
  eigeni.eigenValues = Dim<3>::Vector(eigeni.eigenValues(i),
                                      eigeni.eigenValues(j),
                                      eigeni.eigenValues(k));
  eigeni.eigenVectors = Dim<3>::Tensor(eigeni.eigenVectors(0,i), eigeni.eigenVectors(0,j), eigeni.eigenVectors(0,k),
                                       eigeni.eigenVectors(1,i), eigeni.eigenVectors(1,j), eigeni.eigenVectors(1,k),
                                       eigeni.eigenVectors(2,i), eigeni.eigenVectors(2,j), eigeni.eigenVectors(2,k));
  ENSURE((eigeni.eigenValues(0) >= eigeni.eigenValues(1)) &&
         (eigeni.eigenValues(1) >= eigeni.eigenValues(2)));
}

//------------------------------------------------------------------------------
// Determine the effective rotational transformation from the given velocity
// gradient.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
effectiveRotation(const Dim<1>::Tensor&) {
  return Dim<1>::Tensor::one;
}

inline
Dim<2>::Tensor
effectiveRotation(const Dim<2>::Tensor& DvDx) {
  const double theta = DvDx.xy() - DvDx.yx();
  return Dim<2>::Tensor(cos(theta), sin(theta),
                        -sin(theta), cos(theta));
}

inline
Dim<3>::Tensor
effectiveRotation(const Dim<3>::Tensor& DvDx) {
  const Dim<3>::Vector theta(DvDx.yz() - DvDx.zy(),
                             DvDx.zx() - DvDx.xz(),
                             DvDx.xy() - DvDx.yx());
  return (Dim<3>::Tensor(cos(theta.x()), sin(theta.x()), 0.0,
                         -sin(theta.x()), cos(theta.x()), 0.0,
                         0.0, 0.0, 1.0)*
          Dim<3>::Tensor(cos(theta.y()), sin(theta.y()), 0.0,
                         0.0, 1.0, 0.0,
                         -sin(theta.y()), 0.0, cos(theta.y()))*
          Dim<3>::Tensor(1.0, 0.0, 0.0,
                         0.0, cos(theta.z()), sin(theta.z()),
                         0.0, -sin(theta.z()), cos(theta.z())));
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamagePolicy<Dimension>::
ProbabilisticDamagePolicy(const bool damageInCompression,  // allow damage in compression
                          const double kWeibull,           // coefficient in Weibull power-law
                          const double mWeibull):          // exponenent in Weibull power-law
  FieldUpdatePolicy<Dimension, SymTensor>({SolidFieldNames::strain}),
  mDamageInCompression(damageInCompression),
  mkWeibull(kWeibull),
  mmWeibull(mWeibull) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamagePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::tensorDamage);
  auto& stateField = state.field(key, SymTensor::zero);

  const auto Dtiny = 0.01;
  const auto Dtiny1 = 1.0/(FastMath::CubeRootHalley2(1.0 - Dtiny) - FastMath::CubeRootHalley2(Dtiny));

  // Alias for shorter call building State Field keys
  auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };

  // Get the state fields.
  const auto& strain = state.field(buildKey(SolidFieldNames::effectiveStrainTensor), SymTensor::zero);
  const auto& DDDt = derivs.field(buildKey(this->prefix() + SolidFieldNames::scalarDamage), 0.0);
  const auto& localDvDx = derivs.field(buildKey(HydroFieldNames::internalVelocityGradient), Tensor::zero);
  const auto& numFlaws = state.field(buildKey(SolidFieldNames::numFlaws), 0);
  const auto& minFlaw = state.field(buildKey(SolidFieldNames::minFlaw), 0.0);
  const auto& maxFlaw = state.field(buildKey(SolidFieldNames::maxFlaw), 0.0);
  const auto& V0 = state.field(buildKey(SolidFieldNames::initialVolume), 0.0);

  // Check if porosity is active for this material
  const auto usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));
  Field<Dimension, Scalar>* alpha0Ptr = nullptr;
  Field<Dimension, Scalar>* alphaPtr = nullptr;
  Field<Dimension, Scalar>* DalphaDtPtr = nullptr;
  if (usePorosity) {
    alpha0Ptr = &state.field(buildKey(SolidFieldNames::porosityAlpha0), 0.0);
    alphaPtr = &state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    DalphaDtPtr = &derivs.field(buildKey(IncrementBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityAlpha), 0.0);
  }

  // Iterate over the internal nodes.
  const auto ni = stateField.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    auto& Di = stateField(i);
    CHECK(Di.eigenValues().minElement() >= 0.0 and
          fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));

    // First apply the rotational term to the current damage.
    const auto spin = localDvDx(i).SkewSymmetric();
    const auto spinCorrection = (spin*Di - Di*spin).Symmetric();
    Di += multiplier*spinCorrection;
    Di = max(1.0e-5, min(1.0 - 2.0e-5, Di));
    {
      const auto maxValue = Di.eigenValues().maxElement();
      if (maxValue > 1.0) Di /= maxValue;
    }
    CHECK(Di.eigenValues().minElement() >= 0.0 and
          fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));

    // The tensor strain on this node.
    auto eigeni = strain(i).eigenVectors();

    // If we're allowing damage in compression, force all strains to be abs value.
    if (mDamageInCompression) abs_in_place(eigeni.eigenValues);

    // We want to go over these from max to min.
    sortEigen(eigeni);

    // Iterate over the eigen values/vectors of the strain.
    for (auto j = 0; j < Dimension::nDim; ++j) {

      // The direction of the strain, and projected current (starting) damage.
      const auto strainDirection = eigeni.eigenVectors.getColumn(j);
      CHECK(fuzzyEqual(strainDirection.magnitude2(), 1.0, 1.0e-10));
      const auto D0 = max(0.0, min(1.0, (Di * strainDirection).magnitude()));
      CHECK(D0 >= 0.0 && D0 <= 1.0);
      const auto strainj = std::max(0.0, eigeni.eigenValues(j)*safeInvVar(1.0 - D0)); //   + plasticStrain(i))/(fDi*fDi + 1.0e-20);

      if (D0 < 1.0) {

        // Find how many flaws are activated in this direction
        CHECK(maxFlaw(i) > 0.0);
        CHECK(maxFlaw(i) > minFlaw(i));
        const auto fractionFlawsActive = (numFlaws(i) == 0 or strainj <= minFlaw(i) ? 0.0 :
                                          strainj >= maxFlaw(i) ? 1.0 :
                                          mkWeibull*V0(i)*(pow(strainj, mmWeibull) - pow(minFlaw(i), mmWeibull))/numFlaws(i));
        CHECK(fractionFlawsActive >= 0.0 and fractionFlawsActive <= 1.0);

        // Choose the allowed range of D.
        const auto D013 = FastMath::CubeRootHalley2(D0);
        double D13min, D13max;
        if (multiplier >= 0.0) {
          D13min = D013;
          D13max = std::max(D013, FastMath::CubeRootHalley2(fractionFlawsActive));
        } else {
          D13min = 0.0;
          D13max = D013;
        }
        CHECK(D13max >= D13min);
        CHECK(D13min >= 0.0 && D13max <= 1.0);

        // Increment the damage due to tensile stresses
        auto D113 = std::max(D13min, std::min(D13max, D013 + multiplier*fractionFlawsActive*DDDt(i)));

        // In the presence of porosity, we can also accumulate damage from crushing pores
        if (usePorosity) {
          CHECK((*alpha0Ptr)(i) >= 1.0);
          const auto alpha0 = (*alpha0Ptr)(i);
          const auto alpha = (*alphaPtr)(i);
          const auto DalphaDti = std::min(0.0, (*DalphaDtPtr)(i));   // Only allowed to grow damage, not reduce it.
          const auto phi0 = 1.0 - 1.0/alpha0;
          const auto DD13Dt_p = -FastMath::CubeRootHalley2(phi0)*safeInv(3.0 * pow(1.0 - (alpha - 1.0)*safeInv(alpha0 - 1.0) + Dtiny, 2.0/3.0) * (alpha0 - 1.0))*Dtiny1*DalphaDti;
          CHECK2(DD13Dt_p >= 0.0, "bad DD13Dt_p: " << DD13Dt_p);
          D113 = std::min(1.0, D113 + multiplier*DD13Dt_p);
        }

        // Increment the damage tensor.
        const auto D1 = std::max(0.0, std::min(1.0, FastMath::cube(D113)));
        const auto deltaD = std::abs(D1 - D0)*sgn(multiplier);
        CHECK2(deltaD*multiplier >= 0.0, D0 << " " << D1 << " " << multiplier << " " << deltaD*multiplier);
        Di += deltaD*strainDirection.selfdyad();
      }
    }

    // Enforce bounds on the damage.
    const auto Dvals = Di.eigenValues();
    if (Dvals.minElement() < 0.0 or
        Dvals.maxElement() > 1.0) {
      const auto Deigen = Di.eigenVectors();
      Di = constructSymTensorWithBoundedDiagonal(Deigen.eigenValues,
                                                 1.0e-5,
                                                 1.0);
      Di.rotationalTransform(Deigen.eigenVectors);
    }
    ENSURE(Di.eigenValues().minElement() >= 0.0 and
           fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    //     ENSURE(fuzzyGreaterThanOrEqual(Di.eigenValues().minElement(), 0.0, 1.0e-5) &&
    //            fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    //     ENSURE(fuzzyGreaterThanOrEqual(Di.eigenValues().minElement(), 0.0) &&
    //            fuzzyLessThanOrEqual(Di.Trace(), 1.0));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ProbabilisticDamagePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ProbabilisticDamagePolicy<Dimension>* rhsPtr = dynamic_cast<const ProbabilisticDamagePolicy<Dimension>*>(&rhs);
  if (rhsPtr == nullptr) {
    return false;
  } else {
    return true;
  }
}

}

