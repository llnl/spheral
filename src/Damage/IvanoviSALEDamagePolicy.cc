//---------------------------------Spheral++----------------------------------//
// IvanoviSALEDamagePolicy
//
// The Ivanov damage update policyt, hopefully close to how it's implemented in iSALE.
// This damage model is most appropriate for rocky materials.
//
// Refs:
// 
// Collins, G. S., Melosh, H. J., & Ivanov, B. A. (2004). Modeling damage and deformation in impact simulations.
//   Meteoritics & Planetary Science, 39(2), 217–231. http://doi.wiley.com/10.1111/j.1945-5100.2004.tb00337.x
//
// Raducan, S. D., Davison, T. M., Luther, R., & Collins, G. S. (2019). The role of asteroid strength, porosity and
//   internal friction in impact momentum transfer. Icarus. https://doi.org/10.1016/J.ICARUS.2019.03.040
//
// Lundborg, N. (1967). The strength-size relation of granite. International Journal of Rock Mechanics and
//   Mining Sciences & Geomechanics Abstracts, 4(3):269–272.f
//
// Created by JMO, Mon Jun 28 12:14:37 PDT 2021
//----------------------------------------------------------------------------//
#include "IvanoviSALEDamagePolicy.hh"
#include "oneMinusDamage.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
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
  return Dim<1>::Tensor::one();
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
IvanoviSALEDamagePolicy<Dimension>::
IvanoviSALEDamagePolicy(const double minPlasticFailure,             // minimum plastic strain for failure
                        const double plasticFailurePressureSlope,   // slope for critical plastic strain
                        const double plasticFailurePressureOffset,  // intercept for critical plastic strain
                        const double tensileFailureStress):         // threshold for tensile failure
  FieldUpdatePolicy<Dimension, SymTensor>({SolidFieldNames::strain}),
  mEpsPfb(minPlasticFailure),
  mB(plasticFailurePressureSlope),
  mPc(plasticFailurePressureOffset),
  mTensileFailureStress(tensileFailureStress) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamagePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::tensorDamage);
  auto& stateField = state.field(key, SymTensor::zero());

  // Get the state fields.
  const auto strainKey = State<Dimension>::buildFieldKey(SolidFieldNames::effectiveStrainTensor, nodeListKey);
  const auto Skey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  const auto Pkey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const auto DdamageDtKey = State<Dimension>::buildFieldKey(this->prefix() + SolidFieldNames::scalarDamage, nodeListKey);
  const auto DvDxKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, nodeListKey);
  CHECK(state.registered(strainKey));
  CHECK(state.registered(Skey));
  CHECK(state.registered(Pkey));
  CHECK(derivs.registered(DdamageDtKey));
  CHECK(derivs.registered(DvDxKey));
  const auto& strain = state.field(strainKey, SymTensor::zero());
  const auto& S = state.field(Skey, SymTensor::zero());
  const auto& P = state.field(Pkey, 0.0);
  const auto& DDDt = derivs.field(DdamageDtKey, 0.0);
  const auto& localDvDx = derivs.field(DvDxKey, Tensor::zero());

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

    // First we look for the plastic strain based damage.  We are using the PsuedoPlasticStrain definition for strain,
    // which is a reasonable tensor generalization of the plastic strain used in iSALE.
    auto ps_eigeni = strain(i).eigenVectors();
    abs_in_place(ps_eigeni.eigenValues);
    sortEigen(ps_eigeni);

    // The threshold failure plastic strain.
    const auto epsf = max(mEpsPfb, mB*(P(i) - mPc));
    CHECK(epsf > 0.0);

    // Iterate over the plastic strain elements and increment the shear damage.
    for (auto jdim = 0; jdim < Dimension::nDim; ++jdim) {
      const auto strainj = ps_eigeni.eigenValues(jdim);
      if (strainj > epsf) {
        
        // The direction of the strain, and projected current (starting) damage.
        const auto strainDirection = ps_eigeni.eigenVectors.getColumn(jdim);
        CHECK(fuzzyEqual(strainDirection.magnitude2(), 1.0, 1.0e-10));
        const auto D0 = max(0.0, min(1.0, (Di * strainDirection).magnitude()));
        CHECK(D0 >= 0.0 && D0 <= 1.0);
        const auto D1 = min(1.0, max(D0, strainj/epsf));

        // Increment the damage in this direction by the plastic shear damage.
        Di += (D1 - D0)*strainDirection.selfdyad();
      }
    }

    // Now apply any tensile damage.
    const auto stressi = S(i) - P(i)*SymTensor::one();
    auto tensile_eigeni = stressi.eigenVectors();
    sortEigen(tensile_eigeni);

    // Apply the tensile failure.
    for (auto jdim = 0; jdim < Dimension::nDim; ++jdim) {
      const auto strainj = tensile_eigeni.eigenValues(jdim);
      if (strainj > mTensileFailureStress) {

        // The direction of the strain, and projected current (starting) damage.
        const auto strainDirection = ps_eigeni.eigenVectors.getColumn(jdim);
        CHECK(fuzzyEqual(strainDirection.magnitude2(), 1.0, 1.0e-10));
        const auto D0 = max(0.0, min(1.0, (Di * strainDirection).magnitude()));
        CHECK(D0 >= 0.0 && D0 <= 1.0);

        // Increment the damage.
        const auto D013 = FastMath::CubeRootHalley2(D0);
        const auto D113 = D013 + multiplier*DDDt(i);
        const auto D1 = FastMath::cube(D113);
        CHECK((D1 - D0)*multiplier >= 0.0);

        // Increment the damage tensor.
        Di += (D1 - D0)*strainDirection.selfdyad();

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
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
IvanoviSALEDamagePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IvanoviSALEDamagePolicy<Dimension>* rhsPtr = dynamic_cast<const IvanoviSALEDamagePolicy<Dimension>*>(&rhs);
  if (rhsPtr == nullptr) {
    return false;
  } else {
    return true;
  }
}

}

