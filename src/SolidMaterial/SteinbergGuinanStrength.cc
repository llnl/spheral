//---------------------------------Spheral++----------------------------------//
// SteinbergGuinanStrength -- Implements the Steinberg-Guinan strength model.
//   ** Need reference **
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#include "SteinbergGuinanStrength.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"
#include "SolidEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SteinbergGuinanStrength<Dimension>::
SteinbergGuinanStrength(const SolidEquationOfState<Dimension>& eos,
                        const double G0,     
                        const double Gmax,     
                        const double A,      
                        const double B,      
                        const double Y0,     
                        const double Ymax,
                        const double Yp,
                        const double beta,
                        const double gamma0, 
                        const double nhard,
                        const NinthOrderPolynomialFit& coldEnergyFit,
                        const NinthOrderPolynomialFit& meltEnergyFit):
  StrengthModel<Dimension>(),
  mEOSPtr(&eos),
  mG0(G0),
  mGmax(Gmax),
  mA(A),
  mB(B),
  mY0(Y0),
  mYmax(Ymax),
  mYp(Yp),
  mbeta(beta),
  mgamma0(gamma0),
  mnhard(nhard),
  mColdEnergyFit(coldEnergyFit),
  mMeltEnergyFit(meltEnergyFit) {
}

//------------------------------------------------------------------------------
// Constructor (backwards compatible one without Gmax).
//------------------------------------------------------------------------------
template<typename Dimension>
SteinbergGuinanStrength<Dimension>::
SteinbergGuinanStrength(const SolidEquationOfState<Dimension>& eos,
                        const double G0,     
                        const double A,      
                        const double B,      
                        const double Y0,     
                        const double Ymax,
                        const double Yp,
                        const double beta,
                        const double gamma0, 
                        const double nhard,
                        const NinthOrderPolynomialFit& coldEnergyFit,
                        const NinthOrderPolynomialFit& meltEnergyFit):
  StrengthModel<Dimension>(),
  mEOSPtr(&eos),
  mG0(G0),
  mGmax(1e100),
  mA(A),
  mB(B),
  mY0(Y0),
  mYmax(Ymax),
  mYp(Yp),
  mbeta(beta),
  mgamma0(gamma0),
  mnhard(nhard),
  mColdEnergyFit(coldEnergyFit),
  mMeltEnergyFit(meltEnergyFit) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SteinbergGuinanStrength<Dimension>::
~SteinbergGuinanStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& specificThermalEnergy,
             const Field<Dimension, Scalar>& pressure,
             const Field<Dimension, SymTensor>& damage) const {
  const auto n = density.numInternalElements();
  if ((mG0 > 0.0) &&
      (mA == 0.0) &&
      (mB == 0.0) &&
      (mbeta == 0.0) &&
      (mnhard == 0.0)) {
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
      shearModulus(i) = (1.0 - Di)*mG0;
    }

  } else {
    Field<Dimension, Scalar> T("temperature", density.nodeList());
    this->computeTemperature(T, density, specificThermalEnergy);
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto eta = mEOSPtr->boundedEta(density(i));
      CHECK(distinctlyGreaterThan(eta, 0.0));
      const auto fmelt = meltAttenuation(density(i), specificThermalEnergy(i));
      shearModulus(i) = min(mGmax, 
                            mG0*max(1.0e-10,
                                    fmelt*(1.0 + mA*pressure(i)/FastMath::CubeRootHalley2(eta) - mB*T(i))));
      const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
      shearModulus(i) = (1.0 - Di)*shearModulus(i) + Di*fmelt*mG0;
      CHECK(shearModulus(i) > 0.0);
    }
  }
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& specificThermalEnergy,
              const Field<Dimension, Scalar>& pressure,
              const Field<Dimension, Scalar>& plasticStrain,
              const Field<Dimension, Scalar>& /*plasticStrainRate*/,
              const Field<Dimension, SymTensor>& damage) const {
  const auto n = density.numInternalElements();
  if ((mG0 > 0.0) &&
      (mA == 0.0) &&
      (mB == 0.0) &&
      (mbeta == 0.0) &&
      (mnhard == 0.0)) {

#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
      yieldStrength(i) = (1.0 - Di)*mY0;
    }
    
  } else {

    this->shearModulus(yieldStrength, density, specificThermalEnergy, pressure, damage);
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto eta = mEOSPtr->boundedEta(density(i));
      CONTRACT_VAR(eta);
      CHECK(distinctlyGreaterThan(eta, 0.0));
      const auto Yhard = min(mYmax, mY0*pow(1.0 + mbeta*(plasticStrain(i) + mgamma0), mnhard));
      const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
      yieldStrength(i) = (1.0 - Di)*Yhard*yieldStrength(i)/mG0;
    }
  }
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
soundSpeed(Field<Dimension, Scalar>& soundSpeed,
           const Field<Dimension, Scalar>& density,
           const Field<Dimension, Scalar>& specificThermalEnergy,
           const Field<Dimension, Scalar>& pressure,
           const Field<Dimension, Scalar>& fluidSoundSpeed,
           const Field<Dimension, SymTensor>& damage) const {
  Field<Dimension, Scalar> mu("shear modulus", density.nodeList());
  this->shearModulus(mu, density, specificThermalEnergy, pressure, damage);
  const auto n = density.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    CHECK(density(i) > 0.0);
    const auto cs2 = fluidSoundSpeed(i)*fluidSoundSpeed(i) + std::abs(4.0/3.0 * mu(i)) / density(i);
    CHECK(cs2 > 0.0);
    soundSpeed(i) = std::sqrt(cs2); // * std::max(0.0, 1.0 - damage(i).eigenValues().maxElement());
  }
}

//------------------------------------------------------------------------------
// Compute the melt specific energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
meltSpecificEnergy(Field<Dimension, Scalar>& meltSpecificEnergy,
                   const Field<Dimension, Scalar>& density,
                   const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  const auto rho0 = mEOSPtr->referenceDensity();
  CHECK(rho0 > 0.0);
  const auto n = meltSpecificEnergy.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto mu = mEOSPtr->boundedEta(density(i)) - 1.0;
    CHECK(mu >= -1.0);
    meltSpecificEnergy(i) = mMeltEnergyFit(mu)/rho0;   // Note converted to specific energy.
  }
}

//------------------------------------------------------------------------------
// Compute the cold specific energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
coldSpecificEnergy(Field<Dimension, Scalar>& coldSpecificEnergy,
                   const Field<Dimension, Scalar>& density,
                   const Field<Dimension, Scalar>& /*specificEnergy*/) const {
  const auto rho0 = mEOSPtr->referenceDensity();
  CHECK(rho0 > 0.0);
  const auto n = coldSpecificEnergy.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto mu = mEOSPtr->boundedEta(density(i)) - 1.0;
    CHECK(mu >= -1.0);
    coldSpecificEnergy(i) = mColdEnergyFit(mu)/rho0;   // Note converted to specific energy.
  }
}

//------------------------------------------------------------------------------
// Compute the melt attenuation.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
meltAttenuation(const double density, const double specificThermalEnergy) const {
  //const double tiny = 1.0e-5;
  const double rho0 = mEOSPtr->referenceDensity();
  const double mu = mEOSPtr->boundedEta(density) - 1.0;
  CHECK(rho0 > 0.0);
  CHECK(mu >= -1.0);
  const double emelt = mMeltEnergyFit(mu)/rho0;   // Note converted to specific thermal energy.

  double result;
  if (fuzzyGreaterThanOrEqual(specificThermalEnergy, emelt)) {
    result = 0.0;
  } else {
    CHECK(distinctlyLessThan(specificThermalEnergy, emelt));
    result = exp(-mYp*specificThermalEnergy/(emelt - specificThermalEnergy));
  }

  return result;
}

//------------------------------------------------------------------------------
// Compute the Steinberg-Guinan effective temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
computeTemperature(Field<Dimension, Scalar>& temperature,
                   const Field<Dimension, Scalar>& density,
                   const Field<Dimension, Scalar>& specificThermalEnergy) const {
  Field<Dimension, Scalar> eps1("new energy", density.nodeList());
  const auto rho0 = mEOSPtr->referenceDensity();
  const auto n = density.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto mu = mEOSPtr->boundedEta(density(i)) - 1.0;
    CHECK(mu >= -1.0);
    eps1(i) = specificThermalEnergy(i) - mColdEnergyFit(mu)/rho0;   // Note converted to specific thermal energy.
  }
  mEOSPtr->setTemperature(temperature, density, eps1);
  temperature -= 300.0;

  // Field<Dimension, Scalar> cV("specific heat", density.nodeList());
  // mEOSPtr->setSpecificHeat(cV, density, specificThermalEnergy);
  // const double rho0 = mEOSPtr->referenceDensity();
  // CHECK(rho0 > 0.0);
  // for (unsigned i = 0; i != density.numInternalElements(); ++i) {
  //   CHECK(cV(i) > 0.0);
  //   const double emelt = mMeltEnergyFit(mu)/rho0;
  //   const double eps = max(0.0, min(emelt, specificThermalEnergy(i)));
  //   temperature(i) = max(0.0, eps - mColdEnergyFit(mu))*rho0/cV(i) + mRefTempOffset;
  //   temperature(i) = max(0.0, eps - mColdEnergyFit(mu) + mRefEpsOffset)*rho0/cV(i);
  // }
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
G0() const {
  return mG0;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Gmax() const {
  return mGmax;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Y0() const {
  return mY0;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Ymax() const {
  return mYmax;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Yp() const {
  return mYp;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
beta() const {
  return mbeta;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
gamma0() const {
  return mgamma0;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
nhard() const {
  return mnhard;
}

template<typename Dimension>
const NinthOrderPolynomialFit&
SteinbergGuinanStrength<Dimension>::
coldEnergyFit() const {
  return mColdEnergyFit;
}

template<typename Dimension>
const NinthOrderPolynomialFit&
SteinbergGuinanStrength<Dimension>::
meltEnergyFit() const {
  return mMeltEnergyFit;
}

}
