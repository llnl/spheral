//---------------------------------Spheral++----------------------------------//
// LinearPolynomialEquationOfState -- An equation of state approximated by a
// linear polynomial, i.e.:
//
//   P(rho, e) = A0 + A1*mu + a2*mu^2 + a3*mu^3 + (B0 + B1*mu + B2*mu^2)*e
//   mu = rho/rho0 - 1.0
//
// Created by JMO, Thu May  5 16:07:36 PDT 2005
//----------------------------------------------------------------------------//
#include "LinearPolynomialEquationOfState.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
LinearPolynomialEquationOfState<Dimension>::
LinearPolynomialEquationOfState(const double referenceDensity,
                                const double etamin,
                                const double etamax,
                                const double a0,
                                const double a1,
                                const double a2,
                                const double a3,
                                const double b0,
                                const double b1,
                                const double b2,
                                const double atomicWeight,
                                const PhysicalConstants& constants,
                                const double externalPressure,
                                const double minimumPressure,
                                const double maximumPressure,
                                const double minimumPressureDamage,
                                const MaterialPressureMinType minPressureType):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minimumPressureDamage,
                                  minPressureType,
                                  externalPressure),
  mA0(a0),
  mA1(a1),
  mA2(a2),
  mA3(a3),
  mB0(b0),
  mB1(b1),
  mB2(b2),
  mAtomicWeight(atomicWeight),
  mCv(3.0 * constants.molarGasConstant() / atomicWeight),
  mGamma(mB0 + 1.0) {
  REQUIRE(distinctlyGreaterThan(mAtomicWeight, 0.0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LinearPolynomialEquationOfState<Dimension>::
~LinearPolynomialEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    pressure(i) = std::get<0>(this->pressureAndDerivs(massDensity(i), specificThermalEnergy(i)));
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto stuff = this->pressureAndDerivs(massDensity(i), specificThermalEnergy(i));
    pressure(i) = std::get<0>(stuff);
    dPdu(i) = std::get<1>(stuff);
    dPdrho(i) = std::get<2>(stuff);
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  REQUIRE(valid());
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& /*massDensity*/,
	      const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a 
// polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    bulkModulus(i)=this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearPolynomialEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    entropy(i) = std::get<0>(pressureAndDerivs(massDensity(i), specificThermalEnergy(i)))*safeInvVar(pow(massDensity(i), mGamma));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
std::tuple<typename Dimension::Scalar, typename Dimension::Scalar, typename Dimension::Scalar>
LinearPolynomialEquationOfState<Dimension>::
pressureAndDerivs(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  if (fuzzyEqual(eta, this->etamin())) return std::make_tuple(0.0, 0.0, 0.0);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  return std::make_tuple(this->applyPressureLimits(mA0 + mA1*mu + mA2*mu*mu + mA3*mu*mu*mu +                   // P
                                                   (mB0 + mB1*mu + mB2*mu*mu)*specificThermalEnergy),
                         mB0 + mB1*mu + mB2*mu*mu,                                                             // \partial P/\partial u
                         (mA1 + 2.0*mA2*mu + 3.0*mA3*mu*mu + (mB1 + 2.0*mB2*mu)*specificThermalEnergy)/rho0);  // \partial P/\partial rho
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is a *hokey* definition -- have to do better if we ever really care
// about the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return specificThermalEnergy/mCv + 300;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar temperature) const {
  REQUIRE(valid());
  return (temperature - 300.0)*mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  REQUIRE(valid());
  return mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double c2 = computeDPDrho(massDensity, specificThermalEnergy);
  ENSURE(c2 >= 0.0);
  return sqrt(c2);
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
gamma(const Scalar /*massDensity*/,
      const Scalar /*specificThermalEnergy*/) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return std::get<0>(this->pressureAndDerivs(massDensity, specificThermalEnergy))*safeInvVar(pow(massDensity, mGamma));
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho) at constant entropy.  This is the form
// used for the sound speed.
// 
// This turns out to be 
// \partial P       \partial P   |          P     \partial P   |
// ------------   = -------------|      + ------  -------------|
// \partial \rho    \partial \rho|_\eps   \rho^2  \partial \eps|_\rho
//------------------------------------------------------------------------------
template<typename Dimension>
double
LinearPolynomialEquationOfState<Dimension>::
computeDPDrho(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double rho = rho0*eta;
  const double dPdrho_eps = std::abs(mA1 + 2.0*mA2*mu + 3.0*mA3*mu*mu +
                                     (mB1 + 2.0*mB2*mu)*specificThermalEnergy)/rho0;
  const double Prho2 = std::get<0>(this->pressureAndDerivs(massDensity, specificThermalEnergy))/(rho*rho);
  const double dPdeps_rho = mB0 + mB1*mu + mB2*mu*mu;
  const double result = std::max(0.0, dPdrho_eps + Prho2*dPdeps_rho);

  ENSURE(result >= 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
LinearPolynomialEquationOfState<Dimension>::valid() const {
  return (SolidEquationOfState<Dimension>::valid() && 
          mAtomicWeight > 0.0 &&
          mCv > 0.0 &&
          mGamma > 0.0);
}

}
