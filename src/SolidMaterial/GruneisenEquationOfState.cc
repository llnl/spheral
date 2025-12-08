//---------------------------------Spheral++----------------------------------//
// GruneisenEquationOfState -- Gruneisen  equation of state.
//
// Created by MLF & JMO Wed Sep 17 13:32:54 EDT 2003
// Reference: Equation of State and Strength of Properties of Selected Materials
//            Daniel J. Steinberg, UCRL-MA-106439, February 13, 1991
//----------------------------------------------------------------------------//

#include "GruneisenEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given Gruneisen constants.
//------------------------------------------------------------------------------
template<typename Dimension>
GruneisenEquationOfState<Dimension>::
GruneisenEquationOfState(const double referenceDensity,
                         const double etamin,
                         const double etamax,
                         const double C0, 
                         const double S1,
                         const double S2,
                         const double S3,
                         const double gamma0,
                         const double b,
                         const double atomicWeight,
                         const PhysicalConstants& constants,
                         const double externalPressure,
                         const double minimumPressure,
                         const double maximumPressure,
                         const double minimumPressureDamage,
                         const MaterialPressureMinType minPressureType):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  min(etamax, max(0.99*S1/(S1 - 1.0), 2.0)),
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minimumPressureDamage,
                                  minPressureType,
                                  externalPressure),
  mC0(C0),
  mS1(S1),
  mS2(S2),
  mS3(S3),
  mgamma0(gamma0),
  mb(b),
  mAtomicWeight(atomicWeight),
  mCv(0.0),
  mEnergyMultiplier(1.0) {
  REQUIRE(distinctlyGreaterThan(mAtomicWeight, 0.0));
//   mCv = 3.0 * 1000.0*Constants::ElectronCharge*Constants::NAvogadro / mAtomicWeight;
  mCv = 3.0 * constants.molarGasConstant() / mAtomicWeight;
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GruneisenEquationOfState<Dimension>::~GruneisenEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GruneisenEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
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
GruneisenEquationOfState<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
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
GruneisenEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    temperature(i) = this->temperature(massDensity(i),specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GruneisenEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
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
GruneisenEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  CHECK(valid());
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GruneisenEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i),specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
GruneisenEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    gamma(i) = this->gamma(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
GruneisenEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    bulkModulus(i) = this->bulkModulus(massDensity(i),specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GruneisenEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    entropy(i) = std::get<0>(pressureAndDerivs(massDensity(i), specificThermalEnergy(i)))*safeInvVar(pow(massDensity(i), gamma(massDensity(i), specificThermalEnergy(i))));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure with partial derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
std::tuple<typename Dimension::Scalar, typename Dimension::Scalar, typename Dimension::Scalar>
GruneisenEquationOfState<Dimension>::
pressureAndDerivs(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const double tiny = 1.0e-20;
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double K0 = rho0*mC0*mC0;
  const double eps = mEnergyMultiplier*specificThermalEnergy;

  //TODO double check branching and apply eta convention in appropriate branch
  if (mu <= 0.0 or eps < 0.0) {
    return std::make_tuple(this->applyPressureLimits(K0*mu + mgamma0*rho0*eps),    // P
                           mgamma0*rho0,                                           // \partial P/\partial eps
                           K0);                                                    // \partial P/\partial rho

  } else {
    CHECK(eta >= -1.0);
    const double etaInv = safeInvVar(eta, tiny);
    const double ack = mu*etaInv;
    const double thpt1 = mu*mu*etaInv;
    const double thpt2 = thpt1*mu*etaInv;
    const double A = K0*mu*(eta - 0.5*mgamma0*mu - 0.5*mb*mu*mu);
    const double dAdrho = mC0*mC0*(eta + (1.0 - mgamma0)*mu - 1.5*mb*mu*mu);
    const double D = 1.0 - (mS1 - 1.0)*mu - mS2*thpt1 - mS3*thpt2;
    const double Dinv = safeInvVar(D, tiny);
    const double dDdrho = (1.0 - mS1 - mS2*ack*(2.0 - ack) - mS3*ack*ack*(3.0 - 2.0*ack))/rho0;
    return std::make_tuple(this->applyPressureLimits(A*Dinv*Dinv + (mgamma0 + mb*mu)*rho0*eps), // P
                           (mgamma0 + mb*mu)*rho0,                                              // \partial P/\partial eps
                           (D*dAdrho - 2.0*A*dDdrho)*FastMath::pow3(Dinv) + mb*eps);            // \partial P/\partial rho
    // return this->applyPressureLimits((K0*mu*(1.0 + (1.0 - 0.5*mgamma0)*mu - 0.5*mb*mu*mu)*Dinv*Dinv + 
    //                                   (mgamma0 + mb*mu)*rho0*eps));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// Need a real temperature relation here
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const double xmu = this->boundedEta(massDensity) - 1.0; 
  CHECK(xmu>-1.);
  const double gamma1=(mgamma0 + mb*xmu) / (1. + xmu) - 1.; 
  const double kB = mConstants.kB();
  const double mp = mConstants.protonMass();
  return gamma1*mAtomicWeight*mp/kB*specificThermalEnergy;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
// Need a real energy relation here
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  CHECK(valid());
  const double xmu = this->boundedEta(massDensity) - 1.; 
  CHECK(xmu>-1.);
  const double gamma1=(mgamma0 + mb*xmu) / (1. + xmu) - 1.; 
  const double kB = mConstants.kB();
  return kB/(gamma1)*temperature;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  CHECK(valid());
  return mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  const double dPdrho = computeDPDrho(massDensity, specificThermalEnergy);
  CHECK(dPdrho >= 0.0);
  return sqrt(dPdrho);
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::gamma(const Scalar massDensity,
                                           const Scalar /*specificThermalEnergy*/) const {
  const double xmu = this->boundedEta(massDensity) - 1.;
  CHECK(xmu!=-1.);
  return (mgamma0 + mb*xmu) / (1. + xmu);
}

//------------------------------------------------------------------------------
// Get the bulk modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {

  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return std::get<0>(this->pressureAndDerivs(massDensity, specificThermalEnergy))*
    safeInvVar(pow(massDensity, gamma(massDensity, specificThermalEnergy)));
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho).
// 
// This turns out to be 
// \partial P       \partial P   |          P     \partial P   |
// ------------   = -------------|      + ------  -------------|
// \partial \rho    \partial \rho|_\eps   \rho^2  \partial \eps|_\rho
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GruneisenEquationOfState<Dimension>::
computeDPDrho(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
  const double eta = this->boundedEta(massDensity),
               rho0 = this->referenceDensity(),
               rho = rho0*eta;
  const auto [P, dPdeps, dPdrho] = this->pressureAndDerivs(massDensity, specificThermalEnergy);  // Note, requires C++17
  const auto dPdrho_ad = dPdrho + dPdeps*P/(rho*rho);
  return std::abs(dPdrho_ad);

  // CHECK(valid());
  // const double tiny = 1.0e-20;
  // const double eta = this->boundedEta(massDensity);
  // const double mu = eta - 1.0;
  // const double rho0 = this->referenceDensity();
  // const double rho = rho0*eta;
  // const double eps = mEnergyMultiplier*specificThermalEnergy;

  // double ack;
  // if (mu <= 0.0) {
  //   ack = 1.0;
  // } else {
  //   const double N = (1.0 + (1.0 - 0.5*mgamma0 - 0.5*mb*mu)*mu)*mu;
  //   const double x = mu/(1.0 + mu);
  //   const double D = 1.0 - (mS1 - 1.0 + (mS2 + mS3*x)*x)*mu;
  //   const double dNdmu = 1.0 + (2.0 - mgamma0 - 1.5*mb*mu)*mu;
  //   const double dDdmu = - (mS1 - 1.0 + (mS2*(mu + 2.0) + mS3*(mu + 3.0)*x)*x/(1.0 + mu));
  //   const double Dinv = 1.0/(sgn(D)*max(abs(D), tiny));
  //   ack = max(0.0, (dNdmu*D - 2.0*N*dDdmu)*FastMath::cube(Dinv));
  // }
  // CHECK(ack >= 0.0);
  // const double dpdrho_cold = mC0*mC0*ack;

  // // Put the whole thing together, depending on the thermal energy.
  // if (mu <= 0.0 or eps < 0.0) {
  //   return dpdrho_cold;
  // } else {
  //   const double Prho2 = std::get<0>(this->pressureAndDerivs(massDensity, specificThermalEnergy))/(rho*rho);
  //   return std::max(0.0, dpdrho_cold + max(0.0, mb*eps + mb*Prho2));
  // }
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GruneisenEquationOfState<Dimension>::valid() const {
  return (mC0 > 0.0 &&
          mS1 > 0.0 &&
          mCv > 0.0);
}

}

