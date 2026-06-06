//---------------------------------Spheral++----------------------------------//
// OsborneEquationOfState -- Osborne  equation of state.
//
// Reference: PAGOSA Physics manual, LA-14425-M
//----------------------------------------------------------------------------//

#include "OsborneEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given Osborne constants.
//------------------------------------------------------------------------------
template<typename Dimension>
OsborneEquationOfState<Dimension>::
OsborneEquationOfState(const double referenceDensity,
                       const double etamin,
                       const double etamax,
                       const double a1,
                       const double a2pos,
                       const double a2neg,
                       const double b0,
                       const double b1,
                       const double b2pos,
                       const double b2neg,
                       const double c0,
                       const double c1,
                       const double c2pos,
                       const double c2neg,
                       const double E0,
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
  mA1(a1),
  mA2pos(a2pos),
  mA2neg(a2neg),
  mB0(b0),
  mB1(b1),
  mB2pos(b2pos),
  mB2neg(b2neg),
  mC0(c0),
  mC1(c1),
  mC2pos(c2pos),
  mC2neg(c2neg),
  mE0(E0),
  mAtomicWeight(atomicWeight),
  mCv(3.0 * constants.molarGasConstant() / mAtomicWeight) {
  VERIFY(distinctlyGreaterThan(mAtomicWeight, 0.0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
OsborneEquationOfState<Dimension>::~OsborneEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const double rho0 = this->referenceDensity();
  const auto n = pressure.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const double eta = this->boundedEta(massDensity(i));
    const double mu = eta - 1.0;
    const double E = rho0*specificThermalEnergy(i);
    const double a2 = mu > 0.0 ? mA2pos : mA2neg;
    const double b2 = mu > 0.0 ? mB2pos : mB2neg;
    const double c2 = mu > 0.0 ? mC2pos : mC2neg;
    pressure(i) = this->applyPressureLimits((mA1*mu + a2*mu*mu +
                                             (mB0 + mB1*mu + b2*mu*mu)*E +
                                             (mC0 + mC1*mu + c2*mu*mu)*E*E)*safeInvVar(E + mE0));
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {

  CHECK(valid());
  const double rho0 = this->referenceDensity();
  const auto n = Pressure.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const double eta = this->boundedEta(massDensity(i));
    const double mu = eta - 1.0;
    const double E = rho0*specificThermalEnergy(i);
    const double a2 = mu > 0.0 ? mA2pos : mA2neg;
    const double b2 = mu > 0.0 ? mB2pos : mB2neg;
    const double c2 = mu > 0.0 ? mC2pos : mC2neg;
    Pressure(i) = this->applyPressureLimits((mA1*mu + a2*mu*mu +
                                             (mB0 + mB1*mu + b2*mu*mu)*E +
                                             (mC0 + mC1*mu + c2*mu*mu)*E*E)*safeInvVar(E + mE0));
    dPdu(i) = rho0*(E + mE0)*(mB0 + mB1*mu + b2*mu*mu + 2.0*E*(mC0 + mC1*mu + c2*mu*mu) - Pressure(i))*safeInvVar(FastMath::pow2(E + mE0));
    dPdrho(i) = (mA1 + 2.0*a2*mu + (mB1 + 2.0*b2*mu)*E + (mC1 + 2.0*c2*mu)*E*E)*safeInvVar(rho0*(E + mE0));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& /*temperature*/,
               const Field<Dimension, Scalar>& /*massDensity*/,
               const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  VERIFY2(false, "temperature unimplemented for the Osborne equation of state.");
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& /*specificThermalEnergy*/,
                         const Field<Dimension, Scalar>& /*massDensity*/,
                         const Field<Dimension, Scalar>& /*temperature*/) const {
  VERIFY2(false, "specific thermal energy unimplemented for the Osborne equation of state.");
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
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
OsborneEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = soundSpeed.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const double dPdrhoi = DPDrho(massDensity(i), specificThermalEnergy(i));
    CHECK(dPdrhoi >= 0.0);
    soundSpeed(i) = sqrt(dPdrhoi);
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  CHECK(mCv > 0.0);
  const auto n = gamma.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const double eta = this->boundedEta(massDensity(i)),
                 rho0 = this->referenceDensity(),
                 rho = rho0*eta,
                 nDen = rho/mAtomicWeight;
    gamma(i) = 1.0 + mConstants.molarGasConstant()*nDen/mCv;
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const auto n = bulkModulus.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    bulkModulus(i) = massDensity(i) * DPDrho(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  this->setPressure(entropy, massDensity, specificThermalEnergy);
  const auto n = massDensity.numElements();
#pragma omp parallel for
  for (auto i = 0u; i < n ; ++i) {
    const double eta = this->boundedEta(massDensity(i)),
                 rho0 = this->referenceDensity(),
                 rho = rho0*eta,
                 nDen = rho/mAtomicWeight,
                 gamma = 1.0 + mConstants.molarGasConstant()*nDen/mCv;
    entropy(i) *= safeInvVar(pow(massDensity(i), gamma));
  }
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
double
OsborneEquationOfState<Dimension>::
DPDrho(const double massDensity,
       const double specificThermalEnergy) const {
  //const double tiny = 1.0e-20;
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double rho = rho0*eta;
  const double E = rho0*specificThermalEnergy;
  const double a2 = mu > 0.0 ? mA2pos : mA2neg;
  const double b2 = mu > 0.0 ? mB2pos : mB2neg;
  const double c2 = mu > 0.0 ? mC2pos : mC2neg;
  const double P = (mA1*mu + a2*mu*mu +
                    (mB0 + mB1*mu + b2*mu*mu)*E +
                    (mC0 + mC1*mu + c2*mu*mu)*E*E)*safeInvVar(E + mE0);
  const double DPDrho_e = (mA1 + 2.0*a2*mu + E*(mB1 + 2.0*b2*mu + mC1*E))*safeInvVar(rho0*(E + mE0));
  const double DPDeps_r = rho0*(mB0 + mB1*mu + b2*mu*mu + 2.0*E*(mC0 + mC1*mu) - P)*safeInv(E + mE0);
  return max(0.0, DPDrho_e + P/(rho*rho)*DPDeps_r);
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
OsborneEquationOfState<Dimension>::valid() const {
  return true;
}

}
