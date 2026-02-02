//---------------------------------Spheral++----------------------------------//
// TillotsonEquationOfState -- This equation of state is designed to represent
// metallic materials over a range of pressure and density -- spanning solid, 
// liquid, and vapor states.
//
// Tillotson 1962
//
// Created by JMO, Wed Mar 16 23:31:17 PDT 2011
//----------------------------------------------------------------------------//
#include "TillotsonEquationOfState.hh"
#include "Field/Field.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
TillotsonEquationOfState<Dimension>::
TillotsonEquationOfState(const double referenceDensity,
                         const double etamin,
                         const double etamax,
                         const double etamin_solid,
                         const double etamax_solid,
                         const double a,
                         const double b,
                         const double A,
                         const double B,
                         const double alpha,
                         const double beta,
                         const double eps0,
                         const double epsLiquid,
                         const double epsVapor,
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
  mEtaMinSolid(etamin_solid),
  mEtaMaxSolid(etamax_solid),
  ma(a),
  mb(b),
  mA(A),
  mB(B),
  malpha(alpha),
  mbeta(beta),
  meps0(eps0),
  mepsLiquid(epsLiquid),
  mepsVapor(epsVapor),
  mAtomicWeight(atomicWeight),
  mCv(3.0 * constants.molarGasConstant() / atomicWeight),
  mdPdRhoMin(0.0) {
  VERIFY(distinctlyGreaterThan(mAtomicWeight/constants.molarGasConstant(),0.0));
  // Tillotson can compute non-physical negative sound speeds in expansion, so
  // we protect from this with a minimum allowed value (1% of reference sound speed).
  mdPdRhoMin = 1e-4*computeDPDrho(referenceDensity, 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TillotsonEquationOfState<Dimension>::
~TillotsonEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = pressure.size();
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
TillotsonEquationOfState<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = Pressure.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto stuff = this->pressureAndDerivs(massDensity(i), specificThermalEnergy(i));
    Pressure(i) = std::get<0>(stuff);
    dPdu(i) = std::get<1>(stuff);
    dPdrho(i) = std::get<2>(stuff);
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
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
TillotsonEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
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
TillotsonEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
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
TillotsonEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    gamma(i) = this->gamma(massDensity(i),specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a 
// polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
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
TillotsonEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    entropy(i) = this->entropy(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure and it's derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
std::tuple<typename Dimension::Scalar, typename Dimension::Scalar, typename Dimension::Scalar>
TillotsonEquationOfState<Dimension>::
pressureAndDerivs(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const {
  const double eta = this->boundedEta(massDensity),
               mu = eta - 1.0,
               rho0 = this->referenceDensity(),
               rho = rho0*eta,
               eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.

  // Tillotson defines three fundamental pressure regimes:
  //   P1 - solid, compression.
  //   P2 - solid, expansion.
  //   P4 - gaseous, expansion.
  // A forth category (P3) is interpreted as a mixture of gaseous and solid
  // phases, and is interpolated between P2 and P4.

  const auto ack = safeInvVar(1.0 + eps*safeInvVar(meps0*eta*eta, 1.0e-10), 1.0e-10);
  const auto phi = mb*ack;
  const auto dphidu = -phi*ack*safeInvVar(meps0*eta*eta);
  const auto dphidrho = 2.0*phi*ack*eps*safeInvVar(rho0*meps0*eta*eta*eta);

  const auto chi = safeInvVar(eta, 1.0e-10) - 1.0;
  const auto dchidrho = -safeInvVar(rho0*eta*eta);

  // Now determine which regime we're in.
  double P, dPdu, dPdrho;
  if (mu >= 0.0) {

    // Regime 1: compression, solid.
    P = (ma + phi)*rho*eps + mA*mu + mB*mu*mu;
    dPdu = (ma + phi)*rho + rho*eps*dphidu;
    dPdrho = (ma + phi)*eps + rho*eps*dphidrho + (mA + 2.0*mB*mu)*safeInvVar(rho0);

  } else if (eps <= mepsLiquid) {

    // Regime 2: expansion, solid : same as 1, only if rho>cutoff density.
    // Simply set B = 0 as described in Saito et al. 2008
    P = (ma + phi)*rho*eps + mA*mu;
    dPdu = (ma + phi)*rho + rho*eps*dphidu;
    dPdrho = (ma + phi)*eps + rho*eps*dphidrho + mA*safeInvVar(rho0);

  } else if (eps >= mepsVapor) {

    // Regime 4: expansion, vapor.
    const auto expalpha = exp(-malpha*chi*chi);
    const auto expbeta = exp(-mbeta*chi);
    P = ma*rho*eps + (phi*rho*eps + mA*mu*expbeta)*expalpha;
    dPdu = ma*rho + rho*(phi + eps*dphidu)*expalpha;
    dPdrho = ma*eps + (-(phi*rho*eps + mA*mu*expbeta)*2.0*malpha*chi*dchidrho +
                       (phi + rho*dphidrho)*eps + mA*(safeInvVar(rho0) - mu*mbeta*dchidrho)*expbeta)*expalpha;
  
  } else {

    // Regime 3: expansion, partial vapor.
    // Treated here as a linear combination of the solid and gaseous phases.
    // Following <strike>Saito et al. we compute P2 and P4 at the epsLiquid and
    // epsVapor specific energies</strike> Melosh (personal communication) we
    // compute P2 and P4 at the given energy.
    // Regimes 1 and 2 (solid)
    const auto P2 = (ma + phi)*rho*eps + mA*mu;             //  + mB*mu*mu;
    const auto dP2du = (ma + phi)*rho + rho*eps*dphidu;
    const auto dP2drho = (ma + phi)*eps + rho*eps*dphidrho; // + (mA + 2.0*mB*mu)*safeInvVar(rho0);

    // Regime 4 (vapor)
    const auto expalpha = exp(-malpha*chi*chi);
    const auto expbeta = exp(-mbeta*chi);
    const auto P4 = ma*rho*eps + (phi*rho*eps + mA*mu*expbeta)*expalpha;
    const auto dP4du = ma*rho + rho*(phi + eps*dphidu)*expalpha;
    const auto dP4drho = ma*eps - ((phi*rho*eps + mA*mu*expbeta)*2.0*malpha*chi*dchidrho +
                                   (phi*eps* + rho*eps*dphidrho + mA*(safeInvVar(rho0) - mu*mbeta*dchidrho)*expbeta)*expalpha);

    const auto denomInv = safeInvVar(mepsVapor - mepsLiquid);
    const auto feps = (eps - mepsLiquid)*denomInv;
    P = P2 + (P4 - P2)*feps;
    dPdu = dP2du + (P4 - P2)*denomInv + (dP4du - dP2du)*feps;
    dPdrho = dP2drho + (dP4drho - dP2drho)*feps;

  }

  // That's it.  Note, if we run into the pressure limits we also zero out the derivatives.
  const auto Plim = this->applyPressureLimits(P);
  if (P == Plim) {
    return std::make_tuple(Plim, dPdu, dPdrho);
  } else {
    return std::make_tuple(Plim, 0.0, 0.0);
  }
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is a *hokey* definition -- have to do better if we ever really care
// about the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar specificThermalEnergy) const {
  const double eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.
  return eps/mCv + 300;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar temperature) const {
  return (temperature - 300.0)*mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  return mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  const double c2 = computeDPDrho(massDensity, specificThermalEnergy);
  //ENSURE(c2 >= 0.0);
  return sqrt(max(0.0, c2));
}

//------------------------------------------------------------------------------
// Get gamma.
// We assume Cp = Cv for the solid phase, Cp = Cv + R in the gaseous phase,
// and blend the two for liquid.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  // const double eta = this->boundedEta(massDensity),
  //              rho0 = this->referenceDensity(),
  //              rho = rho0*eta,
  //              nDen = rho/mAtomicWeight;
  // CHECK(mCv > 0.0);
  // return 1.0 + mConstants.molarGasConstant()*nDen/mCv;

  const double eta = this->boundedEta(massDensity),
               mu = eta - 1.0,
               //rho0 = this->referenceDensity(),
               //rho = rho0*eta,
               eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.
  if (mu >= 0.0) {

    // Regime 1: compression, solid.
    return 1.0;

  } else if (eps <= mepsLiquid) {

    // Regime 2: expansion, solid : same as 1, but only if rho>cutoff density
    return 1.0;

  } else if (eps >= mepsVapor) {
  
    // Regime 4: expansion, gaseous.
    return 1.0 + mConstants.molarGasConstant()/mCv;

  } else {

    // Regime 3: expansion, liquid.
    // Treated here as a linear combination of the solid and gaseous phases.
    const double gammaVapor = 1.0 + mConstants.molarGasConstant()/mCv;
    return 1.0 + (gammaVapor - 1.0)*(eps - mepsLiquid)/(mepsVapor - mepsLiquid);

  }
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  return std::get<0>(this->pressureAndDerivs(massDensity, specificThermalEnergy))*safeInvVar(pow(massDensity, gamma(massDensity, specificThermalEnergy)));
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho) (adiabatic)
// 
// This turns out to be 
// \partial P       \partial P   |          P     \partial P   |
// ------------   = -------------|      + ------  -------------|
// \partial \rho    \partial \rho|_\eps   \rho^2  \partial \eps|_\rho
//------------------------------------------------------------------------------
template<typename Dimension>
double
TillotsonEquationOfState<Dimension>::
computeDPDrho(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
  const double eta = this->boundedEta(massDensity),
               rho0 = this->referenceDensity(),
               rho = rho0*eta;
  const auto [P, dPdeps, dPdrho] = this->pressureAndDerivs(massDensity, specificThermalEnergy);  // Note, requires C++17
  const auto dPdrho_ad = dPdrho + dPdeps*P/(rho*rho);
  return std::max(dPdrho_ad, mdPdRhoMin);
}

}
