//---------------------------------Spheral++----------------------------------//
// SingularityEquationOfState
//
// The base class for Singularity EOS
//----------------------------------------------------------------------------//
#include "SingularityEquationOfState.hh"

#include <iostream>

#include "Field/Field.hh"
#include "Material/PhysicalConstants.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SingularityEquationOfState<Dimension>::
SingularityEquationOfState(SingularityEOSType eos,
                           const PhysicalConstants& constants,
                           const double minDensityCGS,
                           const double referenceDensity,
                           const double etamin,
                           const double etamax,
                           const double minimumPressure,
                           const double maximumPressure,
                           const double minimumPressureDamage,
                           const MaterialPressureMinType minPressureType,
                           const double externalPressure):
  SolidEquationOfState<Dimension>(referenceDensity, etamin, etamax, constants,
                                  minimumPressure, maximumPressure, minimumPressureDamage,
                                  minPressureType, externalPressure),
  mSEOS(eos),
  mCGS(0.01, 0.001, 1.0, 1.0),
  mToCGS(constants, mCGS),
  mFromCGS(mCGS, constants),
  mCsMin(1),                                                     // cm/s
  mCsMax(0.01 * mCGS.c()),                                       // cm/s
  mCvMin(1e-10),                                                  // erg/(g*K)
  mCvMax(1e15),                                                   // erg/(g*K)
  mRhoMin(minDensityCGS),                                         // g/cc
  mRhoMax(1000) {                                                  // g/cc
  CHECK2(mCsMax <= mCGS.c(), "mCsMax=" << mCsMax << " c=" << mCGS.c());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SingularityEquationOfState<Dimension>::
~SingularityEquationOfState() {
  mSEOS.Finalize();
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& pressure,
                     Field<Dimension, Scalar>& dPdu,
                     Field<Dimension, Scalar>& dPdrho,
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
    dPdu(i) = this->gamma(massDensity(i), specificThermalEnergy(i)) * massDensity(i);
    dPdrho(i) = this->bulkModulus(massDensity(i), specificThermalEnergy(i)) * safeInvVar(massDensity(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    specificHeat(i) = this->specificHeat(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    gamma(i) = this->gamma(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho). 
// This is just the specific heat ratio times pressure for a gamma law gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    bulkModulus(i) = this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SingularityEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (auto i = 0u; i < massDensity.numElements(); ++i) {
    entropy(i) = this->entropy(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  const auto specificThermalEnergyCGS = mToCGS.specificEnergy() * specificThermalEnergy;
  const auto valCGS = mSEOS.PressureFromDensityInternalEnergy(massDensityCGS,
                                                              specificThermalEnergyCGS);
  const auto val = mFromCGS.pressure() * valCGS;
  CHECK2(std::isfinite(val), "P=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  return this->applyPressureLimits(val);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  const auto specificThermalEnergyCGS = mToCGS.specificEnergy() * specificThermalEnergy;
  const auto valCGS = mSEOS.TemperatureFromDensityInternalEnergy(massDensityCGS,
                                                                 specificThermalEnergyCGS);
  const auto val =  mFromCGS.temperature() * valCGS;
  CHECK2(std::isfinite(val), "T=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  CHECK2(val >= 0, "T=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  return val;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  CHECK(valid());
  const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  const auto temperatureCGS = mToCGS.temperature() * temperature;
  const auto valCGS = mSEOS.InternalEnergyFromDensityTemperature(massDensityCGS,
                                                                 temperatureCGS);
  const auto val =  mFromCGS.specificEnergy() * valCGS;
  CHECK2(std::isfinite(val), "e=" << val << ", rho=" << massDensity << ", T=" << temperature);
  return val;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  CHECK(valid());
  const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  const auto temperatureCGS = mToCGS.temperature() * temperature;
  const auto valCGS = cvLim(mSEOS.SpecificHeatFromDensityTemperature(massDensityCGS,
                                                                     temperatureCGS));
  const auto val =  mFromCGS.specificHeat() * valCGS;
  CHECK2(std::isfinite(val), "cv=" << val << ", rho=" << massDensity << ", T=" << temperature);
  CHECK2(val >= 0, "cv=" << val << ", rho=" << massDensity << ", T=" << temperature);
  return val;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const auto massDensityLim = mFromCGS.massDensity() * rhoLim(mToCGS.massDensity() * massDensity);
  const auto val =  std::sqrt(this->bulkModulus(massDensity, specificThermalEnergy) *
                              safeInvVar(massDensityLim));
  // Don't need to limit sound speed or convert since that's done in bulkModulus
  CHECK2(std::isfinite(val) && val >= 0 && val <= this->constants().c(),
         "sound speed outside limits: cs=" << val << ", c=" << this->constants().c() << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  return val;
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  const auto specificThermalEnergyCGS = mToCGS.specificEnergy() * specificThermalEnergy;
  const auto val =  1.0 +
    mSEOS.GruneisenParamFromDensityInternalEnergy(massDensityCGS,
                                                  specificThermalEnergyCGS);
  CHECK2(std::isfinite(val), "gamma=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  CHECK2(val >= 0, "gamma=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  return val;
}

//------------------------------------------------------------------------------
// Calculate an individual bulk modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  const auto specificThermalEnergyCGS = mToCGS.specificEnergy() * specificThermalEnergy;
  const auto valCGS = mSEOS.BulkModulusFromDensityInternalEnergy(massDensityCGS,
                                                                 specificThermalEnergyCGS);
  const auto val =  mFromCGS.bulkModulus() * bmLim(massDensityCGS, valCGS);
  CHECK2(std::isfinite(val), "bm=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  CHECK2(val >= 0, "bm=" << val << ", rho=" << massDensity << ", e=" << specificThermalEnergy);
  return val;
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SingularityEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  return 0.0;
  // Not available for many EOS options, so just comment it out for now
  // CHECK(valid());
  // const auto massDensityCGS = rhoLim(mToCGS.massDensity() * massDensity);
  // const auto specificThermalEnergyCGS = eLim(massDensityCGS,
  //                                            mToCGS.specificEnergy() * specificThermalEnergy,
  //                                            SingularityLimitLevel::ALL);
  // const auto valCGS = mSEOS.EntropyFromDensityInternalEnergy(massDensityCGS,
  //                                                            specificThermalEnergyCGS);
  // const auto val =  mFromCGS.entropy() * valCGS;
  // CHECK(std::isfinite(val) && val >= 0);
  // return val;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SingularityEquationOfState<Dimension>::valid() const {
  return true;
}

} // namespace Spheral
