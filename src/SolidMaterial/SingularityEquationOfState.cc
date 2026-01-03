//---------------------------------Spheral++----------------------------------//
// SingularityEquationOfState
//
// The base class for Singularity EOS
//----------------------------------------------------------------------------//
#include "SingularityEquationOfState.hh"

#include <iomanip>
#include <iostream>

#include "Distributed/Process.hh"
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
                           const SingularityLimitLevel level,
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
  mLevel(level),
  mCGS(0.01, 0.001, 1.0, 1.0),
  mToCGS(constants, mCGS),
  mFromCGS(mCGS, constants),
   // For bulk modulus to be finite, cs < sqrt(limit / rho) with rho=1000 in CGS here
  mCsMin(0.0),
  mCsMax(std::sqrt(std::numeric_limits<double>::max() / 1000)),
  mRhoMin(1e-12),
  mRhoMax(1000) {
  // Set up grid of densities for calculating non-garbage results
  if (mLevel == SingularityLimitLevel::NONE) {
    // Only need the first element
    mDensitiesCGS = {0.0};
    mMinEnergiesCGS = {std::numeric_limits<double>::lowest()};
    mMinTemperaturesCGS = {std::numeric_limits<double>::lowest()};
  }
  else {
    const auto rhoSmallMin = mRhoMin, rhoMin = 0.2, rhoMax = 100.0, rhoBigMax = 200.0;
    const auto multRhoSmall = 10.0, dRho = 0.2, dRhoBig = 5.0;
    mDensitiesCGS.push_back(0.0); // 0
    for (auto rho = rhoSmallMin; rho <= rhoMin; rho *= multRhoSmall) {
      // 1.0e-N with N=-12..-1
      mDensitiesCGS.push_back(rho);
    }
    for (auto rho = rhoMin; rho <= rhoMax; rho += dRho) {
      // 0.2 to 100 in steps of 0.2
      mDensitiesCGS.push_back(rho);
    }
    for (auto rho = rhoMax + dRhoBig; rho <= rhoBigMax; rho += dRhoBig) {
      // 105 to 200 in steps of 5
      mDensitiesCGS.push_back(rho);
    }
    // Absurd limit on top
    mDensitiesCGS.push_back(mRhoMax);

    // Bounds from 0 to 100 keV
    const auto tMax = 100 / (8.61733e-8/*kB (keV/K)*/); // 100 keV into K
    const auto eMin = 1e-10; // erg/g
    const auto eMax = mSEOS.InternalEnergyFromDensityTemperature(2.0, // density of 1.0 for reference
                                                                 tMax);

    // Set up bounds for calculation
    const auto step = 1.1;
    const auto eMult = 2.0; // multiply result by this to be safe
    mCsMin = 1; // cm/s
    mCsMax = 0.01 * mCGS.c(); // cm/s

    // Calculate min energy and temperature per density that doesn't give garbage sound speed results
    const auto numDensities = mDensitiesCGS.size();
    auto avE = 0.0;
    auto avT = 0.0;
    auto minE = std::numeric_limits<double>::max();
    auto minT = std::numeric_limits<double>::max();
    auto maxE = 0.0;
    auto maxT = 0.0;
    mMinEnergiesCGS.resize(numDensities);
    mMinTemperaturesCGS.resize(numDensities);
    for (auto i = 0u; i < numDensities; ++i) {
      const auto rho = mDensitiesCGS[i];
      auto found = false;
      auto cs = 2.0 * mCGS.c();
      auto e = eMin;

      // Stop iterations once max energy is reached
      while (e < eMax) {
        const auto bm = mSEOS.BulkModulusFromDensityInternalEnergy(rho, e);
        cs = std::sqrt(bm / rho);

        // 1. Sound speed should be below limit
        // 2. Sound speed should be going up as energy increases, not down
        if (std::isfinite(cs) and cs < mCsMax) {
          found = true;
          break;
        }
      
        e *= step;
      }

      mMinEnergiesCGS[i] = found ? eMult * e : eMin;
      mMinTemperaturesCGS[i] = mSEOS.TemperatureFromDensityInternalEnergy(rho, mMinEnergiesCGS[i]);
      avE += mMinEnergiesCGS[i];
      avT += mMinTemperaturesCGS[i];
      minE = std::min(minE, mMinEnergiesCGS[i]);
      minT = std::min(minT, mMinTemperaturesCGS[i]);
      maxE = std::max(maxE, mMinEnergiesCGS[i]);
      maxT = std::max(maxT, mMinTemperaturesCGS[i]);
    }
    avE /= numDensities;
    avT /= numDensities;

    // Print out the values
    if (Process::getRank() == 0) {
      std::cout << std::scientific << std::setprecision(2)
                << "Singularity sound speed garbage detection cutoff (min,max,av):"
                << std::endl
                << "\te=(" << minE << "," << maxE << "," << avE << ") "
                << "\tT=(" << minT << "," << maxT << "," << avT << ")"
                << std::endl;
    }

    CHECK(mCsMax <= mCGS.c());
  }
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
  const auto specificThermalEnergyCGS = eLim(massDensityCGS,
                                             mToCGS.specificEnergy() * specificThermalEnergy,
                                             SingularityLimitLevel::ALL);
  const auto valCGS = mSEOS.PressureFromDensityInternalEnergy(massDensityCGS,
                                                              specificThermalEnergyCGS);
  const auto val = mFromCGS.pressure() * valCGS;
  CHECK(std::isfinite(val));
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
  const auto specificThermalEnergyCGS = eLim(massDensityCGS,
                                             mToCGS.specificEnergy() * specificThermalEnergy,
                                             SingularityLimitLevel::ALL);
  const auto valCGS = mSEOS.TemperatureFromDensityInternalEnergy(massDensityCGS,
                                                                 specificThermalEnergyCGS);
  const auto val =  mFromCGS.temperature() * valCGS;
  CHECK(std::isfinite(val));
  CHECK(val >= 0);
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
  const auto temperatureCGS = tLim(massDensityCGS,
                                   mToCGS.temperature() * temperature,
                                   SingularityLimitLevel::ALL);
  const auto valCGS = mSEOS.InternalEnergyFromDensityTemperature(massDensityCGS,
                                                                 temperatureCGS);
  const auto val =  mFromCGS.specificEnergy() * valCGS;
  CHECK(std::isfinite(val));
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
  const auto temperatureCGS = tLim(massDensityCGS,
                                   mToCGS.temperature() * temperature,
                                   SingularityLimitLevel::ALL);
  const auto valCGS = mSEOS.SpecificHeatFromDensityTemperature(massDensityCGS,
                                                               temperatureCGS);
  const auto val =  mFromCGS.specificHeat() * valCGS;
  CHECK(std::isfinite(val));
  CHECK(val >= 0);
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
  const auto specificThermalEnergyCGS = eLim(massDensityCGS,
                                             mToCGS.specificEnergy() * specificThermalEnergy,
                                             SingularityLimitLevel::ALL);
  const auto val =  1.0 +
    mSEOS.GruneisenParamFromDensityInternalEnergy(massDensityCGS,
                                                  specificThermalEnergyCGS);
  CHECK(std::isfinite(val));
  CHECK(val >= 0);
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
  const auto specificThermalEnergyCGS = eLim(massDensityCGS,
                                             mToCGS.specificEnergy() * specificThermalEnergy,
                                             SingularityLimitLevel::CS);
  const auto valCGS = mSEOS.BulkModulusFromDensityInternalEnergy(massDensityCGS,
                                                                 specificThermalEnergyCGS);
  const auto val =  mFromCGS.bulkModulus() * bmLim(massDensityCGS,
                                                   valCGS,
                                                   SingularityLimitLevel::CS);
  CHECK(std::isfinite(val));
  CHECK(val >= 0);
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
