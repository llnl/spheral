//---------------------------------Spheral++----------------------------------//
// GammaLawGas -- The gamma law gas equation of state.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#include "GammaLawGas.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"

#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given gamma and mu.
//------------------------------------------------------------------------------
template<typename Dimension>
GammaLawGas<Dimension>::
GammaLawGas(const double gamma,
            const double mu,
            const PhysicalConstants& constants,
            const double minimumPressure,
            const double maximumPressure,
            const MaterialPressureMinType minPressureType,
            const double externalPressure):
  EquationOfState<Dimension>(constants, minimumPressure, maximumPressure, minPressureType, externalPressure),
  mGamma(gamma),
  mMolecularWeight(mu) {
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GammaLawGas<Dimension>::
~GammaLawGas() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    Pressure(i) = pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    Pressure(i) = pressure(massDensity(i), specificThermalEnergy(i));
    dPdu(i) = mGamma1*massDensity(i);
    dPdrho(i) = mGamma1*specificThermalEnergy(i);
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  CHECK(valid());
  const auto kB = mConstants.kB();
  const auto mp = mConstants.protonMass();
  double Cv = kB/(mGamma1*mMolecularWeight*mp);
  specificHeat = Cv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (size_t i = 0u; i != soundSpeed.size(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
              const Field<Dimension, Scalar>& /*massDensity*/,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  CHECK(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho). 
// This is just the specific heat ratio times pressure for a gamma law gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    bulkModulus(i) = this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaLawGas<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    entropy(i) = pressure(massDensity(i), specificThermalEnergy(i))*safeInvVar(pow(massDensity(i), mGamma));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return this->applyPressureLimits(mGamma1*massDensity*specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const auto kB = mConstants.kB();
  const auto mp = mConstants.protonMass();
  return mGamma1*mMolecularWeight*mp/kB*specificThermalEnergy;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar temperature) const {
  CHECK(valid());
  const auto kB = mConstants.kB();
  const auto mp = mConstants.protonMass();
  return kB/(mGamma1*mMolecularWeight*mp)*temperature;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  CHECK(valid());
  const auto kB = mConstants.kB();
  const auto mp = mConstants.protonMass();
  return kB/(mGamma1*mMolecularWeight*mp);
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
soundSpeed(const Scalar /*massDensity*/,
           const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return sqrt(max(0.0, mGamma*mGamma1*specificThermalEnergy));
}

//------------------------------------------------------------------------------
// Calculate an isentropic bulk modulus.  
// This is just the specific heat ratio times pressure for a gamma law gas.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return mGamma*pressure(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return this->pressure(massDensity, specificThermalEnergy)*safeInvVar(pow(massDensity, mGamma));
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::gamma(const Scalar /*massDensity*/,
                              const Scalar /*specificThermalEnergy*/) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Get and set gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::gamma() const {
  return mGamma;
}

template<typename Dimension>
void
GammaLawGas<Dimension>::gamma(typename Dimension::Scalar gam) {
  mGamma = gam;
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Get and set the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
GammaLawGas<Dimension>::molecularWeight() const {
  return mMolecularWeight;
}

template<typename Dimension>
void
GammaLawGas<Dimension>::molecularWeight(typename Dimension::Scalar mu) {
  mMolecularWeight = mu;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GammaLawGas<Dimension>::valid() const {
  return (mGamma > 0.0 &&
          mMolecularWeight > 0.0 &&
          mGamma1 == mGamma - 1.0);
}

}
