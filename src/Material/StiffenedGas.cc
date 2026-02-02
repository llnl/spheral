//---------------------------------Spheral++----------------------------------//
// StiffenedGas -- The gamma law gas equation of state.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#include "StiffenedGas.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"

#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given gamma and mu.
//------------------------------------------------------------------------------
template<typename Dimension>
StiffenedGas<Dimension>::
StiffenedGas(const double gamma,
             const double P0,
             const double Cv,
             const PhysicalConstants& constants,
             const double minimumPressure,
             const double maximumPressure,
             const MaterialPressureMinType minPressureType,
             const double externalPressure):
  EquationOfState<Dimension>(constants, minimumPressure, maximumPressure, minPressureType, externalPressure),
  mGamma(gamma),
  mP0(P0),
  mCv(Cv) {
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StiffenedGas<Dimension>::
~StiffenedGas() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
              const Field<Dimension, Scalar>& /*massDensity*/,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  CHECK(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (isentropic) rho*c^2
//------------------------------------------------------------------------------
template<typename Dimension>
void
StiffenedGas<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    bulkModulus(i) = this->bulkModulus(massDensity(i),specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StiffenedGas<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    entropy(i) = this->entropy(massDensity(i),specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return this->applyPressureLimits(mGamma1*massDensity*specificThermalEnergy-mGamma*mP0);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return 1.0/mCv*specificThermalEnergy;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar temperature) const {
  CHECK(valid());
  return mCv*temperature;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::
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
StiffenedGas<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return sqrt(max(0.0, mGamma*mGamma1*(specificThermalEnergy-mP0/massDensity)));
}

//------------------------------------------------------------------------------
// isentropic bulk modulus (rho*c^2)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return max(0.0, mGamma*mGamma1*(massDensity*specificThermalEnergy-mP0));
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return (mP0+this->pressure(massDensity, specificThermalEnergy))*safeInvVar(pow(massDensity, mGamma));
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
StiffenedGas<Dimension>::gamma(const Scalar /*massDensity*/,
                              const Scalar /*specificThermalEnergy*/) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StiffenedGas<Dimension>::valid() const {
  return (mGamma > 0.0 &&
          mGamma1 == mGamma - 1.0);
}

}
