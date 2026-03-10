#include "Geodyn.hh"
#include "Field/Field.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>
using std::vector;

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Geodyn<Dimension>::
Geodyn(const PhysicalConstants& constants,
       const double minimumPressure,
       const double maximumPressure,
       const MaterialPressureMinType minPressureType):
  Physics<Dimension>(),
  SolidEquationOfState<Dimension>(1.0, 1.0, 1.0,
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minPressureType),
  StrengthModel<Dimension>(),
  mGeodynUnits(0.01,    // cm expressed as meters.
               0.001,   // g expressed in kg.
               1.0),    // sec in secs.
  mRhoConv(1.0),
  mTconv(1.0),
  mPconv(1.0),
  mEconv(1.0),
  mCVconv(1.0),
  mVelConv(1.0) {
  
  // Build our unit conversion factors.
  const double lconv = mGeodynUnits.unitLengthMeters() / constants.unitLengthMeters(),
               mconv = mGeodynUnits.unitMassKg()       / constants.unitMassKg(),
               tconv = mGeodynUnits.unitTimeSec()      / constants.unitTimeSec();
  mRhoConv = mconv/(lconv*lconv*lconv);
  mPconv = mconv/(lconv*tconv*tconv);
  mTconv = 1.0;   // Kelvin?
  mEconv = FastMath::square(lconv/tconv);
  mCVconv = mEconv/mTconv;
  mVelConv = lconv/tconv;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Geodyn<Dimension>::
~Geodyn() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
setEntropy(Field<Dimension, Scalar>& bulkModulus,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
}

//------------------------------------------------------------------------------
// valid
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Geodyn<Dimension>::
valid() const {
  return true;
}

//------------------------------------------------------------------------------
// shear modulus
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& specificThermalEnergy,
             const Field<Dimension, Scalar>& pressure) const {
}

//------------------------------------------------------------------------------
// yield strength
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& specificThermalEnergy,
              const Field<Dimension, Scalar>& pressure,
              const Field<Dimension, Scalar>& plasticStrain,
              const Field<Dimension, Scalar>& plasticStrainRate) const {
}

//------------------------------------------------------------------------------
// sound speed with strength
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
soundSpeed(Field<Dimension, Scalar>& soundSpeed,
           const Field<Dimension, Scalar>& density,
           const Field<Dimension, Scalar>& specificThermalEnergy,
           const Field<Dimension, Scalar>& pressure,
           const Field<Dimension, Scalar>& fluidSoundSpeed) const {
}

//------------------------------------------------------------------------------
// evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
}

//------------------------------------------------------------------------------
// dt
//------------------------------------------------------------------------------
template<typename Dimension>
typename Geodyn<Dimension>::TimeStepType
Geodyn<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return TimeStepType(1.0e100, "Geodyn no timestep vote.");
}

//------------------------------------------------------------------------------
// registerState
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
}

//------------------------------------------------------------------------------
// registerDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
Geodyn<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

}

}
