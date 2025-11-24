//---------------------------------Spheral++----------------------------------//
// UnitConversion
//
// Convert between two unit systems
//----------------------------------------------------------------------------//

#include "UnitConversion.hh"

namespace Spheral {

UnitConversion::
UnitConversion(const PhysicalConstants& unitsFrom,
               const PhysicalConstants& unitsTo):
  mUnitsFrom(unitsFrom),
  mUnitsTo(unitsTo) {
  // From these units (comments to right of code for common-sense check)
  const auto l0 = unitsFrom.unitLengthMeters(); // 0.01 (cm)
  const auto m0 = unitsFrom.unitMassKg();       // 0.001 (g)
  const auto t0 = unitsFrom.unitTimeSec();
  const auto k0 = unitsFrom.unitTemperatureKelvin();
  const auto c0 = unitsFrom.unitChargeCoulomb();

  // To these units
  const auto l = unitsTo.unitLengthMeters(); // 1000 (km)
  const auto m = unitsTo.unitMassKg();       // 1.0 (kg)
  const auto t = unitsTo.unitTimeSec();
  const auto k = unitsTo.unitTemperatureKelvin();
  const auto c = unitsTo.unitChargeCoulomb();

  // Conversion factors
  const double lC = l0 / l; // 0.01/1000=0.00001 cm/km
  const double mC = m0 / m; // 0.001/1.0=0.001 g/kg
  const double tC = t0 / t;
  const double kC = k0 / k;
  const double cC = c0 / c;

  // Store these
  mLength = lC; // From cm to km: multiply by 0.00001
  mMass = mC;   // From g to kg:  multiply by 0.001
  mTime = tC;
  mTemperature = kC;
  mCharge = cC;

  // Computed factors
  mMassDensity = mC / (lC * lC * lC); // 1e-3/(1e-5)^3=10^12 g/cm^3 / (kg/km^3)
  mSpecificEnergy = (lC * lC) / (tC * tC);
  mEnergyDensity = mC / (lC * tC * tC);
  mSpecificHeat = (lC * lC) / (tC * tC * kC);
  mSpeed = lC / tC;
  mVolume = 1.0 / (lC * lC * lC);
}

}
