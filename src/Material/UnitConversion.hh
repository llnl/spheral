//---------------------------------Spheral++----------------------------------//
// UnitConversion
//
// Convert between two unit systems
//----------------------------------------------------------------------------//

#ifndef __Spheral_UnitConversion_hh__
#define __Spheral_UnitConversion_hh__

#include "PhysicalConstants.hh"

namespace Spheral {

class UnitConversion {
public:

  // Constructor
  UnitConversion(const PhysicalConstants& unitsFrom,
                 const PhysicalConstants& unitsTo);

  // Fundamental quantities
  double length() const { return mLength; }
  double mass() const { return mMass; }
  double time() const { return mTime; }
  double temperature() const { return mTemperature; }
  double charge() const { return mCharge; }

  // Derived quantities
  double massDensity() const { return mMassDensity; }
  double specificEnergy() const { return mSpecificEnergy; }
  double energyDensity() const { return mEnergyDensity; }
  double pressure() const { return mEnergyDensity; }
  double bulkModulus() const { return mEnergyDensity; }
  double specificHeat() const { return mSpecificHeat; }
  double entropy() const { return mSpecificHeat; }
  double speed() const { return mSpeed; }
  double volume() const { return mVolume; }

private:
  // Input
  PhysicalConstants mUnitsFrom;
  PhysicalConstants mUnitsTo;

  // Precomputed unit data
  double mLength, mMass, mTime, mTemperature, mCharge;
  double mMassDensity, mSpecificEnergy, mEnergyDensity, mSpecificHeat, mSpeed, mVolume;
};

}

#endif

