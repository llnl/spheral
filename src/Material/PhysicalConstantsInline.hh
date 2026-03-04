namespace Spheral {

//------------------------------------------------------------------------------
// Unit length in meters.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitLengthMeters() const {
  return mUnitLm;
}

//------------------------------------------------------------------------------
// Unit mass in kg.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitMassKg() const {
  return mUnitMkg;
}

//------------------------------------------------------------------------------
// Unit time in sec.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitTimeSec() const {
  return mUnitTsec;
}

//------------------------------------------------------------------------------
// Unit temperature in kelvin.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitTemperatureKelvin() const {
  return mUnitTeK;
}

//------------------------------------------------------------------------------
// Unit charge in coulomb.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitChargeCoulomb() const {
  return mUnitCcou;
}

//------------------------------------------------------------------------------
// Unit energy in J.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitEnergyJ() const {
  return UnitEnergyJ;
}

//------------------------------------------------------------------------------
// Unit mass density in kg/m^3
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitMassDensity() const {
    return UnitMassDensity;
}

//------------------------------------------------------------------------------
// Proton mass.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::protonMass() const {
  return ProtonMass;
}

//------------------------------------------------------------------------------
// Electron mass.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::electronMass() const {
  return ElectronMass;
}

//------------------------------------------------------------------------------
// Electron charge.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::electronCharge() const {
  return ElectronCharge;
}

//------------------------------------------------------------------------------
// Gravitational constant.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::G() const {
  return GGravity;
}

//------------------------------------------------------------------------------
// Speed of light.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::c() const {
  return cLight;
}

//------------------------------------------------------------------------------
// Boltzmann constant.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::kB() const {
  return kBoltzmann;
}

//------------------------------------------------------------------------------
// Avagadro's constant.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::Navogadro() const {
  return NAvogadro;
}

//------------------------------------------------------------------------------
// molar constant
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::molarGasConstant() const {
  return MolarGasConstant;
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::kelvinsToEnergyPerMole() const {
  return KelvinsToEnergyPerMole;
}
    
//------------------------------------------------------------------------------
// Stefan-Boltzmann's constant
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::stefanBoltzmannConstant() const {
    return Sigma;
}

//------------------------------------------------------------------------------
// Black body constant
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::blackBodyConstant() const {
    return BlackBody;
}

//------------------------------------------------------------------------------
// Planck constant
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::planckConstant() const {
    return Planck;
}

//------------------------------------------------------------------------------
// Vacuum permeability (magnetic constant)
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::vacuumPermeability() const {
    return Mu0;
}

}
