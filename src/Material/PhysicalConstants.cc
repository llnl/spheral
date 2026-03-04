//---------------------------------Spheral++----------------------------------//
// PhysicalConstants -- Choose the physical units for a given Spheral run.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#include "PhysicalConstants.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
PhysicalConstants::
PhysicalConstants(const double unitLm,
                  const double unitMkg,
                  const double unitTsec,
                  const double unitTeK,
                  const double unitCcou): // Coulomb > Ampere
  mUnitLm(unitLm),
  mUnitMkg(unitMkg),
  mUnitTsec(unitTsec),
  mUnitTeK(unitTeK),
  mUnitCcou(unitCcou),
  UnitEnergyJ(unitMkg * FastMath::square(unitLm / unitTsec)),
  ProtonMass(mpMKS/unitMkg),
  ElectronMass(meMKS/unitMkg),
  ElectronCharge(qeMKS / unitCcou),
  GGravity(GMKS/(unitLm/unitMkg*FastMath::square(unitLm/unitTsec))),
  cLight(cMKS/(unitLm/unitTsec)),
  kBoltzmann(kBMKS*unitTeK/(unitMkg*FastMath::square(unitLm/unitTsec))),
  MolarGasConstant(RgasMKS*unitTeK/(unitMkg*FastMath::square(unitLm/unitTsec))),
  KelvinsToEnergyPerMole(unitMkg*FastMath::square(unitLm/unitTsec)/(kBMKS*unitTeK)*NAvogadro),
  UnitMassDensity(unitMkg/(unitLm*unitLm*unitLm)),
  Sigma(StefanBoltzmannMKS*unitTeK*unitTeK*unitTeK*unitTeK/unitMkg*unitTsec*unitTsec*unitTsec),
  BlackBody(4*StefanBoltzmannMKS*unitTeK*unitTeK*unitTeK*unitTeK/cMKS*unitTsec*unitTsec*unitLm/unitMkg),
  Planck(PlanckMKS*unitTsec/(unitMkg*unitLm*unitLm)),
  Mu0(Mu0MKS*unitCcou*unitCcou/(unitMkg*unitLm)) {
}

//------------------------------------------------------------------------------
// Define the base quantities -- we'll use the MKS SI standard.
// Source -- https://physics.nist.gov/cuu/Constants/index.html
// Old values are commented to the right of new values
//------------------------------------------------------------------------------
const double PhysicalConstants::mpMKS =     1.67262192369e-27;//1.672621777  // kg
const double PhysicalConstants::meMKS =     9.1093837015e-31;// 9.10938291  // kg
const double PhysicalConstants::qeMKS =     1.602176634e-19;//1.602176565 // Coulombs
const double PhysicalConstants::GMKS =      6.67430e-11; //6.67384     // N*m^2/kg^2
const double PhysicalConstants::cMKS =      2.99792458e8;    // m/s
const double PhysicalConstants::kBMKS =     1.380649e-23;//1.3806488   // J/K
const double PhysicalConstants::RgasMKS =   8.314462618;//8.3144621       // J/mole/K
const double PhysicalConstants::NAvogadro = 6.02214076e23;//6.02214129   // mol^-1
const double PhysicalConstants::StefanBoltzmannMKS = 5.670374419e-8;//5.67   // W/m^2/K^4
const double PhysicalConstants::PlanckMKS = 6.62607015e-34; // J*s
const double PhysicalConstants::Mu0MKS = 1.25663706127e-6; // N*A^-2

}
