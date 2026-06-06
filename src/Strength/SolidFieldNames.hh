//---------------------------------Spheral++----------------------------------//
// SolidFieldNames -- A collection of standard Field names for solid materials.
//
// Created by JMO, Wed Sep 8 11:05:49 2004
//----------------------------------------------------------------------------//
#ifndef _Spheral_SolidFieldNames_
#define _Spheral_SolidFieldNames_

#include <string>

namespace Spheral {

struct SolidFieldNames {
  const inline static std::string deviatoricStress = "deviatoric stress";
  const inline static std::string deviatoricStressTT = "deviatoric stress theta theta";
  const inline static std::string plasticStrain = "plastic strain";
  const inline static std::string plasticStrainRate = "plastic strain rate";
  const inline static std::string scalarDamage = "scalar damage";
  const inline static std::string tensorDamage = "tensor damage";
  const inline static std::string damageCoupling = "damage coupling";
  const inline static std::string strain = "strain";
  const inline static std::string strainTensor = "tensor strain";
  const inline static std::string effectiveStrainTensor = "effective tensor strain";
  const inline static std::string bulkModulus = "bulk modulus";
  const inline static std::string shearModulus = "shear modulus";
  const inline static std::string YoungsModulus = "Youngs modulus";
  const inline static std::string longitudinalSoundSpeed = "longitudinal sound speed";
  const inline static std::string yieldStrength = "yield strength";
  const inline static std::string flaws = "flaws";
  const inline static std::string numFlaws = "num flaws";
  const inline static std::string minFlaw = "minimum flaw";
  const inline static std::string maxFlaw = "maximum flaw";
  const inline static std::string initialVolume = "initial volume";
  const inline static std::string randomGenerator = "random generator";
  const inline static std::string porositySolidDensity = "porosity solid mass density";
  const inline static std::string porosityAlpha = "porosity alpha";
  const inline static std::string porosityStrain = "porosity strain";
  const inline static std::string porosityAlpha0 = "initial porosity alpha";
  const inline static std::string porosityc0 = "initial porosity sound speed";
  const inline static std::string fDSjutzi = "f deviatoric stress factor";
  const inline static std::string fragmentIDs = "fragment index";
  const inline static std::string particleTypes = "particle type";
  const inline static std::string meltSpecificEnergy = "melt specific energy";
  const inline static std::string mask = "mask";
  const inline static std::string damagedPressure = "damaged pressure";
};

}

#endif
