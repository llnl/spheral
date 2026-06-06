//---------------------------------Spheral++----------------------------------//
// iSALEROCKYieldStrength -- Implements a pressure dependent yield strength
// model appropriate for geological materials.
//
// This is based on the ROCK model in iSALE.  Described in a few places, though
// it seems to keep changing slightly.
// 
//    See Collins, Melosh, Ivanov, 2004 Appendix, MAPS
//    Raducan et al. 2020 (projectile shape effects paper)
//
// Created by JMO, Mon Jul 12 15:05:00 PDT 2021
//----------------------------------------------------------------------------//
#include "iSALEROCKStrength.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"
#include "SolidEquationOfState.hh"
#include "Field/Field.hh"
#include "NodeList/SolidNodeList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
iSALEROCKStrength<Dimension>::
iSALEROCKStrength(const StrengthModel<Dimension>& shearModulusModel,
                  const double Yi0,                                   // Intact strength at zero pressure     
                  const double Yiinf,                                 // Intact strength at infinite pressure 
                  const double fi,                                    // Intact internal friction coefficient 
                  const double Yd0,                                   // Damaged strength at zero pressure    
                  const double Ydinf,                                 // Damaged strength at infinite pressure
                  const double fd):                                   // Damaged internal friction coefficient
  StrengthModel<Dimension>(),
  mShearModulusModel(shearModulusModel),
  mYi0(Yi0),
  mYiinf(Yiinf),
  mfi(fi),
  mYd0(Yd0),
  mYdinf(Ydinf),
  mfd(fd) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
iSALEROCKStrength<Dimension>::
~iSALEROCKStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
iSALEROCKStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& specificThermalEnergy,
             const Field<Dimension, Scalar>& pressure,
             const Field<Dimension, SymTensor>& damage) const {
  mShearModulusModel.shearModulus(shearModulus, density, specificThermalEnergy, pressure, damage);
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
iSALEROCKStrength<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
              const Field<Dimension, Scalar>& pressure,
              const Field<Dimension, Scalar>& /*plasticStrain*/,
              const Field<Dimension, Scalar>& /*plasticStrainRate*/,
              const Field<Dimension, SymTensor>& damage) const {

  const auto n = density.numInternalElements();
  const auto YdiffInv = safeInvVar(mYiinf - mYi0);
#pragma omp for
  for (unsigned i = 0; i < n; ++i) {
    const auto Pi = std::max(0.0, pressure(i));
    const auto Yi = mYi0 + mfi*Pi/(1.0 + mfi*Pi*YdiffInv);
    const auto Yd = std::min(mYdinf, mYd0 + mfd*Pi);
    CHECK(Yi >= 0.0 and Yd >= 0.0);
    const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
    yieldStrength(i) = (1.0 - Di)*Yi + Di*Yd;
    CHECK(yieldStrength(i) >= 0.0);
  }
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
iSALEROCKStrength<Dimension>::
soundSpeed(Field<Dimension, Scalar>& soundSpeed,
           const Field<Dimension, Scalar>& density,
           const Field<Dimension, Scalar>& specificThermalEnergy,
           const Field<Dimension, Scalar>& pressure,
           const Field<Dimension, Scalar>& fluidSoundSpeed,
           const Field<Dimension, SymTensor>& damage) const {
  mShearModulusModel.soundSpeed(soundSpeed, density, specificThermalEnergy, pressure, fluidSoundSpeed, damage);
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
const StrengthModel<Dimension>&
iSALEROCKStrength<Dimension>::
shearModulusModel() const {
  return mShearModulusModel;
}

template<typename Dimension>
double
iSALEROCKStrength<Dimension>::
Yi0() const {
  return mYi0;
}

template<typename Dimension>
double
iSALEROCKStrength<Dimension>::
Yiinf() const {
  return mYiinf;
}

template<typename Dimension>
double
iSALEROCKStrength<Dimension>::
fi() const {
  return mfi;
}

template<typename Dimension>
double
iSALEROCKStrength<Dimension>::
Yd0() const {
  return mYd0;
}

template<typename Dimension>
double
iSALEROCKStrength<Dimension>::
Ydinf() const {
  return mYdinf;
}

template<typename Dimension>
double
iSALEROCKStrength<Dimension>::
fd() const {
  return mfd;
}

}
