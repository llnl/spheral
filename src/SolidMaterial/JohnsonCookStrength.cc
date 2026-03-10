//---------------------------------Spheral++----------------------------------//
// JohnsonCookStrength -- Implements the Johnson-Cook strength model.
//   ** Need reference **
//
// Created by JMO, Tue Jul 28 13:39:50 PDT 2015
//----------------------------------------------------------------------------//
#include "JohnsonCookStrength.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"
#include "SolidEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookStrength<Dimension>::
JohnsonCookStrength(const SolidEquationOfState<Dimension>& eos,
                    const StrengthModel<Dimension>& shearModulusModel,
                    const double A,
                    const double B,
                    const double C,
                    const double C4,
                    const double m,
                    const double nhard,
                    const double epsdot0,
                    const double epsdotmin,
                    const double Tmelt,
                    const double Troom,
                    const double mu0,
                    const bool shearModulusScaling):
  mEOSPtr(&eos),
  mShearModulusModelPtr(&shearModulusModel),
  mA(A),
  mB(B),
  mC(C),
  mC4(C4),
  mm(m),
  mnhard(nhard),
  mEpsdot0(epsdot0),
  mEpsdotmin(epsdotmin),
  mTmelt(Tmelt),
  mTroom(Troom),
  mmu0(mu0),
  mShearModulusScaling(shearModulusScaling) {
  VERIFY2(mEpsdot0 > 0.0,
          "JohnsonCookStrength ERROR: reference strain-rate must be greater than zero.");
  VERIFY2(mTmelt > mTroom,
          "JohnsonCookStrength ERROR: Tmelt must be greater than or equal Troom.");
  VERIFY2((not shearModulusScaling) or mu0 > 0.0,
          "JohnsonCookStrength ERROR: require mu0 >= 0.0 if using shear modulus scaling.");
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookStrength<Dimension>::
~JohnsonCookStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& specificThermalEnergy,
             const Field<Dimension, Scalar>& pressure,
             const Field<Dimension, SymTensor>& damage) const {
  mShearModulusModelPtr->shearModulus(shearModulus, density, specificThermalEnergy, pressure, damage);
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookStrength<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& specificThermalEnergy,
              const Field<Dimension, Scalar>& pressure,
              const Field<Dimension, Scalar>& plasticStrain,
              const Field<Dimension, Scalar>& plasticStrainRate,
              const Field<Dimension, SymTensor>& damage) const {
  Field<Dimension, Scalar> T("temperature", yieldStrength.nodeList());
  mEOSPtr->setTemperature(T, density, specificThermalEnergy);
  const auto n = yieldStrength.numInternalElements();
#pragma omp for
  for (auto i = 0u; i < n; ++i) {
    const auto Tstar = std::max(0.0, std::min(1.0, T(i) - mTroom)/(mTmelt - mTroom));
    const auto fmelt = std::max(0.0, std::min(1.0, 1.0 - pow(Tstar, mm)));
    yieldStrength(i) = 
      ((mA + mB*pow(plasticStrain(i), mnhard))*
       (1.0 + mC*log(max(mEpsdotmin, plasticStrainRate(i))/mEpsdot0))*fmelt +
       mC4*pressure(i));
    const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
    yieldStrength(i) = (1.0 - Di)*yieldStrength(i);
  }

  // Optionally scale by the relative shear modulus.
  if (mShearModulusScaling) {
    Field<Dimension, Scalar> mu("shear modulus", yieldStrength.nodeList());
    mShearModulusModelPtr->shearModulus(mu, density, specificThermalEnergy, pressure, damage);
    for (auto i = 0u; i != yieldStrength.numInternalElements(); ++i) {
      yieldStrength(i) *= mu(i)*safeInvVar(mmu0);
    }
  }
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookStrength<Dimension>::
soundSpeed(Field<Dimension, Scalar>& soundSpeed,
           const Field<Dimension, Scalar>& density,
           const Field<Dimension, Scalar>& specificThermalEnergy,
           const Field<Dimension, Scalar>& pressure,
           const Field<Dimension, Scalar>& fluidSoundSpeed,
           const Field<Dimension, SymTensor>& damage) const {
  mShearModulusModelPtr->soundSpeed(soundSpeed, density, specificThermalEnergy, pressure, fluidSoundSpeed, damage);
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
C() const {
  return mC;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
C4() const {
  return mC4;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
m() const {
  return mm;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
nhard() const {
  return mnhard;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
epsdot0() const {
  return mEpsdot0;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
epsdotmin() const {
  return mEpsdotmin;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
Tmelt() const {
  return mTmelt;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
Troom() const {
  return mTroom;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
mu0() const {
  return mmu0;
}

template<typename Dimension>
bool
JohnsonCookStrength<Dimension>::
shearModulusScaling() const {
  return mShearModulusScaling;
}

}
