//---------------------------------Spheral++----------------------------------//
// ArtificialViscosityView -- The base class for the ArtificialViscosity view class
// Spheral++.
//
// Created by LDO, Wed Oct 15 14:16:43 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosityView__
#define __Spheral_ArtificialViscosityView__

#include <utility>
#include <typeindex>

namespace Spheral {

template<typename Dimension>
class ArtificialViscosityView {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor
  SPHERAL_HOST_DEVICE
  ArtificialViscosityView() = default;
  SPHERAL_HOST_DEVICE
  ArtificialViscosityView(const Scalar Clinear,
                          const Scalar Cquadratic,
                          const bool   BalsaraShearCorrection = false,
                          const Scalar Epsilon2 = 1.0E-2,
                          const Scalar NegligibleSoundSpeed = 1.0E-10) :
    mClinear(Clinear),
    mCquadratic(Cquadratic),
    mBalsaraShearCorrection(BalsaraShearCorrection),
    mEpsilon2(Epsilon2),
    mNegligibleSoundSpeed(NegligibleSoundSpeed) {}

  // Apparently virtual destructors keep a vptr/vtable in the
  // device view even after QPiij is not virtual
  SPHERAL_HOST_DEVICE
  ~ArtificialViscosityView() = default;

  //...........................................................................
  // Methods
  // Calculate the curl of the velocity given the stress tensor.
  SPHERAL_HOST_DEVICE
  Scalar curlVelocityMagnitude(const Tensor& DvDx) const;

  // Find the Balsara shear correction multiplier
  SPHERAL_HOST_DEVICE
  Scalar calcBalsaraShearCorrection(const Tensor& DvDx,
                                    const SymTensor& H,
                                    const Scalar& cs) const;

  SPHERAL_HOST_DEVICE Scalar Cl()                     const { return mClinear; }
  SPHERAL_HOST_DEVICE Scalar Cq()                     const { return mCquadratic; }
  SPHERAL_HOST_DEVICE Scalar epsilon2()               const { return mEpsilon2; }
  SPHERAL_HOST_DEVICE bool   balsaraShearCorrection() const { return mBalsaraShearCorrection; }
  SPHERAL_HOST_DEVICE Scalar negligibleSoundSpeed()   const { return mNegligibleSoundSpeed; }
  SPHERAL_HOST_DEVICE void Cl(Scalar x)                     { mClinear = x; }
  SPHERAL_HOST_DEVICE void Cq(Scalar x)                     { mCquadratic = x; }
  SPHERAL_HOST_DEVICE void epsilon2(Scalar x)               { mEpsilon2 = x; }
  SPHERAL_HOST_DEVICE void balsaraShearCorrection(bool x)   { mBalsaraShearCorrection = x; }
  SPHERAL_HOST_DEVICE void negligibleSoundSpeed(Scalar x)   { REQUIRE(x > 0.0); mNegligibleSoundSpeed = x; }

protected:
  Scalar mClinear = 0.0;
  Scalar mCquadratic = 0.0;
  // Switch for the Balsara shear correction.
  bool mBalsaraShearCorrection = false;
  // Parameters for the Q limiter.
  Scalar mEpsilon2 = 1.0E-2;
  Scalar mNegligibleSoundSpeed = 1.0E-10;
};

}

#include "ArtificialViscosityViewInline.hh"

#endif
