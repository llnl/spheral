//---------------------------------Spheral++----------------------------------//
// ANEOS -- An interface to the ANEOS equation of state from 
// Melosh amony others.  This is a C++ wrapper around calls to the underlying
// ANEOS fortran library.  The user must provide the ancillary file containing
// parameters for ANEOS in its expected format.  Also make sure to keep units
// consistent in this separate input file with what you give the C++ EOS here!
//
// Created by JMO, Tue Apr 23 14:55:28 PDT 2013
//----------------------------------------------------------------------------//
#include "ANEOS.hh"
#include "Field/Field.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectRoot.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralMessage.hh"

#include <iostream>
#include <ctime>
using std::vector;
using std::string;
using std::pair;


// Fortran baby!
extern "C" {
  void aneos_initialize_(char* filename, int* num, int* izetl);
  void call_aneos_(int* matnum, double* T, double* rho, 
                   double* P, double* E, double* S, double* CV, double* DPDT, double* DPDR, double* cs);
  double get_aneos_atomicweight_(int* matnum);
  double get_aneos_referencedensity_(int* matnum);
}

namespace Spheral {

// //------------------------------------------------------------------------------
// // Define an inline common method to handle calling ANEOS1.
// //------------------------------------------------------------------------------
// inline
// void call_ANEOS1(double& T, double& rho, 
//                  double& P, double& E, double& S, double& CV, double& DPDT, double& DPDR,
//                  int& mat) {
//   anesqt_.sqts[0] = std::sqrt(T);  // What the?  Some stuff is passed in in the common block
//   aneos1_(&T, &rho, 
//           &P, &E, &S, &CV, &DPDT, &DPDR,
//           &mat);
// }

namespace { // anonymous

using InterpolatorType = CubicHermiteInterpolator;
using BiInterpolatorType = BiCubicInterpolator;

//------------------------------------------------------------------------------
// A functor to compute eps(rho, T) for use building the interpolation table
//------------------------------------------------------------------------------
class epsFunc {
public:
  epsFunc(int matNum, double rhoConv, double Tconv, double Econv):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mEconv(Econv) {};
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = y/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return eps * mEconv;
  }

public:
  mutable int mMatNum;
  double mRhoConv, mTconv, mEconv;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// A functor to compute T(rho, eps) for use building the interpolation table
//------------------------------------------------------------------------------
class Tfunc {
public:
  Tfunc(double Tmin, 
        double Tmax,
        const InterpolatorType& epsMinInterp,
        const InterpolatorType& epsMaxInterp,
        const BiInterpolatorType& epsInterp,
        const epsFunc& Feps):
    mTmin(Tmin),
    mTmax(Tmax),
    mEpsMinInterp(epsMinInterp),
    mEpsMaxInterp(epsMaxInterp),
    mEpsInterp(epsInterp) {}

  double operator()(const double rho, const double eps) const {
    if (eps < mEpsMinInterp(rho)) {
      // cerr << " **> BAIL low: " << eps << " " << mEpsMinInterp(rho) << " -> " << mTmin << endl;
      return mTmin;
    } else if (eps > mEpsMaxInterp(rho)) {
      // cerr << " **> BAIL high: " << eps << " " << mEpsMaxInterp(rho) << " -> " << mTmax << endl;
      return mTmax;
    } else {
      const auto FT = Trho_func(rho, eps, mEpsInterp);
      const double FTmin = FT(mTmin), FTmax = FT(mTmax);
      const double result = (FTmin*FTmax > 0.0 ?
                             (abs(FTmin) < abs(FTmax) ? mTmin : mTmax) :
                             bisectRoot(FT, mTmin, mTmax, 1.0e-10, 1.0e-10, 200u));
      // cerr << " **> (" << rho << " " << eps << ") [" << mEpsMinInterp(rho) << " " << mEpsMaxInterp(rho) << "] " << FT(mTmin) << " " << FT(mTmax) << " -> " << result << endl;
      return result;
    }
  }

private:
  double mTmin, mTmax;
  const InterpolatorType& mEpsMinInterp, mEpsMaxInterp;
  const BiInterpolatorType& mEpsInterp;

  // We need to make a single argument functor for eps(T) given a fixed rho
  class Trho_func {
    double mrho, meps;
    const BiInterpolatorType& mEpsInterp;
  public:
    Trho_func(const double rho,
              const double eps,
              const BiInterpolatorType& epsInterp):
      mrho(rho),
      meps(eps),
      mEpsInterp(epsInterp) {}
    double operator()(const double T) const { return mEpsInterp(mrho, T) - meps; }
  };
};

//------------------------------------------------------------------------------
// A wrapper around the T(rho, eps) interpolator with safe limiting
//------------------------------------------------------------------------------
class Textrapolator {
public:
  Textrapolator(const InterpolatorType& epsMinInterp,
                const InterpolatorType& epsMaxInterp,
                const BiInterpolatorType& Tinterp):
    mEpsMinInterp(epsMinInterp),
    mEpsMaxInterp(epsMaxInterp),
    mTinterp(Tinterp) {}

  double operator()(const double rho, const double eps) const {
    const auto eps0 = mEpsMinInterp(rho);
    const auto eps1 = mEpsMaxInterp(rho);
    if (eps < eps0) {
      return mTinterp(rho, eps0);
    } else if (eps > eps1) {
      return mTinterp(rho, eps1);
    } else {
      return mTinterp(rho, eps);
    }
  }

  double dTdeps(const double rho, const double eps) const {
    const auto eps0 = mEpsMinInterp(rho);
    const auto eps1 = mEpsMaxInterp(rho);
    if (eps < eps0) {
      return mTinterp.prime_y(rho, eps0);
    } else if (eps > eps1) {
      return mTinterp.prime_y(rho, eps1);
    } else {
      return mTinterp.prime_y(rho, eps);
    }
  }

private:
  const InterpolatorType& mEpsMinInterp, mEpsMaxInterp;
  const BiInterpolatorType& mTinterp;
};

//------------------------------------------------------------------------------
// P(rho, eps)
//------------------------------------------------------------------------------
class Pfunc {
public:
  Pfunc(int matNum, double rhoConv, double Tconv, double Pconv,
        const Textrapolator& Textra):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mTextra(Textra) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = mTextra(x, y)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (P != P) P = 0.0;   // Hack for ANEOS giving us a NaN
    return P * mPconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const Textrapolator& mTextra;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// dPdeps(rho, eps)
//------------------------------------------------------------------------------
class dPdeps_func {
public:
  dPdeps_func(int matNum, double rhoConv, double Tconv, double Pconv,
              const Textrapolator& Textra):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mTextra(Textra) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = mTextra(x, y)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (DPDT != DPDT) DPDT = 0.0;   // Hack for ANEOS giving us a NaN
    return DPDT * mPconv/mTconv * mTextra.dTdeps(x, y);  // Last term is (partial T)/(partial eps) in Spheral units
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const Textrapolator& mTextra;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// dPdrho(rho, eps)
//------------------------------------------------------------------------------
class dPdrho_func {
public:
  dPdrho_func(int matNum, double rhoConv, double Tconv, double Pconv,
              const Textrapolator& Textra):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mTextra(Textra) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = mTextra(x, y)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (DPDR != DPDR) DPDR = 0.0;   // Hack for ANEOS giving us a NaN
    return DPDR * mPconv/mRhoConv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const Textrapolator& mTextra;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// cV(rho, T)
//------------------------------------------------------------------------------
class cVfunc {
public:
  cVfunc(int matNum, double rhoConv, double Tconv, double cVconv):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mCVconv(cVconv) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = y/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (cV != cV) cV = 1.0e-30;   // Hack for ANEOS giving us a NaN
    return cV * mCVconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mCVconv;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// cs(rho, eps)
//------------------------------------------------------------------------------
class csfunc {
public:
  csfunc(int matNum, double rhoConv, double Tconv, double velConv,
         const Textrapolator& Textra):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mVelConv(velConv),
    mTextra(Textra) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = mTextra(x, y)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (cs != cs) cs = 1.0e-30;   // Hack for ANEOS giving us a NaN
    return cs * mVelConv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mVelConv;
  const Textrapolator& mTextra;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// DPDrho(rho, eps)
//------------------------------------------------------------------------------
class Kfunc {
public:
  Kfunc(int matNum, double rhoConv, double Tconv, double Pconv,
        const Textrapolator& Textra):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mTextra(Textra) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = mTextra(x, y)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (DPDR != DPDR) DPDR = 0.0;   // Hack for ANEOS giving us a NaN
    return std::abs(rho * DPDR * mPconv);
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const Textrapolator& mTextra;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// s(rho, eps) (entropy)
//------------------------------------------------------------------------------
class sfunc {
public:
  sfunc(int matNum, double rhoConv, double Tconv, double Sconv,
        const Textrapolator& Textra):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mSconv(Sconv),
    mTextra(Textra) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = mTextra(x, y)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    if (S != S) S = 0.0;   // Hack for ANEOS giving us a NaN
    return S * mSconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mSconv;
  const Textrapolator& mTextra;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

} // anonymous

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ANEOS<Dimension>::
ANEOS(const int materialNumber,
      const unsigned numRhoVals,
      const unsigned numTvals,
      const double rhoMin,
      const double rhoMax,
      const double Tmin,
      const double Tmax,
      const PhysicalConstants& constants,
      const double externalPressure,
      const double minimumPressure,
      const double maximumPressure,
      const double minimumPressureDamage,
      const MaterialPressureMinType minPressureType,
      const bool useInterpolation):
  SolidEquationOfState<Dimension>(get_aneos_referencedensity_(const_cast<int*>(&materialNumber)),  // not in the right units yet!
                                  0.0,                                           // dummy etamin
                                  std::numeric_limits<double>::max(),            // dummy etamax
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minimumPressureDamage,
                                  minPressureType,
                                  externalPressure),
  mUseInterpolation(useInterpolation),
  mMaterialNumber(materialNumber),
  mNumRhoVals(numRhoVals),
  mNumTvals(numTvals),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax),
  mTmin(Tmin),
  mTmax(Tmax),
  mEpsMin(std::numeric_limits<double>::max()),
  mEpsMax(std::numeric_limits<double>::min()),
  mEpsMinInterp(std::make_shared<InterpolatorType>()),
  mEpsMaxInterp(std::make_shared<InterpolatorType>()),
  mEpsInterp(std::make_shared<BiInterpolatorType>()),
  mTinterp(std::make_shared<BiInterpolatorType>()),
  mPinterp(std::make_shared<BiInterpolatorType>()),
  mCVinterp(std::make_shared<BiInterpolatorType>()),
  mCSinterp(std::make_shared<BiInterpolatorType>()),
  mKinterp(std::make_shared<BiInterpolatorType>()),
  mSinterp(std::make_shared<BiInterpolatorType>()),
  mDPDepsInterp(std::make_shared<BiInterpolatorType>()),
  mDPDRinterp(std::make_shared<BiInterpolatorType>()),
  mANEOSunits(0.01,   // cm expressed as meters.
              0.001,  // g expressed in kg.
              1.0),   // sec in secs.
  mRhoConv(1.0),
  mTconv(1.0),
  mPconv(1.0),
  mEconv(1.0),
  mCVconv(1.0),
  mVelConv(1.0),
  mSconv(1.0),
  mAtomicWeight(0.0) {
  VERIFY2(numRhoVals > 1,
          "ANEOS ERROR : specify numRhoVals > 1");
  VERIFY2(numTvals > 1,
          "ANEOS ERROR : specify numTvals > 1");
  VERIFY2(rhoMin < rhoMax,
          "ANEOS ERROR : specify rhoMin < rhoMax");
  VERIFY2(Tmin < Tmax,
          "ANEOS ERROR : specify Tmin < Tmax");
  VERIFY2(Tmin > 0.0,
          "ANEOS ERROR : specify Tmin > 0.0");

  // // Convert temperature range to log space.
  // mTmin = log(mTmin);
  // mTmax = log(mTmax);
  
  // Look up the atomic weight.
  mAtomicWeight = get_aneos_atomicweight_(&mMaterialNumber);
  VERIFY2(mAtomicWeight > 0.0, 
          "ANEOS ERROR : bad atomic weight for material " << mMaterialNumber << " : " << mAtomicWeight);

  // Build our unit conversion factors.  After looking through the ANEOS source some it appears to me 
  // that they use mostly CGS units, except for temperatures which are in eV.
  const double lconv = mANEOSunits.unitLengthMeters() / constants.unitLengthMeters(),
               mconv = mANEOSunits.unitMassKg() / constants.unitMassKg(),
               tconv = mANEOSunits.unitTimeSec() / constants.unitTimeSec();
  mRhoConv = mconv/(lconv*lconv*lconv);
  mPconv = mconv/(lconv*tconv*tconv);
  mTconv = 1.160452e4; // eV/kB (CRC 2013)
  mEconv = FastMath::square(lconv/tconv);
  mCVconv = mEconv/mTconv;
  mVelConv = lconv/tconv;
  mSconv = mEconv/(mconv*mTconv);

  // Fix the reference density.
  this->referenceDensity(this->referenceDensity() * mRhoConv);

  // Build the interpolator for the mininum specific energy along the min density line.
  // Also find the maximum specific energy.
  const auto Feps = epsFunc(mMaterialNumber, mRhoConv, mTconv, mEconv);
  const auto drho = (mRhoMax - mRhoMin)/(mNumRhoVals - 1u);
  const auto dT = (mTmax - mTmin)/(mNumTvals - 1u);
  CHECK(drho > 0.0);
  vector<double> epsMinVals(mNumRhoVals), epsMaxVals(mNumRhoVals);
  for (auto i = 0u; i < mNumRhoVals; ++i) {

    // ANEOS seems to often have ranges of constant return values for low temperatures.  We're
    // looking for where the EOS starts to actually have structure.
    auto j = 0u;
    auto eps0 = Feps(mRhoMin + i*drho, mTmin);
    while (j++ < mNumTvals and
           abs(eps0*safeInv(Feps(mRhoMin + i*drho, mTmin + j*dT)) - 1.0) < 1.0e-8) {
      eps0 = Feps(mRhoMin + i*drho, mTmin + j*dT);
    }
    epsMinVals[i] = eps0;
    epsMaxVals[i] = Feps(mRhoMin + i*drho, mTmax);
    mEpsMin = min(mEpsMin, epsMinVals[i]);
    mEpsMax = max(mEpsMax, epsMaxVals[i]);
  }
  mEpsMinInterp->initialize(mTmin, mTmax, epsMinVals);
  mEpsMaxInterp->initialize(mTmin, mTmax, epsMaxVals);
  mEpsMinInterp->makeMonotonic();
  mEpsMaxInterp->makeMonotonic();

  // Build the interpolation function for eps(rho, T)
  auto t0 = clock();
  mEpsInterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                     mTmin, mTmax,
                                                     mNumRhoVals, mNumTvals, Feps);
  SpheralMessage("ANEOS: Time to build epsInterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  // Now the hard inversion method for looking up T(rho, eps)
  t0 = clock();
  const auto Ftemp = Tfunc(mTmin, mTmax, *mEpsMinInterp, *mEpsMaxInterp, *mEpsInterp, Feps);
  mTinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                   mEpsMin, mEpsMax,
                                                   mNumRhoVals, mNumTvals, Ftemp);
  SpheralMessage("ANEOS: Time to build Tinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  // And finally the interpolators for most of our derived quantities
  t0 = clock();
  const auto Textra = Textrapolator(*mEpsMinInterp, *mEpsMaxInterp, *mTinterp);
  const auto Fpres = Pfunc(mMaterialNumber, mRhoConv, mTconv, mPconv, Textra);
  mPinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                   mEpsMin, mEpsMax,
                                                   mNumRhoVals, mNumTvals, Fpres);
  SpheralMessage("ANEOS: Time to build Pinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  const auto Fcv = cVfunc(mMaterialNumber, mRhoConv, mTconv, mCVconv);
  mCVinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                    mTmin, mTmax,
                                                    mNumRhoVals, mNumTvals, Fcv);
  SpheralMessage("ANEOS: Time to build CVinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  const auto Fcs = csfunc(mMaterialNumber, mRhoConv, mTconv, mVelConv, Textra);
  mCSinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                    mEpsMin, mEpsMax,
                                                    mNumRhoVals, mNumTvals, Fcs);
  SpheralMessage("ANEOS: Time to build CSinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  const auto FK = Kfunc(mMaterialNumber, mRhoConv, mTconv, mPconv, Textra);
  mKinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                   mEpsMin, mEpsMax,
                                                   mNumRhoVals, mNumTvals, FK);
  SpheralMessage("ANEOS: Time to build Kinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  const auto Fs = sfunc(mMaterialNumber, mRhoConv, mTconv, mSconv, Textra);
  mSinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                   mEpsMin, mEpsMax,
                                                   mNumRhoVals, mNumTvals, Fs);
  SpheralMessage("ANEOS: Time to build Sinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  const auto Fdpdeps = dPdeps_func(mMaterialNumber, mRhoConv, mTconv, mPconv, Textra);
  mDPDepsInterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                        mEpsMin, mEpsMax,
                                                        mNumRhoVals, mNumTvals, Fdpdeps);
  SpheralMessage("ANEOS: Time to build DPDUinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  const auto Fdpdrho = dPdrho_func(mMaterialNumber, mRhoConv, mTconv, mPconv, Textra);
  mDPDRinterp = std::make_shared<BiInterpolatorType>(mRhoMin, mRhoMax,
                                                      mEpsMin, mEpsMax,
                                                      mNumRhoVals, mNumTvals, Fdpdrho);
  SpheralMessage("ANEOS: Time to build DPDRinterp: " << double(clock() - t0)/CLOCKS_PER_SEC);
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
    dPdu(i) = (*mDPDepsInterp)(massDensity(i), specificThermalEnergy(i));      // We don't bother with the non-interpolation option for the derivs
    dPdrho(i) = (*mDPDRinterp)(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    specificHeat(i) = this->specificHeat(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    gamma(i) = this->gamma(massDensity(i), specificThermalEnergy(i));
  }
  // int n = massDensity.numElements();
  // if (n > 0) {
  //   Field<Dimension, Scalar> T("temperature", gamma.nodeList()),
  //                           cv("cv", gamma.nodeList());
  //   this->setTemperature(T, massDensity, specificThermalEnergy);
  //   this->setSpecificHeat(cv, massDensity, specificThermalEnergy);
  //   Scalar nDen;
  //   for (int i = 0; i != n; ++i) {
  //     nDen = massDensity(i)/mAtomicWeight;
  //     gamma(i) = 1.0 + mConstants.molarGasConstant()*nDen*safeInvVar(cv(i));
  //   }
  // }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    bulkModulus(i)=this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    entropy(i)=this->entropy(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  double P;
  if (mUseInterpolation) {
    P = (*mPinterp)(massDensity, specificThermalEnergy);
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    P *= mPconv;
  }
  return this->applyPressureLimits(P);
}

//------------------------------------------------------------------------------
// Calculate individual pressure and it's derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
std::tuple<typename Dimension::Scalar, typename Dimension::Scalar, typename Dimension::Scalar>
ANEOS<Dimension>::
pressureAndDerivs(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const {
  double P, dPdU, dPdR;
  if (mUseInterpolation) {
    P = (*mPinterp)(massDensity, specificThermalEnergy);
    dPdU = (*mDPDepsInterp)(massDensity, specificThermalEnergy);
    dPdR = (*mDPDRinterp)(massDensity, specificThermalEnergy);
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double eps, S, cV, dPdT, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &dPdT, &dPdR, &cs);
    P *= mPconv;
    dPdU = dPdT * mPconv/mEconv;
    dPdR *= mPconv/mRhoConv;
  }
  return std::make_tuple(this->applyPressureLimits(P), dPdU, dPdR);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  if (mUseInterpolation) {
    const auto Textra = Textrapolator(*mEpsMinInterp, *mEpsMaxInterp, *mTinterp);
    return Textra(massDensity, specificThermalEnergy);
  } else {
    const auto Feps = epsFunc(mMaterialNumber, mRhoConv, mTconv, mEconv);
    const auto Ftemp = Tfunc(mTmin, mTmax, *mEpsMinInterp, *mEpsMaxInterp, *mEpsInterp, Feps);
    return Ftemp(massDensity, specificThermalEnergy);
  }
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  if (mUseInterpolation) {
    return (*mEpsInterp)(massDensity, temperature);
  } else {
    auto T = temperature/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return eps * mEconv;
  }    
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  if (mUseInterpolation) {
    return max(1.0e-20, (*mCVinterp)(massDensity, temperature));
  } else {
    auto T = temperature/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return max(1.0e-20, cV * mCVconv);
  }    
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  if (mUseInterpolation) {
    return (*mCSinterp)(massDensity, specificThermalEnergy);
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return cs * mVelConv;
  }    
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  const auto Ti = (*mTinterp)(massDensity, specificThermalEnergy);     // this->temperature(massDensity, specificThermalEnergy);
  const auto cvi = this->specificHeat(massDensity, Ti);
  // return 1.0 + mConstants.molarGasConstant()*safeInvVar(cvi);
  const auto nDen = massDensity/mAtomicWeight;
  return 1.0 + mConstants.molarGasConstant()*nDen*safeInv(cvi, 1.0e-10);
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  if (mUseInterpolation) {
    const auto eps0 = (*mEpsMinInterp)(massDensity);
    const auto eps1 = (*mEpsMaxInterp)(massDensity);
    double K;
    if (specificThermalEnergy <= eps0) {
      K = (*mKinterp)(massDensity, eps0);
    } else if (specificThermalEnergy >= eps1) {
      K = (*mKinterp)(massDensity, eps1);
    } else {
      K = (*mKinterp)(massDensity, specificThermalEnergy);
    }
    return std::abs(K);
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return  std::abs(rho*DPDR) * mPconv;
  }
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  if (mUseInterpolation) {
    return max(1.0e-20, (*mSinterp)(massDensity, specificThermalEnergy));
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return max(1.0e-20, S * mSconv);
  }    
}

//------------------------------------------------------------------------------
// valid
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ANEOS<Dimension>::
valid() const {
  return true;
}

//------------------------------------------------------------------------------
// Access the attributes.
//------------------------------------------------------------------------------
template<typename Dimension>
int
ANEOS<Dimension>::
materialNumber() const {
  return mMaterialNumber;
}

template<typename Dimension>
unsigned
ANEOS<Dimension>::
numRhoVals() const {
  return mNumRhoVals;
}

template<typename Dimension>
unsigned
ANEOS<Dimension>::
numTvals() const {
  return mNumTvals;
}

template<typename Dimension>
double
ANEOS<Dimension>::
rhoMin() const {
  return mRhoMin;
}

template<typename Dimension>
double
ANEOS<Dimension>::
rhoMax() const {
  return mRhoMax;
}

template<typename Dimension>
double
ANEOS<Dimension>::
Tmin() const {
  return mTmin;
}

template<typename Dimension>
double
ANEOS<Dimension>::
Tmax() const {
  return mTmax;
}

template<typename Dimension>
double
ANEOS<Dimension>::
epsMin() const {
  return mEpsMin;
}

template<typename Dimension>
double
ANEOS<Dimension>::
epsMax() const {
  return mEpsMax;
}

template<typename Dimension>
bool
ANEOS<Dimension>::
useInterpolation() const {
  return mUseInterpolation;
}

//------------------------------------------------------------------------------
// Atomic weight.
//------------------------------------------------------------------------------
template<typename Dimension>
double
ANEOS<Dimension>::
atomicWeight() const {
  return mAtomicWeight;
}

}

