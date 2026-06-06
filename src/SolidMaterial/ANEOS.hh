//---------------------------------Spheral++----------------------------------//
// ANEOS -- An interface to the ANEOS equation of state from 
// Melosh amony others.  This is a C++ wrapper around calls to the underlying
// ANEOS fortran library.  The user must provide the ancillary file containing
// parameters for ANEOS in its expected format.  Also make sure to keep units
// consistent in this separate input file with what you give the C++ EOS here!
//
// Created by JMO, Tue Apr 23 14:55:28 PDT 2013
//----------------------------------------------------------------------------//
#ifndef ANEOS_HH
#define ANEOS_HH

#include "SolidMaterial/SolidEquationOfState.hh"
#include "Utilities/CubicHermiteInterpolator.hh"
#include "Utilities/BiCubicInterpolator.hh"

#include <memory>
#include <tuple>

// Forward declarations.
namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class ANEOS: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
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
        const bool useInterpolation);
  ANEOS(const ANEOS& rhs) = default;
  ~ANEOS() = default;

  // We require any equation of state to define the following properties.
  virtual void setPressure(Field<Dimension, Scalar>& Pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,           // set pressure
                                    Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                                    Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                                    const Field<Dimension, Scalar>& massDensity,
                                    const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const override;

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const override;

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  // We also want the equivalent functions for individual calculations.
  Scalar pressure(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const;

  std::tuple<Scalar, Scalar, Scalar> pressureAndDerivs(const Scalar massDensity,
                                                       const Scalar specificThermalEnergy) const;

  Scalar temperature(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar specificThermalEnergy(const Scalar massDensity,
                               const Scalar temperature) const;

  Scalar specificHeat(const Scalar massDensity,
                      const Scalar temperature) const;

  Scalar soundSpeed(const Scalar massDensity,
                    const Scalar specificThermalEnergy) const;

  // Get the effective gamma (ratio of specific heats) for this eos.
  Scalar gamma(const Scalar massDensity,
               const Scalar specificThermalEnergy) const;

  // Get the bulk modulus.
  Scalar bulkModulus(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar entropy(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const;

  // The valid method.
  virtual bool valid() const override;

  // Access local variables used to lookup eps based on T.
  int materialNumber() const;
  unsigned numRhoVals() const;
  unsigned numTvals() const;
  double rhoMin() const;
  double rhoMax() const;
  double Tmin() const;
  double Tmax() const;
  double epsMin() const;
  double epsMax() const;
  bool useInterpolation() const;

  double atomicWeight() const;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mUseInterpolation;
  int mMaterialNumber;
  unsigned mNumRhoVals, mNumTvals;
  double mRhoMin, mRhoMax, mTmin, mTmax, mEpsMin, mEpsMax;

  using InterpolatorType = CubicHermiteInterpolator;
  using BiInterpolatorType = BiCubicInterpolator;
  std::shared_ptr<InterpolatorType> mEpsMinInterp, mEpsMaxInterp;
  std::shared_ptr<BiInterpolatorType> mEpsInterp, mTinterp, mPinterp, mCVinterp, mCSinterp, mKinterp, mSinterp, mDPDepsInterp, mDPDRinterp;

  // ANEOS internal units.
  PhysicalConstants mANEOSunits;

  // Units conversion from ANEOS.
  double mRhoConv, mTconv, mPconv, mEconv, mCVconv, mVelConv, mSconv;

  // Atomic weight.
  double mAtomicWeight;

  // Disallow default constructor
  ANEOS();

  using EquationOfState<Dimension>::mConstants;
};

}

#endif
