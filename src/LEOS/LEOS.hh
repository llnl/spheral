//---------------------------------Spheral++----------------------------------//
// LEOS
//
// Wrap up the Livermore Equation Of State (LEOS) package as an option in 
// Spheral.
//
// Created by JMO, Wed May  8 16:29:41 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_LEOS__
#define __Spheral_LEOS__

#include "SolidMaterial/SolidEquationOfState.hh"
#include "Utilities/DBC.hh"

#include <tuple>

namespace Spheral {

template<typename Dimension>
class LEOS: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
  LEOS(const int materialNumber,
       const PhysicalConstants& constants,
       const double externalPressure,
       const double minimumPressure,
       const double maximumPressure,
       const double minimumPressureDamage,
       const MaterialPressureMinType minPressureType,
       const std::string dbname = "leos",
       const std::string leosFileFormat = "",
       const double atomicWeight = 0.0);
  virtual ~LEOS() = default;

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

  // We can also provide the melt temperature using a similar Field interface.
  virtual void setMeltTemperature(Field<Dimension, Scalar>& meltTemperature,
                                  const Field<Dimension, Scalar>& massDensity,
                                  const Field<Dimension, Scalar>& specificThermalEnergy) const;

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

  Scalar gamma(const Scalar massDensity,
               const Scalar specificThermalEnergy) const;

  Scalar bulkModulus(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar entropy(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const;

  Scalar meltTemperature(const Scalar massDensity,
                         const Scalar specificThermalEnergy) const;

  // We override the base class reference density methods.
  using SolidEquationOfState<Dimension>::referenceDensity;
  virtual void referenceDensity(const double /*x*/) override { VERIFY2(false, "LEOS does not allow setting the reference density."); }

  // Use LEOS to do this reverse lookup for specific thermal energy as a function of pressure
  virtual Scalar specificThermalEnergyForPressure(const Scalar Ptarget,
                                                  const Scalar rho,
                                                  const Scalar epsMin,
                                                  const Scalar epsMax,
                                                  const Scalar epsTol,
                                                  const Scalar Ptol,
                                                  const unsigned maxIterations,
                                                  const bool verbose = false) const override;

  // Attributes.
  int materialNumber() const                                  { return mMaterialNumber; }
  std::string databaseName() const                            { return mDatabaseName; }
  double referenceTemperature() const;
  double atomicWeight() const                                 { return mAtomicWeight; }
  double K0() const                                           { return mK0; }
  std::string descriptor() const;
  
  // Forbidden methods
  LEOS() = delete;
  LEOS(const LEOS&) = delete;
  LEOS& operator=(const LEOS&) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  int mMaterialNumber;
  std::string mDatabaseName;
  double mRhoConv, mTconv, mPconv, mEconv, mCVconv, mVelConv, mSconv;
  double mAtomicWeight, mK0;

  // LEOS internal units.
  PhysicalConstants mLEOSunits;

  using EquationOfState<Dimension>::mConstants;
};

}

#endif
