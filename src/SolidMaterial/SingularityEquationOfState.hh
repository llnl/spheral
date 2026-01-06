//---------------------------------Spheral++----------------------------------//
// SingularityEquationOfState
//
// The base class for Singularity EOS
//----------------------------------------------------------------------------//
#ifndef __Spheral_SingularityEquationOfState_hh__
#define __Spheral_SingularityEquationOfState_hh__

#include "Material/UnitConversion.hh"
#include "SolidEquationOfState.hh"

#include <singularity-eos/eos/eos.hpp>

namespace Spheral {

using SingularityEOSType =
  singularity::Variant<singularity::IdealGas,
                       singularity::Vinet,
                       singularity::SpinerEOSDependsRhoSie,
                       singularity::SpinerEOSDependsRhoT,
                       singularity::EOSPAC>;

enum SingularityLimitLevel {
  NONE = 0,
  CS = 1,
  ALL = 2,
};

template<typename Dimension>
class SingularityEquationOfState: public SolidEquationOfState<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructor
  SingularityEquationOfState(SingularityEOSType eos,
                             const PhysicalConstants& constants,
                             const SingularityLimitLevel level,
                             const double referenceDensity,
                             const double etamin,
                             const double etamax,
                             const double minimumPressure,
                             const double maximumPressure,
                             const double minimumPressureDamage,
                             const MaterialPressureMinType minPressureType,
                             const double externalPressure);
  // Destructor
  ~SingularityEquationOfState();

  // We require any equation of state to define the following properties.
  virtual void setPressure(Field<Dimension, Scalar>& pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setPressureAndDerivs(Field<Dimension, Scalar>& pressure,
                                    Field<Dimension, Scalar>& dPdu,
                                    Field<Dimension, Scalar>& dPdrho,
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

  virtual bool valid() const override;

protected:

  // Temporary until we figure out the table issues
  double rhoLim(const Scalar massDensity) const;
  double eLim(const Scalar massDensity,
              const Scalar specificThermalEnergy,
              const SingularityLimitLevel level) const;
  double tLim(const Scalar massDensity,
              const Scalar temperature,
              const SingularityLimitLevel level) const;
  double csLim(const Scalar soundSpeed,
               SingularityLimitLevel level) const;
  double bmLim(const Scalar massDensity,
               const Scalar bulkModulus,
               const SingularityLimitLevel level) const;
  
private:
  // Equation of state
  SingularityEOSType mSEOS;

  // Limit level
  SingularityLimitLevel mLevel;
  
  // Conversion classes
  PhysicalConstants mCGS;
  UnitConversion mToCGS;
  UnitConversion mFromCGS;

  // Version 1.8 of Singularity doesn't have robust inverses
  static constexpr double mVerySmol = 100 * std::numeric_limits<double>::min();

  // Until we figure out what's going on with the tables
  double mCsMin, mCsMax, mRhoMin, mRhoMax;
  std::vector<double> mDensitiesCGS, mMinEnergiesCGS, mMinTemperaturesCGS;
};

}

#include "SingularityEquationOfStateInline.hh"

#endif
