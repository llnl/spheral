//---------------------------------Spheral++----------------------------------//
// GruneisenEquationOfState -- Gruneisen  equation of state.
//
// Created by MLF & JMO Wed Sep 17 13:32:54 EDT 2003 
// Reference: Equation of State and Strength of Properties of Selected Materials
//            Daniel J. Steinberg, UCRL-MA-106439, February 13, 1991
//----------------------------------------------------------------------------//
#ifndef GruneisenEquationOfState_HH
#define GruneisenEquationOfState_HH

#include "SolidEquationOfState.hh"

#include <tuple>
#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class GruneisenEquationOfState: public SolidEquationOfState<Dimension> { 

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  GruneisenEquationOfState(const double referenceDensity,
                           const double etamin,
                           const double etamax,
                           const double C0, 
                           const double S1,
                           const double S2,
                           const double S3,
                           const double gamma0,
                           const double b,
                           const double atomicWeight,
                           const PhysicalConstants& constants,
                           const double externalPressure,
                           const double minimumPressure,
                           const double maximumPressure,
                           const double minimumPressureDamage,
                           const MaterialPressureMinType minPressureType);
  virtual ~GruneisenEquationOfState();

  // We require any equation of state to define the following methods for Fields.
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

  Scalar gamma(const Scalar massdensity,
               const Scalar specificThermalEnergy) const;

  Scalar bulkModulus(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar entropy(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const;

  // Access the member data.
  double C0() const;
  double S1() const;
  double S2() const;
  double S3() const;
  double gamma0() const;
  double b() const;
  double atomicWeight() const;
  double Cv() const;

  void C0(double val);
  void S1(double val);
  void S2(double val);
  void S3(double val);
  void gamma0(double val);
  void b(double val);
  void atomicWeight(double val);

  // Option to scale the thermal energy term by.  This is mostly useful for test problems
  // where you want to make the Gruneisen independent of energy.
  double energyMultiplier() const;
  void energyMultiplier(double val);

  // Compute (\partial P)/(\partial rho) for use in sound speed and bulk modulus.
  Scalar computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

  // Equations of state should have a valid test.
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
  double mC0;
  double mS1;
  double mS2;
  double mS3;
  double mgamma0;
  double mb;
  double mAtomicWeight;
  double mCv;
  double mEnergyMultiplier;

  // No default constructor, copying, or assignment.
  GruneisenEquationOfState();
  GruneisenEquationOfState(const GruneisenEquationOfState&);
  GruneisenEquationOfState& operator=(const GruneisenEquationOfState&);

  using EquationOfState<Dimension>::mConstants;
};

}

#include "GruneisenEquationOfStateInline.hh"

#endif
