//---------------------------------Spheral++----------------------------------//
// LinearPolynomialEquationOfState -- An equation of state approximated by a
// linear polynomial, i.e.:
//
//   P(rho, e) = A0 + A1*mu + a2*mu^2 + a3*mu^3 + (B0 + B1*mu + B2*mu^2)*e
//   mu = rho/rho0 - 1.0
//
// Created by JMO, Thu May  5 16:07:36 PDT 2005
//----------------------------------------------------------------------------//
#ifndef LinearPolynomialEquationOfState_HH
#define LinearPolynomialEquationOfState_HH

#include "SolidEquationOfState.hh"

#include <tuple>
#include <float.h>

// Forward declarations.
namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class LinearPolynomialEquationOfState: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  LinearPolynomialEquationOfState(const double referenceDensity,
                                  const double etamin,
                                  const double etamax,
                                  const double a0,
                                  const double a1,
                                  const double a2,
                                  const double a3,
                                  const double b0,
                                  const double b1,
                                  const double b2,
                                  const double atomicWeight,
                                  const PhysicalConstants& constants,
                                  const double externalPressure,
                                  const double minimumPressure,
                                  const double maximumPressure,
                                  const double minimumPressureDamage,
                                  const MaterialPressureMinType minPressureType);
  ~LinearPolynomialEquationOfState();

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

  // Access the member data.
  double a0() const;
  double a1() const;
  double a2() const;
  double a3() const;
  double b0() const;
  double b1() const;
  double b2() const;
  double atomicWeight() const;
  
  void a0(double x);
  void a1(double x);
  void a2(double x);
  void a3(double x);
  void b0(double x);
  void b1(double x);
  void b2(double x);
  void atomicWeight(double x);
  
  // Compute the derivative of the pressure with respect to the density.
  double computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
  double mA0, mA1, mA2, mA3;
  double mB0, mB1, mB2;
  double mAtomicWeight;
  double mCv;
  double mGamma;

  // Disallow default constructor
  LinearPolynomialEquationOfState();

};

}

#include "LinearPolynomialEquationOfStateInline.hh"

#endif
