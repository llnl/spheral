//---------------------------------Spheral++----------------------------------//
// OsborneEquationOfState -- Osborne  equation of state.
//
// Reference: PAGOSA Physics manual, LA-14425-M
//----------------------------------------------------------------------------//
#ifndef OsborneEquationOfState_HH
#define OsborneEquationOfState_HH

#include "SolidEquationOfState.hh"

#include <limits>

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class OsborneEquationOfState: public SolidEquationOfState<Dimension> { 

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  OsborneEquationOfState(const double referenceDensity,
                         const double etamin,
                         const double etamax,
                         const double a1,
                         const double a2pos,
                         const double a2neg,
                         const double b0,
                         const double b1,
                         const double b2pos,
                         const double b2neg,
                         const double c0,
                         const double c1,
                         const double c2pos,
                         const double c2neg,
                         const double E0,
                         const double atomicWeight,
                         const PhysicalConstants& constants,
                         const double externalPressure,
                         const double minimumPressure,
                         const double maximumPressure,
                         const double minimumPressureDamage,
                         const MaterialPressureMinType minPressureType);
  virtual ~OsborneEquationOfState();

  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(Field<Dimension, Scalar>& pressure,
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

  // Access the member data.
  double a1() const;
  double a2pos() const;
  double a2neg() const;
  double b0() const;
  double b1() const;
  double b2pos() const;
  double b2neg() const;
  double c0() const;
  double c1() const;
  double c2pos() const;
  double c2neg() const;
  double E0() const;
  double atomicWeight() const;
  double Cv() const;

  void a1(double val);
  void a2pos(double val);
  void a2neg(double val);
  void b0(double val);
  void b1(double val);
  void b2pos(double val);
  void b2neg(double val);
  void c0(double val);
  void c1(double val);
  void c2pos(double val);
  void c2neg(double val);
  void E0(double val);
  void atomicWeight(double val);

  // Compute an individual value for DPDrho.
  double DPDrho(const double massDensity,
                const double specificThermalEnergy) const;

  // Equations of state should have a valid test.
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
  double mA1, mA2pos, mA2neg, mB0, mB1, mB2pos, mB2neg, mC0, mC1, mC2pos, mC2neg, mE0;
  double mAtomicWeight;
  double mCv;

  // No default constructor, copying, or assignment.
  OsborneEquationOfState();
  OsborneEquationOfState(const OsborneEquationOfState&);
  OsborneEquationOfState& operator=(const OsborneEquationOfState&);

  using EquationOfState<Dimension>::mConstants;
};

}

#include "OsborneEquationOfStateInline.hh"

#endif
