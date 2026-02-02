//---------------------------------Spheral++----------------------------------//
// PolytropicEquationOfState -- A polytrope equation of state.
//
// Created by JMO, Wed Apr  2 14:48:12 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_PolytropicEquationOfState_hh__
#define __Spheral_PolytropicEquationOfState_hh__

#include "EquationOfState.hh"

namespace Spheral {

template<typename Dimension>
class PolytropicEquationOfState: public EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  PolytropicEquationOfState(const double K,
                            const double index,
                            const double mu,
                            const PhysicalConstants& constants,
                            const double minimumPressure,
                            const double maximumPressure,
                            const MaterialPressureMinType minPressureType,
                            const double externalPressure);
  ~PolytropicEquationOfState();

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
  Scalar polytropicConstant() const;
  Scalar polytropicIndex() const;
  Scalar gamma() const;
  virtual Scalar molecularWeight() const override;
  
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mPolytropicConstant;
  Scalar mPolytropicIndex;
  Scalar mGamma;
  Scalar mGamma1;
  Scalar mMolecularWeight;

  // No default constructor, copying, or assignment.
  PolytropicEquationOfState();
  PolytropicEquationOfState(const PolytropicEquationOfState&);
  PolytropicEquationOfState& operator=(const PolytropicEquationOfState&);

  using EquationOfState<Dimension>::mConstants;
};

}

#include "PolytropicEquationOfStateInline.hh"

#endif
