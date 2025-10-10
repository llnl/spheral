//---------------------------------Spheral++----------------------------------//
// Geodyn -- An interface to the Geodyn equation of state.
//
// This class implements three distinct Spheral interfaces:
//   EquationOfState : provides the ordinary Spheral EOS interface
//   StrengthModel   : also answers the Strength questions about shear modulus 
//                     and yield strength
//   Physics         : Geodyn also needs to advance it's own internal state, so
//                     this class implements the Physics interface to support
//                     that.
//
// Created by JMO, Wed Jun  3 22:46:29 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_Geodyn_hh__
#define __Spheral_Geodyn_hh__

#include "SolidMaterial/SolidEquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Physics/Physics.hh"

// Forward declarations.
namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class Geodyn:
    public Physics<Dimension>,
    public SolidEquationOfState<Dimension>,
    public StrengthModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors, destructors.
  // We should add arguments to the constructor to specify the Geodyn model
  Geodyn(const PhysicalConstants& constants,
         const double minimumPressure,
         const double maximumPressure,
         const MaterialPressureMinType minPressureType);
  ~Geodyn();

  // .............................. EquationOfState interface ..............................
  virtual void setPressure(Field<Dimension, Scalar>& Pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const;

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const;

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual bool valid() const;

  // .............................. StrengthModel interface ..............................
  virtual void shearModulus(Field<Dimension, Scalar>& shearModulus,
                            const Field<Dimension, Scalar>& density,
                            const Field<Dimension, Scalar>& specificThermalEnergy,
                            const Field<Dimension, Scalar>& pressure) const;

  virtual void yieldStrength(Field<Dimension, Scalar>& yieldStrength,
                             const Field<Dimension, Scalar>& density,
                             const Field<Dimension, Scalar>& specificThermalEnergy,
                             const Field<Dimension, Scalar>& pressure,
                             const Field<Dimension, Scalar>& plasticStrain,
                             const Field<Dimension, Scalar>& plasticStrainRate) const;

  virtual void soundSpeed(Field<Dimension, Scalar>& soundSpeed,
                          const Field<Dimension, Scalar>& density,
                          const Field<Dimension, Scalar>& specificThermalEnergy,
                          const Field<Dimension, Scalar>& pressure,
                          const Field<Dimension, Scalar>& fluidSoundSpeed) const;

  // .............................. Physics interface ..............................
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);

  // It's useful to have labels for Physics packages.  We'll require this to have
  // the same signature as the restart label.
  virtual std::string label() const { return "Geodyn Spheral interface"; }

private:
  //--------------------------- Private Interface ---------------------------//
  // Tables for the temp->energy lookup.
  std::vector<double> mGeodynState;

  // GEODYN internal units.
  PhysicalConstants mGeodynUnits;

  // Units conversion from Geodyn.
  double mRhoConv, mTconv, mPconv, mEconv, mCVconv, mVelConv;

  // Disallow default constructor
  Geodyn();

  using EquationOfState<Dimension>::mConstants;
};

}

#endif
