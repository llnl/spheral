//---------------------------------Spheral++----------------------------------//
// SolidSPH -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSPH_hh__
#define __Spheral_SolidSPH_hh__

#include <float.h>
#include <string>

#include "SPH/SPH.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosityHandle;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
class FileIO;

template<typename Dimension>
class SolidSPH: public SPH<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SolidSPH(DataBase<Dimension>& dataBase,
           ArtificialViscosityHandle<Dimension>& Q,
           const TableKernel<Dimension>& W,
           const TableKernel<Dimension>& WPi,
           const TableKernel<Dimension>& WGrad,
           const double cfl,
           const bool useVelocityMagnitudeForDt,
           const bool compatibleEnergyEvolution,
           const bool evolveTotalEnergy,
           const bool gradhCorrection,
           const bool XSPH,
           const bool correctVelocityGradient,
           const bool sumMassDensityOverAllNodeLists,
           const MassDensityType densityUpdate,
           const double epsTensile,
           const double nTensile,
           const bool damageRelieveRubble,
           const bool strengthInDamage,
           const Vector& xmin,
           const Vector& xmax);

  // No default constructor, copying, or assignment.
  SolidSPH() = delete;
  SolidSPH(const SolidSPH&) = delete;
  SolidSPH& operator=(const SolidSPH&) = delete;

  // Destructor.
  virtual ~SolidSPH() = default;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;
  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;
  template<typename QType>
  void evaluateDerivativesImpl(const Scalar time,
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase,
                               const State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives,
                               const QType& Q) const;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Gradient kernel
  const TableKernel<Dimension>& GradKernel()                   const { return mGradKernel; }

  // The state field lists we're maintaining.
  const FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const { return mDdeviatoricStressDt; }
  const FieldList<Dimension, Scalar>& bulkModulus()            const { return mBulkModulus; }
  const FieldList<Dimension, Scalar>& shearModulus()           const { return mShearModulus; }
  const FieldList<Dimension, Scalar>& yieldStrength()          const { return mYieldStrength; }
  const FieldList<Dimension, Scalar>& plasticStrain0()         const { return mPlasticStrain0; }

  // Control whether allow damaged material to have stress relieved.
  bool damageRelieveRubble()                                   const { return mDamageRelieveRubble; }
  void damageRelieveRubble(bool x)                                   { mDamageRelieveRubble = x; }

  // Do we allow damaged material to have strength?
  bool strengthInDamage()                                      const { return mStrengthInDamage; }
  void strengthInDamage(bool x)                                      { mStrengthInDamage = x; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label()                                 const override { return "SolidSPH"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mDamageRelieveRubble, mStrengthInDamage;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mGradKernel;   // Gradient kernel

  // Some internal scratch fields.
  FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldList<Dimension, Scalar> mBulkModulus;
  FieldList<Dimension, Scalar> mShearModulus;
  FieldList<Dimension, Scalar> mYieldStrength;
  FieldList<Dimension, Scalar> mPlasticStrain0;

};

}

#endif
