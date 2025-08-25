//---------------------------------Spheral++----------------------------------//
// PSPHHydroBase -- The PSPH/APSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Wed Dec 16 20:52:02 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_PSPHHydroBase_hh__
#define __Spheral_PSPHHydroBase_hh__

#include <string>

#include "SPH.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
class FileIO;

template<typename Dimension>
class PSPH: public SPH<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  PSPH(DataBase<Dimension>& dataBase,
       ArtificialViscosity<Dimension>& Q,
       const TableKernel<Dimension>& W,
       const TableKernel<Dimension>& WPi,
       const double cfl,
       const bool useVelocityMagnitudeForDt,
       const bool compatibleEnergyEvolution,
       const bool evolveTotalEnergy,
       const bool XSPH,
       const bool correctVelocityGradient,
       const bool HopkinsConductivity,
       const bool sumMassDensityOverAllNodeLists,
       const MassDensityType densityUpdate,
       const Vector& xmin,
       const Vector& xmax);

  // No default constructor, copying, or assignment.
  PSPH() = delete;
  PSPH(const PSPH&) = delete;
  PSPH& operator=(const PSPH&) = delete;

  // Destructor.
  virtual ~PSPH() = default;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual
  void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                            State<Dimension>& state,
                                            StateDerivatives<Dimension>& derivs) override;

  // Pre-step initializations.
  virtual 
  void preStepInitialize(const DataBase<Dimension>& dataBase, 
                         State<Dimension>& state,
                         StateDerivatives<Dimension>& derivatives) override;

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

  // Post-state update. For PSPH this is where we recompute the PSPH pressure and corrections.
  virtual 
  bool postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Flag determining if we're applying Hopkins 2014 conductivity.
  bool HopkinsConductivity()                                        const { return mHopkinsConductivity; }
  void HopkinsConductivity(bool x)                                        { mHopkinsConductivity = x; }

  // The state field lists we're maintaining.
  const FieldList<Dimension, Scalar>&    gamma()                    const { return mGamma; }
  const FieldList<Dimension, Scalar>&    PSPHcorrection()           const { return mPSPHcorrection; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label()                                       const override { return "PSPH"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  bool mHopkinsConductivity;

  //PSPH Fields
  FieldList<Dimension, Scalar>    mGamma;
  FieldList<Dimension, Scalar>    mPSPHcorrection;
};

}

#endif
