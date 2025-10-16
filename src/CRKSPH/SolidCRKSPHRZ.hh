//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHRZ -- The CRKSPH/ACRKSPH solid material hydrodynamic
// package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Fri May 13 10:50:36 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidCRKSPHRZ_hh__
#define __Spheral_SolidCRKSPHRZ_hh__

#include "CRKSPH/SolidCRKSPH.hh"

#include <memory>
#include <string>

namespace Spheral {

class SolidCRKSPHRZ: public SolidCRKSPH<Dim<2>> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using SymTensor = Dimension::SymTensor;
  using ThirdRankTensor = Dimension::ThirdRankTensor;
  using FourthRankTensor = Dimension::FourthRankTensor;
  using FifthRankTensor = Dimension::FifthRankTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 2u>;
  using ConstBoundaryIterator = Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHRZ(DataBase<Dimension>& dataBase,
                ArtificialViscosity<Dimension>& Q,
                const RKOrder order,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool XSPH,
                const MassDensityType densityUpdate,
                const double epsTensile,
                const double nTensile,
                const bool damageRelieveRubble);

  // No default constructor, copying, or assignment.
  SolidCRKSPHRZ() = delete;
  SolidCRKSPHRZ(const SolidCRKSPHRZ&) = delete;
  SolidCRKSPHRZ& operator=(const SolidCRKSPHRZ&) = delete;

  // Destructor.
  virtual ~SolidCRKSPHRZ() = default;

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

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
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
                               chai::managed_ptr<QType> Q) const;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Access our state.
  const PairAccelerationsType& pairAccelerations()        const { VERIFY2(mPairAccelerationsPtr, "SolidCRKSPHRZ ERROR: pairAccelerations not initialized on access"); return *mPairAccelerationsPtr; }
  const FieldList<Dimension, Vector>& selfAccelerations() const { return mSelfAccelerations; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidCRKSPHRZ"; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
  FieldList<Dimension, Vector> mSelfAccelerations;
};

}

#endif
