//---------------------------------Spheral++----------------------------------//
// SolidSPHRZ -- The axisymmetric (RZ) SPH/ASPH solid material
//                        hydrodynamic package for Spheral++.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Mon May  9 11:01:51 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSPHRZ_hh__
#define __Spheral_SolidSPHRZ_hh__

#include <memory>
#include <string>

#include "SolidSPH.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SolidSPHRZ: public SolidSPH<Dim<2>> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using SymTensor = Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 2u>;
  using ConstBoundaryIterator = Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SolidSPHRZ(DataBase<Dimension>& dataBase,
             ArtificialViscosity<Dimension>& Q,
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
  SolidSPHRZ() = delete;
  SolidSPHRZ(const SolidSPHRZ&) = delete;
  SolidSPHRZ& operator=(const SolidSPHRZ&) = delete;

  // Destructor.
  virtual ~SolidSPHRZ() = default;

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
                               const QType& Q) const;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Access our state.
  const PairAccelerationsType& pairAccelerations()        const { VERIFY2(mPairAccelerationsPtr, "SPH ERROR: pairAccelerations not initialized on access"); return *mPairAccelerationsPtr; }
  const FieldList<Dimension, Vector>& selfAccelerations() const { return mSelfAccelerations; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label()                    const override { return "SolidSPHRZ"; }

private:
  //--------------------------- Private Interface ---------------------------//
  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
  FieldList<Dimension, Vector> mSelfAccelerations;
};

}

#endif
