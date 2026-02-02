//---------------------------------Spheral++----------------------------------//
// SPHRZ -- An SPH/ASPH hydrodynamic package for Spheral++,
//                   specialized for 2D RZ (cylindrical) geometry.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Fri May  6 16:18:36 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHRZ_hh__
#define __Spheral_SPHRZ_hh__

#include <string>
#include <memory>

#include "SPHBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension, typename Value, size_t numElements> class PairwiseField;

class SPHRZ: public SPHBase<Dim<2> > {

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
  SPHRZ(DataBase<Dimension>& dataBase,
        ArtificialViscosityHandle<Dimension>& Q,
        const TableKernel<Dimension>& W,
        const TableKernel<Dimension>& WPi,
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
        const Vector& xmin,
        const Vector& xmax);

  // No default constructor, copying, or assignment.
  SPHRZ() = delete;
  SPHRZ(const SPHRZ&) = delete;
  SPHRZ& operator=(const SPHRZ&) = delete;

  // Destructor.
  virtual ~SPHRZ() = default;

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

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SPHRZ" ; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
};

}

#endif
