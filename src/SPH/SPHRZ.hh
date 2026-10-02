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
#include "chai/managed_ptr.hpp"

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
  using PairWorkType = PairwiseField<Dimension, Scalar, 2u>;
  using ConstBoundaryIterator = Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SPHRZ(DataBase<Dimension>& dataBase,
        ArtificialViscosity<Dimension>& Q,
        const TableKernel<Dimension>& W,
        const TableKernel<Dimension>& WPi,
        const double cfl,
        const bool useVelocityMagnitudeForDt,
        const bool useNewAccelerationMagnitudeForDt,
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

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual
  void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
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
                               chai::managed_ptr<QType>& Q) const;

  // Provide a hook to be called after all physics packages have had their
  // evaluateDerivatives method called, but before anyone does anything
  // with those derivatives.
  virtual 
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

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
  const PairWorkType& pairWork()                          const { VERIFY2(mPairWorkPtr, "SPH ERROR: pairWork not initialized on access"); return *mPairWorkPtr; }
  const FieldList<Dimension, Scalar>& DmassDensityDtRZ()  const { return mDmassDensityDtRZ; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SPHRZ" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

  // Integrate an estimate of vr/r over a timestep given the initial radial
  // velocity, position, acceleration, smoothing scale, and timestep.
  static double integrate_vr_over_r(double vr,
                                    double r,
                                    double ar,
                                    const double hr,
                                    const double dt);

private:
  //--------------------------- Private Interface ---------------------------//
  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
  std::unique_ptr<PairWorkType> mPairWorkPtr;

  FieldList<Dimension, Scalar> mDmassDensityDtRZ;
};

}

#endif
