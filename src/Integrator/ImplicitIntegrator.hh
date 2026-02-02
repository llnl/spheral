//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrator -- Abstract base class for implicit in time integrators.
//
// Created by JMO, Fri Oct 18 16:34:11 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_ImplicitIntegrator__
#define __Spheral_ImplicitIntegrator__

#include "Integrator/Integrator.hh"

namespace Spheral {

template<typename Dimension>
class ImplicitIntegrator: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PackageIterator = typename std::vector<Physics<Dimension>*>::iterator;
  using ConstPackageIterator = typename std::vector<Physics<Dimension>*>::const_iterator;

  using BoundaryIterator = typename std::vector<Boundary<Dimension>*>::iterator;
  using ConstBoundaryIterator = typename std::vector<Boundary<Dimension>*>::const_iterator;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;

  // Constructors, destructor
  ImplicitIntegrator(DataBase<Dimension>& dataBase,
                     const std::vector<Physics<Dimension>*>& physicsPackages,
                     const Scalar tol = 1.0e-6);
  ImplicitIntegrator& operator=(const ImplicitIntegrator& rhs) = default;
  virtual ~ImplicitIntegrator() = default;

  // Override the step method for our implicit approach
  virtual bool step(Scalar maxTime) override;

  // Provide a method of looping over the physics packages and picking a
  // time step.
  virtual Scalar selectDt(const Scalar dtMin, 
                          const Scalar dtMax,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs) const override;

  // Find the maximum residual difference in the states
  virtual Scalar computeResiduals(const State<Dimension>& state1,
                                  const State<Dimension>& state0,
                                  const bool forceVerbose = false) const;

  // Internal data
  Scalar convergenceTolerance()                  const { return mTol; }
  Scalar maxAllowedDtMultiplier()                const { return mMaxAllowedDtMultiplier; }
  Scalar maxGoodDtMultiplier()                   const { return mMaxGoodDtMultiplier; }
  size_t numExplicitSteps()               const { return mNumExplicitSteps; }
  size_t numImplicitSteps()               const { return mNumImplicitSteps; }

  void convergenceTolerance(const Scalar x)            { mTol = x; }
  void maxAllowedDtMultiplier(const Scalar x)          { mMaxAllowedDtMultiplier = x; }

  //****************************************************************************
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const std::string& pathName)          const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName)       override;
  //****************************************************************************

  // Forbidden methods
  ImplicitIntegrator() = delete;

  using Integrator<Dimension>::step;
  using Integrator<Dimension>::mDtMultiplier;

protected:
  //-------------------------- Protected Interface --------------------------//
  // Override the package dt method to call the implicit version
  virtual TimeStepType dt(const Physics<Dimension>* pkg,
                          const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  Scalar mTol, mMaxAllowedDtMultiplier, mMaxGoodDtMultiplier;
  size_t mNumExplicitSteps, mNumImplicitSteps;
};

}

#endif
