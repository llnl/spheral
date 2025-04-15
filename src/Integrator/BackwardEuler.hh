//---------------------------------Spheral++----------------------------------//
// BackwardEuler -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_BackwardEuler__
#define __Spheral_BackwardEuler__

#include "Integrator/ImplicitIntegrator.hh"
#include <vector>

namespace Spheral {

template<typename Dimension>
class BackwardEuler: public ImplicitIntegrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  BackwardEuler(DataBase<Dimension>& dataBase,
                const std::vector<Physics<Dimension>*> physicsPackages = std::vector<Physics<Dimension>*>(),
                const Scalar beta = 1.0,
                const Scalar ftol = 1.0e-8,
                const Scalar tol = 1.0e-8,     // scsteptol in KINSOL
                const size_t maxIterations = 100u);
  BackwardEuler& operator=(const BackwardEuler& rhs) = default;
  virtual ~BackwardEuler() = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // Access internal state
  Scalar beta()                           const { return mBeta; }
  size_t maxIterations()                  const { return mMaxIters; }
  Scalar ftol()                           const { return mftol; }
  Scalar tM2()                            const { return mtM2; }
  Scalar tM1()                            const { return mtM1; }
  const std::vector<double>& solutionM2() const { return mSolutionM2; }
  const std::vector<double>& solutionM1() const { return mSolutionM1; }

  void beta(const Scalar x)                     { mBeta = x; }
  void maxIterations(const size_t x)            { mMaxIters = x; }
  void ftol(const Scalar x)                     { mftol = x; }

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  //--------------------------- Public Interface ---------------------------//
private:
  Scalar mBeta, mftol, mtM2, mtM1;
  size_t mMaxIters;
  std::vector<double> mSolutionM2, mSolutionM1;

  using Integrator<Dimension>::mDtMultiplier;
  using ImplicitIntegrator<Dimension>::mNumExplicitSteps;
  using ImplicitIntegrator<Dimension>::mNumImplicitSteps;
};

}

#endif
