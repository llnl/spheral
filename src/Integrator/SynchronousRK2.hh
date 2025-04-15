//---------------------------------Spheral++----------------------------------//
// SynchronousRK2 -- Advance the set of Physics packages in time using second
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun  1 22:26:45 PDT 2000
//----------------------------------------------------------------------------//
#ifndef SynchronousRK2_HH
#define SynchronousRK2_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class SynchronousRK2: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  SynchronousRK2(DataBase<Dimension>& dataBase,
                 const std::vector<Physics<Dimension>*>& physicsPackages);
  virtual ~SynchronousRK2() = default;
  SynchronousRK2& operator=(const SynchronousRK2& rhs) = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Fobidden methods
  SynchronousRK2() = delete;
};

}

#endif
