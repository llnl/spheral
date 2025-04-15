//---------------------------------Spheral++----------------------------------//
// ForwardEuler -- Advance the set of Physics packages in time using first
// order Forward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Tue Aug  1 14:44:27 PDT 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_ForwardEuler__
#define __Spheral_ForwardEuler__

#include "Integrator.hh"

namespace Spheral {

template<typename Dimension>
class ForwardEuler: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  ForwardEuler(DataBase<Dimension>& dataBase,
               const std::vector<Physics<Dimension>*>& physicsPackages);
  ForwardEuler& operator=(const ForwardEuler& rhs) = default;
  virtual ~ForwardEuler() = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Forbidden methods
  ForwardEuler() = delete;
};

}

#endif
