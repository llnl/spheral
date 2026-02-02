//---------------------------------Spheral++----------------------------------//
// Second order predictor corrector time integrator.
// Based on the predictor corrector scheme described in Monaghans tensile
// instability paper (Monaghan 2000, JCP 159, 290.)
//
// Created by JMO, Wed Dec  4 21:39:36 PST 2002
//----------------------------------------------------------------------------//
#ifndef PredictorCorrector_HH
#define PredictorCorrector_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class PredictorCorrector: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  PredictorCorrector(DataBase<Dimension>& dataBase,
                     const std::vector<Physics<Dimension>*>& physicsPackages);
  virtual ~PredictorCorrector() = default;
  PredictorCorrector& operator=(const PredictorCorrector& rhs) = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Forbidden methods
  PredictorCorrector() = delete;
};

}

#endif
