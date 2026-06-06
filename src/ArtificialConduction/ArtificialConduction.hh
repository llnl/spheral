//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//

#ifndef ArtificialConduction_HH
#define ArtificialConduction_HH

#include "Physics/Physics.hh"
#include "RK/RKCorrectionParams.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class ArtificialConduction: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors
  ArtificialConduction(const TableKernel<Dimension>& W,
                       const Scalar alphaArCond, const RKOrder ACcorrectionOrder = RKOrder::LinearOrder);

  // Destructor
  virtual ~ArtificialConduction();

  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  //Allow access to the AC correction order.
  RKOrder ACcorrectionOrder() const;
  void ACcorrectionOrder(RKOrder order);

  // Provide default methods for registering and iterating derivatives.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  virtual std::string label() const override { return "Artificial Conduction"; }

  // Accessor Fns

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mKernel;

  // Our derivative field(s).
  FieldList<Dimension, Vector> mGradP;
  FieldList<Dimension, Scalar> mDepsDtArty;
  FieldList<Dimension, Scalar> mVsigMax;
  Scalar mAlphaArCond;
  RKOrder mACcorrectionOrder;

};

}

#endif
