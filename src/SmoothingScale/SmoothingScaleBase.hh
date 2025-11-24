//---------------------------------Spheral++----------------------------------//
// SmoothingScaleBase
//
// Abstract base class for packages that advance the smoothing scale.
//
// Created by JMO, Wed Sep 14 13:27:39 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_SmooothingScaleBase__
#define __Spheral_SmooothingScaleBase__

#include "Geometry/Dimension.hh"
#include "Physics/Physics.hh"
#include "Field/FieldList.hh"
#include "DataOutput/registerWithRestart.hh"
#include "SmoothingScale/ASPHSmoothingScaleUserFilter.hh"
#include "SmoothingScale/ASPHRadialFunctor.hh"

#include <utility>
#include <cmath>
#include <memory>   // std::shared_ptr

namespace Spheral {

class FileIO;

enum class HEvolutionType {
  IdealH = 0,
  IntegrateH = 1,
  FixedH = 2,
};

template<typename Dimension>
class SmoothingScaleBase: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;
  using HidealFilterType = ASPHSmoothingScaleUserFilter<Dimension>;
  using RadialFunctorType = ASPHRadialFunctor<Dimension>;

  // Constructors, destructor.
  SmoothingScaleBase(const HEvolutionType HUpdate,
                     const bool fixShape,
                     const bool radialOnly);
  SmoothingScaleBase() = delete;
  virtual ~SmoothingScaleBase() = default;

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Default smoothing scale methods to not constrain the time step
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override { return TimeStepType(1e100, "SmoothingScale -- no vote"); }

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Optional hook to be called at the beginning of a time step.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Return the maximum state change we care about for checking for convergence in the implicit integration methods.
  // We assume the default limting on position change is good enough for H, so don't add anything extra here.
  virtual ResidualType maxResidual(const DataBase<Dimension>& dataBase, 
                                   const State<Dimension>& state1,
                                   const State<Dimension>& state0,
                                   const Scalar tol) const override { return ResidualType(0.0, "SmoothingScale -- no vote"); }

  // Given the volume and target nperh, compute an effective target hmax
  Scalar hmax(const Scalar Vi, const Scalar nPerh) const;

  // Flag to select how we want to evolve the H tensor.
  // the continuity equation.
  HEvolutionType HEvolution() const                                               { return mHEvolution; }
  void HEvolution(HEvolutionType type)                                            { mHEvolution = type; }

  // Special evolution flags
  bool fixShape() const                                                           { return mFixShape; }
  bool radialOnly() const                                                         { return mRadialOnly; }
  void fixShape(const bool x)                                                     { mFixShape = x; }
  void radialOnly(const bool x)                                                   { mRadialOnly = x; }

  // Optional user functor to manipulate the final ideal H vote
  std::shared_ptr<HidealFilterType> HidealFilter() const                          { return mHidealFilterPtr; }
  void HidealFilter(std::shared_ptr<HidealFilterType> functorPtr)                 { mHidealFilterPtr = functorPtr; }

  // Optional user functor to override the radial unit normal and radius for radialOnly mode
  std::shared_ptr<RadialFunctorType> RadialFunctor() const                        { return mRadialFunctorPtr; }
  void RadialFunctor(std::shared_ptr<RadialFunctorType> functorPtr)               { mRadialFunctorPtr = functorPtr; }

  // Exponent for damping smoothing scale iterations
  double damping() const                                                          { return mDamping; }
  void dampling(const double d)                                                   { mDamping = d; }

  // Our state fields
  const FieldList<Dimension, SymTensor>& Hideal() const                           { return mHideal; }
  const FieldList<Dimension, SymTensor>& DHDt() const                             { return mDHDt; }
  const FieldList<Dimension, Scalar>& radius0() const                             { return mRadius0; }

  //****************************************************************************
  // Methods required for restarting (descendants still need to provide the required "label")
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mFixShape, mRadialOnly;
  double mDamping;
  HEvolutionType mHEvolution;
  FieldList<Dimension, SymTensor> mHideal, mDHDt;
  FieldList<Dimension, Scalar> mRadius0;
  std::shared_ptr<HidealFilterType> mHidealFilterPtr;
  std::shared_ptr<RadialFunctorType> mRadialFunctorPtr;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "SmoothingScaleBaseInline.hh"

#endif
