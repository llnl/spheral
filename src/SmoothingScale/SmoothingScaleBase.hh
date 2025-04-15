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

#include <utility>
#include <cmath>

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

  // Constructors, destructor.
  explicit SmoothingScaleBase(const HEvolutionType HUpdate);
  SmoothingScaleBase() = delete;
  virtual ~SmoothingScaleBase() {};

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
  HEvolutionType HEvolution() const;
  void HEvolution(HEvolutionType type);

  // Our state fields
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, SymTensor>& DHDt() const;

  //****************************************************************************
  // Methods required for restarting (descendants still need to provide the required "label")
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  HEvolutionType mHEvolution;
  FieldList<Dimension, SymTensor> mHideal, mDHDt;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "SmoothingScaleBaseInline.hh"

#endif
