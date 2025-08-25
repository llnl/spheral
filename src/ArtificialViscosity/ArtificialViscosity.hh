//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- A base class for ArtficialViscosity that strips
// off the QPiType template parameter.  This makes a convenient way to break
// that template parameter from spreading into classes that need to consume an
// ArtificialViscosity.
//
// Created by JMO, Fri Dec 13 10:06:12 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosity__
#define __Spheral_ArtificialViscosity__

#include "Physics/Physics.hh"
#include "Field/FieldList.hh"
#include "DataOutput/registerWithRestart.hh"
#include "Utilities/SpheralMessage.hh"

#include <utility>
#include <typeindex>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class Boundary;
class FileIO;

template<typename Dimension>
class ArtificialViscosity: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;

  // Constructors, destructor
  ArtificialViscosity(const Scalar Clinear,
                            const Scalar Cquadratic,
                            const TableKernel<Dimension>& kernel);
  virtual ~ArtificialViscosity() = default;

  // No default constructor, copying, or assignment
  ArtificialViscosity() = delete;
  ArtificialViscosity(const ArtificialViscosity&) = delete;
  ArtificialViscosity& operator=(const ArtificialViscosity&) = delete;

  //...........................................................................
  // Virtual methods we expect ArtificialViscosities to provide
  // Require ArtificialViscosities to specify the type_index of the descendant QPiType
  virtual std::type_index QPiTypeIndex() const = 0;

  // Some AVs need the velocity gradient computed, so they should override this to true
  virtual bool requireVelocityGradient()                                  const { return false; }

  // Update the locally stored velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs);

  //...........................................................................
  // Standard Physics package methods
  // Most ArtificialViscosities will not have an evaluateDerivatives, so by default no-op this
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override {};

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  // initialize and post-state update are our chances to update the velocity gradient if needed
  virtual bool initialize(const Scalar time, 
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase, 
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivatives) override;

  virtual bool postStateUpdate(const Scalar time, 
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase, 
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives) override;

  // Return the maximum state change we care about for checking for convergence in the implicit integration methods.
  // Artificial viscosities default to just relying on the hydro to handle this, so return no vote.
  virtual ResidualType maxResidual(const DataBase<Dimension>& dataBase, 
                                   const State<Dimension>& state1,
                                   const State<Dimension>& state0,
                                   const Scalar tol) const override { return std::make_pair<double, std::string>(0.0, this->label() + " no vote"); }

  //...........................................................................
  // Methods
  // Calculate the curl of the velocity given the stress tensor.
  Scalar curlVelocityMagnitude(const Tensor& DvDx) const;

  // Find the Balsara shear correction multiplier
  Scalar calcBalsaraShearCorrection(const Tensor& DvDx,
                                    const SymTensor& H,
                                    const Scalar& cs) const;

  // Access stored state
  Scalar                              Cl()                                const { return mClinear; }
  Scalar                              Cq()                                const { return mCquadratic; }
  bool                                balsaraShearCorrection()            const { return mBalsaraShearCorrection; }
  Scalar epsilon2()                                                       const { return mEpsilon2; }
  Scalar negligibleSoundSpeed()                                           const { return mNegligibleSoundSpeed; }
  const FieldList<Dimension, Scalar>& maxViscousPressure()                const { return mMaxViscousPressure; }
  const FieldList<Dimension, Scalar>& effViscousPressure()                const { return mEffViscousPressure; }
  const FieldList<Dimension, Tensor>& DvDx()                              const { return mDvDx; }

  bool                                rigorousVelocityGradient()          const { return mRigorousVelocityGradient; }
  const TableKernel<Dimension>&       kernel()                            const { return mWT; }

  void Cl(Scalar x)                                                             { mClinear = x; }
  void Cq(Scalar x)                                                             { mCquadratic = x; }
  void balsaraShearCorrection(bool x)                                           { mBalsaraShearCorrection = x; }
  void epsilon2(Scalar x)                                                       { mEpsilon2 = x; }
  void negligibleSoundSpeed(Scalar x)                                           { REQUIRE(x > 0.0); mNegligibleSoundSpeed = x; }
  void rigorousVelocityGradient(bool x)                                         { mRigorousVelocityGradient = x; }

  // Deprecated options
  bool                               limiter()                            const { DeprecationWarning("ArtificialViscosity::limiter"); return false; }
  void                               limiter(const bool x)                      { DeprecationWarning("ArtificialViscosity::limiter"); }

  //...........................................................................
  // Methods required for restarting.
  virtual std::string label()                                    const override { return "ArtificialViscosity"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

protected:
  //--------------------------- Protected Interface ---------------------------//
  Scalar mClinear;
  Scalar mCquadratic;

  // Switch for the Balsara shear correction.
  bool mBalsaraShearCorrection;

  // Parameters for the Q limiter.
  Scalar mEpsilon2;
  Scalar mNegligibleSoundSpeed;
    
  // Maintain the last max viscous pressure for timestep control
  FieldList<Dimension, Scalar> mMaxViscousPressure;
  FieldList<Dimension, Scalar> mEffViscousPressure;

  // State for maintaining the velocity gradient
  bool mRigorousVelocityGradient;
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Tensor> mM;
  FieldList<Dimension, Tensor> mDvDx;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "ArtificialViscosityInline.hh"

#endif
