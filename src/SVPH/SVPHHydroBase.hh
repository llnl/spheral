//---------------------------------Spheral++----------------------------------//
// SVPHHydroBase -- The SVPH hydrodynamic package for Spheral++.
//
// Created by JMO, Sun Jul 28 20:57:01 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHHydroBase_hh__
#define __Spheral_SVPHHydroBase_hh__

#include <string>

#include "Physics/GenericHydro.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosityHandle;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class SVPHHydroBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SVPHHydroBase(const TableKernel<Dimension>& W,
                ArtificialViscosityHandle<Dimension>& Q,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool XSVPH,
                const bool linearConsistent,
                const MassDensityType densityUpdate,
                const Scalar fcentroidal,
                const Vector& xmin,
                const Vector& xmax);

  // No default constructor, copying, or assignment.
  SVPHHydroBase() = delete;
  SVPHHydroBase(const SVPHHydroBase&) = delete;
  SVPHHydroBase& operator=(const SVPHHydroBase&) = delete;

  // Destructor.
  virtual ~SVPHHydroBase();

  // Tasks we do once on problem startup.
  // An optional hook to initialize once when the problem is starting up.
  // Typically this is used to size arrays once all the materials and NodeLists have
  // been created.  It is assumed after this method has been called it is safe to
  // call Physics::registerState for instance to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
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
                               StateDerivatives<Dimension>& derivatives) const;

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Finalize the hydro at the completion of an integration step.
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) override;
               
  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Flag to choose whether we want to sum for density, or integrate
  // the continuity equation.
  MassDensityType densityUpdate() const;
  void densityUpdate(const MassDensityType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(const bool val);

  // Flag to determine if we're using the XSVPH algorithm.
  bool XSVPH() const;
  void XSVPH(const bool val);

  // Flag to select whether or not to use the linear corrections.
  bool linearConsistent() const;
  void linearConsistent(const bool val);

  // Fraction of centroidal motion to apply each step.
  Scalar fcentroidal() const;
  void fcentroidal(const Scalar val);

  // Optionally we can provide a bounding box for use generating the mesh.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  // Access the stored interpolation kernel
  const TableKernel<Dimension>& kernel() const;

  // The tessellation.
  const Mesh<Dimension>& mesh() const;

  // The state field lists we're maintaining.
  const FieldList<Dimension, Scalar>&    A() const;
  const FieldList<Dimension, Vector>&    B() const;
  const FieldList<Dimension, Tensor>&    gradB() const;
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    volume() const;
  const FieldList<Dimension, Scalar>&    massDensitySum() const;
  const FieldList<Dimension, Vector>&    XSVPHDeltaV() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldList<Dimension, Tensor>&    internalDvDx() const;
  const FieldList<Dimension, std::vector<Vector> >& pairAccelerations() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SVPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The interpolation kernel
  const TableKernel<Dimension>& mKernel;

  // A bunch of switches.
  MassDensityType mDensityUpdate;
  bool mCompatibleEnergyEvolution, mXSVPH, mLinearConsistent;
  Scalar mfcentroidal;

  // Optional bounding box for generating the mesh.
  Vector mXmin, mXmax;

  // The mesh.
  using MeshPtr = std::shared_ptr<Mesh<Dimension> >;
  MeshPtr mMeshPtr;

  // Some internal scratch fields.
  FieldList<Dimension, Scalar>    mA;
  FieldList<Dimension, Vector>    mB;
  FieldList<Dimension, Tensor>    mGradB;
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;

  FieldList<Dimension, Scalar>    mMassDensitySum;

  FieldList<Dimension, Vector>    mXSVPHDeltaV;

  FieldList<Dimension, Vector>    mDxDt;
  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;

  FieldList<Dimension, Scalar>    mVolume;

  FieldList<Dimension, std::vector<Vector> > mPairAccelerations;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "SVPHHydroBaseInline.hh"

#endif
