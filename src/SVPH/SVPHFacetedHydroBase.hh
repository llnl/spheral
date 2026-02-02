//---------------------------------Spheral++----------------------------------//
// SVPHFacetedHydroBase -- The SVPHFaceted hydrodynamic package for Spheral++.
//
// Created by JMO, Sun Jul 28 20:57:01 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHFacetedHydroBase_hh__
#define __Spheral_SVPHFacetedHydroBase_hh__

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
class SVPHFacetedHydroBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors.
  SVPHFacetedHydroBase(const TableKernel<Dimension>& W,
                       ArtificialViscosityHandle<Dimension>& Q,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool XSVPH,
                       const bool linearConsistent,
                       const bool generateVoid,
                       const MassDensityType densityUpdate,
                       const Scalar fcentroidal,
                       const Scalar fcellPressure,
                       const Vector& xmin,
                       const Vector& xmax);

  // No default constructor, copying, or assignment.
  SVPHFacetedHydroBase() = delete;
  SVPHFacetedHydroBase(const SVPHFacetedHydroBase&) = delete;
  SVPHFacetedHydroBase& operator=(const SVPHFacetedHydroBase&) = delete;

  // Destructor.
  virtual ~SVPHFacetedHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

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
               
  // This algorithm does not use node->node connectivity.
  virtual bool requireConnectivity() const override { return false; }

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
  void densityUpdate(MassDensityType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(bool val);

  // Flag to determine if we're using the XSVPH algorithm.
  bool XSVPH() const;
  void XSVPH(bool val);

  // Flag to select whether or not to use the linear corrections.
  bool linearConsistent() const;
  void linearConsistent(bool val);

  // Flag to select whether or not to generate void points in the tessellation.
  bool generateVoid() const;
  void generateVoid(bool val);

  // Fraction of centroidal motion to apply each step.
  Scalar fcentroidal() const;
  void fcentroidal(Scalar val);

  // Fraction of the pressure to take from local cell.
  Scalar fcellPressure() const;
  void fcellPressure(Scalar val);

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
  // const FieldList<Dimension, std::vector<Scalar> >&    A() const;
  // const FieldList<Dimension, std::vector<Vector> >&    B() const;
  // const FieldList<Dimension, std::vector<Tensor> >&    gradB() const;
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    cellPressure() const;
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
  // const FieldList<Dimension, std::vector<Scalar> >& faceMass() const;
  // const FieldList<Dimension, std::vector<Vector> >& faceVelocity() const;
  // const FieldList<Dimension, std::vector<Vector> >& faceAcceleration() const;
  const FieldList<Dimension, std::vector<Vector> >& faceForce() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SVPHFacetedHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The interpolation kernel
  const TableKernel<Dimension>& mKernel;

  // A bunch of switches.
  MassDensityType mDensityUpdate;
  bool mCompatibleEnergyEvolution, mXSVPH, mLinearConsistent, mGenerateVoid;
  Scalar mfcentroidal, mfcellPressure;

  // Optional bounding box for generating the mesh.
  Vector mXmin, mXmax;

  // The mesh.
  using MeshPtr = std::shared_ptr<Mesh<Dimension> >;
  MeshPtr mMeshPtr;

  // Some internal scratch fields.
  // FieldList<Dimension, std::vector<Scalar> >&    mA;
  // FieldList<Dimension, std::vector<Vector> >&    mB;
  // FieldList<Dimension, std::vector<Tensor> >&    mGradB;
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mCellPressure;
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

  // FieldList<Dimension, std::vector<Scalar> >    mFaceMass;
  // FieldList<Dimension, std::vector<Vector> >    mFaceVelocity;
  // FieldList<Dimension, std::vector<Vector> >    mFaceAcceleration;
  FieldList<Dimension, std::vector<Vector> >    mFaceForce;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "SVPHFacetedHydroBaseInline.hh"

#endif
