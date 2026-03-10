//---------------------------------Spheral++----------------------------------//
// SPHBase -- Base class for SPH/ASPH hydrodynamic packages for Spheral++.
// 
// This class contains most of the boilerplate storage for SPH specliazations
// to inherit and use, with the exception of pair-wise accelerations for
// the compatible energy discretization.
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHBase__
#define __Spheral_SPHBase__

#include "Physics/GenericHydro.hh"

#include <string>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosityHandle;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
class FileIO;

template<typename Dimension>
class SPHBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SPHBase(DataBase<Dimension>& dataBase,
          ArtificialViscosityHandle<Dimension>& Q,
          const TableKernel<Dimension>& W,
          const TableKernel<Dimension>& WPi,
          const double cfl,
          const bool useVelocityMagnitudeForDt,
          const bool compatibleEnergyEvolution,
          const bool evolveTotalEnergy,
          const bool gradhCorrection,
          const bool XSPH,
          const bool correctVelocityGradient,
          const bool sumMassDensityOverAllNodeLists,
          const MassDensityType densityUpdate,
          const double epsTensile,
          const double nTensile,
          const Vector& xmin,
          const Vector& xmax);

  // No default constructor, copying, or assignment.
  SPHBase() = delete;
  SPHBase(const SPHBase&) = delete;
  SPHBase& operator=(const SPHBase&) = delete;

  // Destructor.
  virtual ~SPHBase() = default;

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual
  void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                            State<Dimension>& state,
                                            StateDerivatives<Dimension>& derivs) override;

  // Register the state
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Stuff to do post-state updates
  virtual 
  bool postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // A method to fill in the volume in the State, optionally enforcing
  // boundary conditions.
  void updateVolume(State<Dimension>& state,
                    const bool boundaries) const;

  // Flag to choose whether we want to sum for density, or integrate the continuity equation.
  MassDensityType densityUpdate()                                               const { return mDensityUpdate; }
  void densityUpdate(MassDensityType x)                                               { mDensityUpdate = x; }

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution()                                              const { return mCompatibleEnergyEvolution; }
  void compatibleEnergyEvolution(bool x)                                              { mCompatibleEnergyEvolution = x; }

  // Flag controlling if we evolve total or specific energy.
  bool evolveTotalEnergy()                                                      const { return mEvolveTotalEnergy; }
  void evolveTotalEnergy(bool x)                                                      { mEvolveTotalEnergy = x; }

  // Flag to determine if we're using the grad h correction.
  bool gradhCorrection()                                                        const { return mGradhCorrection; }
  void gradhCorrection(bool x)                                                        { mGradhCorrection = x; }

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH()                                                                   const { return mXSPH; }
  void XSPH(bool x)                                                                   { mXSPH = x; }

  // Flag to determine if we're applying the linear correction for the velocity gradient.
  bool correctVelocityGradient()                                                const { return mCorrectVelocityGradient; }
  void correctVelocityGradient(bool x)                                                { mCorrectVelocityGradient = x; }

  // Flag to determine if the sum density definition extends over neighbor NodeLists.
  bool sumMassDensityOverAllNodeLists()                                         const { return mSumMassDensityOverAllNodeLists; }
  void sumMassDensityOverAllNodeLists(bool x)                                         { mSumMassDensityOverAllNodeLists = x; }

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile()                                                       const { return mEpsTensile; }
  void epsilonTensile(Scalar x)                                                       { mEpsTensile = x; }

  Scalar nTensile()                                                             const { return mnTensile; }
  void nTensile(Scalar x)                                                             { mnTensile = x; }

  // Optionally we can provide a bounding box for use generating the mesh
  // for the Voronoi mass density update.
  const Vector& xmin()                                                          const { return mxmin; }
  const Vector& xmax()                                                          const { return mxmax; }
  void xmin(const Vector& x)                                                          { mxmin = x; }
  void xmax(const Vector& x)                                                          { mxmax = x; }

  // Access the stored interpolation kernels.
  const TableKernel<Dimension>& kernel()                                        const { return mKernel; }
  const TableKernel<Dimension>& PiKernel()                                      const { return mPiKernel; }

  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&          timeStepMask()                      const { return mTimeStepMask; }
  const FieldList<Dimension, Scalar>&       pressure()                          const { return mPressure; }
  const FieldList<Dimension, Scalar>&       soundSpeed()                        const { return mSoundSpeed; }
  const FieldList<Dimension, Scalar>&       volume()                            const { return mVolume; }
  const FieldList<Dimension, Scalar>&       omegaGradh()                        const { return mOmegaGradh; }
  const FieldList<Dimension, Scalar>&       entropy()                           const { return mEntropy; }
  const FieldList<Dimension, Scalar>&       massDensityCorrection()             const { return mMassDensityCorrection; }
  const FieldList<Dimension, Scalar>&       massDensitySum()                    const { return mMassDensitySum; }
  const FieldList<Dimension, Scalar>&       normalization()                     const { return mNormalization; }
  const FieldList<Dimension, Scalar>&       XSPHWeightSum()                     const { return mXSPHWeightSum; }
  const FieldList<Dimension, Vector>&       XSPHDeltaV()                        const { return mXSPHDeltaV; }
  const FieldList<Dimension, Tensor>&       M()                                 const { return mM; }
  const FieldList<Dimension, Tensor>&       localM()                            const { return mLocalM; }
  const FieldList<Dimension, Vector>&       DxDt()                              const { return mDxDt; }
  const FieldList<Dimension, Vector>&       DvDt()                              const { return mDvDt; }
  const FieldList<Dimension, Scalar>&       DmassDensityDt()                    const { return mDmassDensityDt; }
  const FieldList<Dimension, Scalar>&       DspecificThermalEnergyDt()          const { return mDspecificThermalEnergyDt; }
  const FieldList<Dimension, Tensor>&       DvDx()                              const { return mDvDx; }
  const FieldList<Dimension, Tensor>&       internalDvDx()                      const { return mInternalDvDx; }
  const FieldList<Dimension, Vector>&       gradRho()                           const { return mGradRho; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SPHBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The interpolation kernels.
  const TableKernel<Dimension>& mKernel;
  const TableKernel<Dimension>& mPiKernel;

  // A bunch of switches.
  MassDensityType mDensityUpdate;
  bool mCompatibleEnergyEvolution, mEvolveTotalEnergy, mGradhCorrection, mXSPH, mCorrectVelocityGradient, mSumMassDensityOverAllNodeLists;

  // Tensile correction.
  Scalar mEpsTensile, mnTensile;

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // Some internal scratch fields.
  FieldList<Dimension, int>                mTimeStepMask;
  FieldList<Dimension, Scalar>             mPressure;
  FieldList<Dimension, Scalar>             mSoundSpeed;
  FieldList<Dimension, Scalar>             mOmegaGradh;
  FieldList<Dimension, Scalar>             mEntropy;

  FieldList<Dimension, Scalar>             mMassDensityCorrection;
  FieldList<Dimension, Scalar>             mMassDensitySum;
  FieldList<Dimension, Scalar>             mNormalization;

  FieldList<Dimension, Scalar>             mWeightedNeighborSum;
  FieldList<Dimension, Vector>             mMassFirstMoment;
  FieldList<Dimension, SymTensor>          mMassSecondMomentEta;
  FieldList<Dimension, SymTensor>          mMassSecondMomentLab;

  FieldList<Dimension, Scalar>             mXSPHWeightSum;
  FieldList<Dimension, Vector>             mXSPHDeltaV;

  FieldList<Dimension, Vector>             mDxDt;
  FieldList<Dimension, Vector>             mDvDt;
  FieldList<Dimension, Scalar>             mDmassDensityDt;
  FieldList<Dimension, Scalar>             mDspecificThermalEnergyDt;
  FieldList<Dimension, Tensor>             mDvDx;
  FieldList<Dimension, Tensor>             mInternalDvDx;
  FieldList<Dimension, Vector>             mGradRho;
  FieldList<Dimension, Tensor>             mM;
  FieldList<Dimension, Tensor>             mLocalM;

  FieldList<Dimension, Scalar>             mVolume;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#endif
