//---------------------------------Spheral++----------------------------------//
// CRKSPHBase -- Base class for the CRKSPH/ACRKSPH hydrodynamic packages
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHBase_hh__
#define __Spheral_CRKSPHBase_hh__

#include "Physics/GenericHydro.hh"
#include "Geometry/CellFaceFlag.hh"
#include "RK/RKCorrectionParams.hh"

#include <string>

namespace Spheral {
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosityHandle;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;
}

namespace Spheral {

template<typename Dimension>
class CRKSPHBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  CRKSPHBase(DataBase<Dimension>& dataBase,
             ArtificialViscosityHandle<Dimension>& Q,
             const RKOrder order,
             const double cfl,
             const bool useVelocityMagnitudeForDt,
             const bool compatibleEnergyEvolution,
             const bool evolveTotalEnergy,
             const bool XSPH,
             const MassDensityType densityUpdate,
             const double epsTensile,
             const double nTensile);

  // No default constructor, copying, or assignment.
  CRKSPHBase() = delete;
  CRKSPHBase(const CRKSPHBase&) = delete;
  CRKSPHBase& operator=(const CRKSPHBase&) = delete;

  // Destructor.
  virtual ~CRKSPHBase() = default;

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
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

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // We require RK corrections
  virtual std::set<RKOrder> requireReproducingKernels() const override;

  // The spatial order
  RKOrder correctionOrder()                                             const { return mOrder; }
  void correctionOrder(RKOrder x)                                             { mOrder = x; }

  // Flag to choose whether we want to sum for density, or integrate
  // the continuity equation.
  MassDensityType densityUpdate()                                       const { return mDensityUpdate; }
  void densityUpdate(MassDensityType x)                                       { mDensityUpdate = x; }

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution()                                      const { return mCompatibleEnergyEvolution; }
  void compatibleEnergyEvolution(bool x)                                      { mCompatibleEnergyEvolution = x; }

  // Flag controlling if we evolve total or specific energy.
  bool evolveTotalEnergy()                                              const { return mEvolveTotalEnergy; }
  void evolveTotalEnergy(bool x)                                              { mEvolveTotalEnergy = x; }

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH()                                                           const { return mXSPH; }
  void XSPH(bool x)                                                           { mXSPH = x; }

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile()                                               const { return mEpsTensile; }
  void epsilonTensile(Scalar x)                                               { mEpsTensile = x; }

  Scalar nTensile()                                                     const { return mnTensile; }
  void nTensile(Scalar x)                                                     { mnTensile = x; }
    
  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask()                 const { return mTimeStepMask; }
  const FieldList<Dimension, Scalar>&    pressure()                     const { return mPressure; }
  const FieldList<Dimension, Scalar>&    soundSpeed()                   const { return mSoundSpeed; }
  const FieldList<Dimension, Scalar>&    entropy()                      const { return mEntropy; }
  const FieldList<Dimension, Vector>&    XSPHDeltaV()                   const { return mXSPHDeltaV; }

  const FieldList<Dimension, Vector>&    DxDt()                         const { return mDxDt; }
  const FieldList<Dimension, Vector>&    DvDt()                         const { return mDvDt; }
  const FieldList<Dimension, Scalar>&    DmassDensityDt()               const { return mDmassDensityDt; }
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt()     const { return mDspecificThermalEnergyDt; }
  const FieldList<Dimension, Tensor>&    DvDx()                         const { return mDvDx; }
  const FieldList<Dimension, Tensor>&    internalDvDx()                 const { return mInternalDvDx; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label()                                  const override { return "CRKSPHBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  // A bunch of switches.
  RKOrder mOrder;
  MassDensityType mDensityUpdate;
  bool mCompatibleEnergyEvolution, mEvolveTotalEnergy, mXSPH;
  Scalar mEpsTensile, mnTensile;

  // Some internal scratch fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldList<Dimension, Scalar>    mEntropy;

  FieldList<Dimension, Vector>    mXSPHDeltaV;
  FieldList<Dimension, Vector>    mDxDt;

  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDHDt;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;

};

}

#endif
