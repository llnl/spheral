//---------------------------------Spheral++----------------------------------//
// GenericRiemannHydro --  pure virtual class for hydros using a Riemann
//                         solver
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_GenericRiemannHydro_hh__
#define __Spheral_GenericRiemannHydro_hh__

#include <string>
#include <memory>

#include "Physics/GenericHydro.hh"
#include "Physics/Physics.hh"
#include "Neighbor/PairwiseField.hh"

namespace Spheral {

enum class GradientType {
  RiemannGradient = 0,
  HydroAccelerationGradient = 1,
  SPHGradient = 2,
  MixedMethodGradient = 3,
  SPHSameTimeGradient = 4,
  SPHUncorrectedGradient = 5,
  NoGradient = 6
};

enum class GSPHEvolutionType {
  IdealH = 0,
  IntegrateH = 1,
  constantNeighborCount = 2
};

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension> class RiemannSolverBase;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
class FileIO;

template<typename Dimension>
class GenericRiemannHydro: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using PairWorkType = PairwiseField<Dimension, Scalar, 2u>;

  // Constructors.
  GenericRiemannHydro(DataBase<Dimension>& dataBase,
                      RiemannSolverBase<Dimension>& riemannSolver,
                      const TableKernel<Dimension>& W,
                      const Scalar epsDiffusionCoeff,
                      const double cfl,
                      const bool useVelocityMagnitudeForDt,
                      const bool compatibleEnergyEvolution,
                      const bool evolveTotalEnergy,
                      const bool XSPH,
                      const bool correctVelocityGradient,
                      const GradientType gradType,
                      const MassDensityType densityUpdate,
                      const double epsTensile,
                      const double nTensile,
                      const Vector& xmin,
                      const Vector& xmax);

  // No default constructor, copying, or assignment.
  GenericRiemannHydro() = delete;
  GenericRiemannHydro(const GenericRiemannHydro&) = delete;
  GenericRiemannHydro& operator=(const GenericRiemannHydro&) = delete;

  // Destructor.
  virtual ~GenericRiemannHydro() = default;

    // We require all Physics packages to provide a method returning their vote
  // for the next time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

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

  
  // This method is called once at the beginning of a timestep, after all state registration.
  // virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
  //                                State<Dimension>& state,
  //                                StateDerivatives<Dimension>& derivs) override;

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  bool initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
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
                  

  // Access the stored riemann solver
  RiemannSolverBase<Dimension>& riemannSolver() const;

  // Access the stored interpolation kernels.
  const TableKernel<Dimension>& kernel() const;

  GradientType gradientType() const;
  void gradientType(GradientType x);

  // Flag for our density update
  MassDensityType densityUpdate() const;
  void densityUpdate(MassDensityType type);

  // setter-getters for our bool switches
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(bool val);

  bool evolveTotalEnergy() const;
  void evolveTotalEnergy(bool val);

  bool useVelocityMagnitudeForDt() const;
  void useVelocityMagnitudeForDt(bool x);

  bool XSPH() const;
  void XSPH(bool val);
  
  bool correctVelocityGradient() const;
  void correctVelocityGradient(bool val);

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile() const;
  void epsilonTensile(Scalar val);

  Scalar nTensile() const;
  void nTensile(Scalar val);

  // diffusion coefficient for specific thermal energy
  Scalar specificThermalEnergyDiffusionCoefficient() const;
  void specificThermalEnergyDiffusionCoefficient(const Scalar x);

  // Also allow access to the CFL timestep safety criteria.
  Scalar cfl() const;
  void cfl(Scalar cfl);

  // Optionally we can provide a bounding box
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    volume() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    normalization() const;
  const FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldList<Dimension, Tensor>&    M() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  
  const PairAccelerationsType&           pairAccelerations() const;
  const PairWorkType&                    pairDepsDt() const;

  const FieldList<Dimension, Vector>&    DrhoDx() const;
  const FieldList<Dimension, Vector>&    riemannDpDx() const;
  const FieldList<Dimension, Tensor>&    riemannDvDx() const;
  const FieldList<Dimension, Vector>&    newRiemannDpDx() const;
  const FieldList<Dimension, Tensor>&    newRiemannDvDx() const;
  
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "GenericRiemannHydro" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  RestartRegistrationType mRestart;

  
private:
  //--------------------------- Private Interface ---------------------------//
  RiemannSolverBase<Dimension>& mRiemannSolver;
  const TableKernel<Dimension>& mKernel;
  GradientType mGradientType;
  MassDensityType mDensityUpdate;
  
  // A bunch of switches.
  bool mCompatibleEnergyEvolution;    
  bool mEvolveTotalEnergy;           
  bool mXSPH;
  bool mCorrectVelocityGradient;
  bool mUseVelocityMagnitudeForDt;
  
  Scalar mEpsTensile, mnTensile;                      
  Scalar mSpecificThermalEnergyDiffusionCoefficient; 
  Scalar mCfl; 
  Vector mxmin, mxmax;

  // Our fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mVolume;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;

  FieldList<Dimension, Scalar>    mNormalization;

  FieldList<Dimension, Scalar>    mXSPHWeightSum;
  FieldList<Dimension, Vector>    mXSPHDeltaV;

  FieldList<Dimension, Tensor>    mM;

  FieldList<Dimension, Vector>    mDxDt;
  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDHDt;

  FieldList<Dimension, Vector>    mDrhoDx;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Vector>    mRiemannDpDx;
  FieldList<Dimension, Tensor>    mRiemannDvDx;
  FieldList<Dimension, Vector>    mNewRiemannDpDx;
  FieldList<Dimension, Tensor>    mNewRiemannDvDx;

  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;    // store pairwise contribution to DvDt for compatible
  std::unique_ptr<PairWorkType>          mPairDepsDtPtr;           // store pairwise contribution to DepsDt for compatible
};

}

#include "GenericRiemannHydroInline.hh"

#endif
