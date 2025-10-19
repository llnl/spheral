//---------------------------------Spheral++----------------------------------//
// MFV -- This is an Arbitrary Eulerian-Lagrangian extension of the
//        MFV approach of Hopkins 2015. Its got several node-motion
//        approaches which promote more regular particle distributions.
//
//        Each of the ALE options defines the velocity of the nodes 
//        differently. The flux that results from the difference
//        between the nodes velocities and the fluid velocity.
//        The velocities are defined as follows for the 
//        NodeMotionTypes:
//
//        1) Eulerian ---- static Nodes
//        2) Lagrangian -- nodal velocity = fluid velocity. (This is
//                         a spheralized version of MFV so there
//                         is some flux between nodes)
//        3) Fickian ------ nodal velocity = fluid velocity + Fician
//                         PST correction
//        4) XSPH -------- nodal velocity = xsph velocity
//
//   Hopkins P.F. (2015) "A New Class of Accurate, Mesh-Free Hydrodynamic 
//   Simulation Methods," MNRAS, 450(1):53-110
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
// TODO:
//   1 backpressure and fician particle shifting
//   2 Eulerian model will still crash on the Noh implosion due to void particles
//   3 Good implementation of Ngb update
//   4 treatment for material interfaces
//---------------------------------------------------------------------------//

#ifndef __Spheral_MFV_hh__
#define __Spheral_MFV_hh__

#include <string>
#include <memory>

#include "GSPH/GenericRiemannHydro.hh"

namespace Spheral {

enum class NodeMotionType {
  Lagrangian = 0,
  Eulerian = 1,
  Fician = 2,
  XSPH = 3,
  BackgroundPressure = 4,
};

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class SmoothingScaleBase;
template<typename Dimension> class TableKernel;
template<typename Dimension> class RiemannSolverBase;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension, typename Value, size_t numElements> class PairwiseField;
class FileIO;

template<typename Dimension>
class MFV: public GenericRiemannHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  typedef typename GenericRiemannHydro<Dimension>::TimeStepType TimeStepType;
  typedef typename GenericRiemannHydro<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  using PairAccelerationsType = typename GenericRiemannHydro<Dimension>::PairAccelerationsType;
  using PairWorkType = typename GenericRiemannHydro<Dimension>::PairWorkType;
  using PairMassFluxType = PairwiseField<Dimension, Scalar, 1u>;

  // Constructors.
  MFV(DataBase<Dimension>& dataBase,
               RiemannSolverBase<Dimension>& riemannSolver,
               const TableKernel<Dimension>& W,
               const Scalar epsDiffusionCoeff,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool compatibleEnergyEvolution,
               const bool evolveTotalEnergy,
               const bool XSPH,
               const bool correctVelocityGradient,
               const double nodeMotionCoefficient,
               const NodeMotionType nodeMotionType,
               const GradientType gradType,
               const MassDensityType densityUpdate,
               const double epsTensile,
               const double nTensile,
               const Vector& xmin,
               const Vector& xmax);

  // No default constructor, copying, or assignment.
  MFV() = delete;
  MFV(const MFV&) = delete;
  MFV& operator=(const MFV&) = delete;

  // Destructor.
  virtual ~MFV() = default;

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

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

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  bool initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) override;
                       
  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void
  firstDerivativesLoop(const typename Dimension::Scalar time,
                       const typename Dimension::Scalar dt,
                       const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state,
                             StateDerivatives<Dimension>& derivatives) const override;
  virtual
  void
  secondDerivativesLoop(const typename Dimension::Scalar time,
                        const typename Dimension::Scalar dt,
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

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  Scalar nodeMotionCoefficient() const;
  void nodeMotionCoefficient(const Scalar x);

  NodeMotionType nodeMotionType() const;
  void nodeMotionType(NodeMotionType x);

  const FieldList<Dimension,Scalar>& DmassDt() const;
  const FieldList<Dimension,Scalar>& DthermalEnergyDt() const;
  const FieldList<Dimension,Vector>& DmomentumDt() const;
  const FieldList<Dimension,Scalar>& DvolumeDt() const;

  const PairMassFluxType& pairMassFlux() const;
  
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "MFV" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************           
private:

  Scalar mNodeMotionCoefficient; 

  NodeMotionType mNodeMotionType;

  FieldList<Dimension, Scalar> mDmassDt;
  FieldList<Dimension, Scalar> mDthermalEnergyDt;
  FieldList<Dimension, Vector> mDmomentumDt;
  FieldList<Dimension, Scalar> mDvolumeDt;

  std::unique_ptr<PairMassFluxType> mPairMassFluxPtr;
};

}

#include "MFVInline.hh"

#endif
