//---------------------------------Spheral++----------------------------------//
// DamageModel -- Base class for the damage physics models.
// This class just provides the basic interface for damage models, and does 
// not fill out the complete physics package interface.
//
// Created by JMO, Thu Sep 29 13:31:57 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamageModel_hh__
#define __Spheral_DamageModel_hh__

#include "Physics/Physics.hh"
#include "Utilities/NodeCoupling.hh"
#include "DataOutput/registerWithRestart.hh"

#include <vector>
#include <memory>

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  template<typename Dimension> class SolidNodeList;
  template<typename Dimension> class DataBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
  template<typename Dimension> class TableKernel;
  class FileIO;
}

namespace Spheral {

enum class EffectiveFlawAlgorithm {
  FullSpectrumFlaws = 0,
  MinFlaw = 1,
  MaxFlaw = 2,
  InverseSumFlaws = 3,
  SampledFlaws = 4,
};

enum class DamageCouplingAlgorithm {
  DirectDamage = 0,
  PairMaxDamage = 1,
  DamageGradient = 2,
  ThreePointDamage = 3,
  TensorPairMaxDamage = 4,
};

template<typename Dimension>
class DamageModel: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;
  using NodeCouplingPtr = typename std::shared_ptr<NodeCoupling>;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;
  using FlawStorageType = Field<Dimension, std::vector<double>>;

  // Constructors, destructor.
  DamageModel(SolidNodeList<Dimension>& nodeList,
              const TableKernel<Dimension>& W,
              const double crackGrowthMultiplier,
              const DamageCouplingAlgorithm damageCouplingAlgorithm);
  virtual ~DamageModel() = default;

  // Compute the generic Grady-Kipp (ala Benz-Asphaug) scalar damage time 
  // derivative.
  virtual 
  void computeScalarDDDt(const DataBase<Dimension>& dataBase,
                         const State<Dimension>& state,
                         const Scalar time,
                         const Scalar dt,
                         Field<Dimension, Scalar>& DDDt) const;

  //...........................................................................
  // Provide a subset of the required physics package interface.
  virtual bool initialize(const Scalar time,
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase,
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivatives) override;

  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  virtual bool requireGhostConnectivity()               const override { return mDamageCouplingAlgorithm == DamageCouplingAlgorithm::ThreePointDamage; }
  virtual bool requireIntersectionConnectivity()        const override { return mComputeIntersectConnectivity; }

  // Return the maximum state change we care about for checking for convergence in the implicit integration methods.
  // Damage model for now defaults to just relying on the hydro to handle this, so return no vote.
  virtual ResidualType maxResidual(const DataBase<Dimension>& dataBase, 
                                   const State<Dimension>& state1,
                                   const State<Dimension>& state0,
                                   const Scalar tol) const override;
  //...........................................................................
  // Access the SolidNodeList we're damaging.
  SolidNodeList<Dimension>& nodeList()                                 { return mNodeList; }
  const SolidNodeList<Dimension>& nodeList()                     const { return mNodeList; }

  // Access the kernel.
  const TableKernel<Dimension>& kernel()                         const { return mW; }

  // Important local parameters.
  double crackGrowthMultiplier()                                 const { return mCrackGrowthMultiplier; }
  DamageCouplingAlgorithm damageCouplingAlgorithm()              const { return mDamageCouplingAlgorithm; }
  const NodeCoupling& nodeCoupling()                             const { return *mNodeCouplingPtr; }

  // Allow the user to specify a set of nodes to be excluded from damage.
  std::vector<int> excludeNodes() const;
  void excludeNodes(std::vector<int> ids);

  // Optionally freeze the damage (no more damage growth)
  bool freezeDamage()                                            const { return mFreezeDamage; }
  void freezeDamage(const bool x)                                      { mFreezeDamage = x; }

  //**************************************************************************
  // Restart methods.
  virtual std::string label()                           const override { return "DamageModel"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //**************************************************************************

  // No default constructor, copying or assignment.
  DamageModel() = delete;
  DamageModel(const DamageModel&) = delete;
  DamageModel& operator=(const DamageModel&) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  SolidNodeList<Dimension>& mNodeList;
  const TableKernel<Dimension>& mW;
  double mCrackGrowthMultiplier;
  DamageCouplingAlgorithm mDamageCouplingAlgorithm;
  Field<Dimension, int> mExcludeNode;
  NodeCouplingPtr mNodeCouplingPtr;
  bool mComputeIntersectConnectivity, mFreezeDamage;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#endif

