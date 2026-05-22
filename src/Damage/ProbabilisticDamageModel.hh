//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamageModel
//
// A damage model based on Weibull statistics that uses volume based
// probabilities per node to decide when damage starts to accrue.  Should
// generate similar results to the classic Benz-Asphaug (Grady-Kipp) model
// without generating explicit flaws.  Also appropriate for use with varying
// resolution materials.
//
// Created by JMO, Tue Apr 13 15:58:08 PDT 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_ProbabilisticDamageModel_hh__
#define __Spheral_ProbabilisticDamageModel_hh__

#include "DamageModel.hh"
#include <memory>

namespace Spheral {

template<typename Dimension>
class ProbabilisticDamageModel: public DamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors, destructor.
  ProbabilisticDamageModel(SolidNodeList<Dimension>& nodeList,
                           const TableKernel<Dimension>& W,
                           const double kWeibull,
                           const double mWeibull,
                           const size_t seed,
                           const size_t minFlawsPerNode,
                           const double crackGrowthMultiplier,
                           const double volumeMultiplier,
                           const DamageCouplingAlgorithm damageCouplingAlgorithm,
                           const TensorStrainAlgorithm strainAlgorithm,
                           const bool damageInCompression,
                           const double criticalDamageThreshold,
                           const Field<Dimension, int>& mask);
  virtual ~ProbabilisticDamageModel() = default;

  //............................................................................
  // Override the Physics package interface.

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  // Compute the derivatives.
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  //............................................................................
  // Accessors for state
  TensorStrainAlgorithm strainAlgorithm() const;
  bool damageInCompression() const;
  double kWeibull() const;
  double mWeibull() const;
  double volumeMultiplier() const;
  double Vmin() const;
  double Vmax() const;
  size_t seed() const;
  size_t minFlawsPerNode() const;
  const Field<Dimension, int>& numFlaws() const;
  const Field<Dimension, Scalar>& minFlaw() const;
  const Field<Dimension, Scalar>& maxFlaw() const;
  const Field<Dimension, Scalar>& initialVolume() const;
  const Field<Dimension, Scalar>& youngsModulus() const;
  const Field<Dimension, Scalar>& longitudinalSoundSpeed() const;
  const Field<Dimension, SymTensor>& strain() const;
  const Field<Dimension, SymTensor>& effectiveStrain() const;
  const Field<Dimension, Scalar>& DdamageDt() const;
  const Field<Dimension, Scalar>& damageTT() const;
  const Field<Dimension, Scalar>& strainTT() const;
  const Field<Dimension, Scalar>& effectiveStrainTT() const;

  // Optionally the user can provide a mask to prevent damage modeling on some points.
  const Field<Dimension, int>& mask() const;
  void mask(const Field<Dimension, int>& val);

  // Optionally ignore timestep votes for material beyond a damage threshold
  double criticalDamageThreshold() const;
  void criticalDamageThreshold(const double val);

  // No default constructor, copying or assignment.
  ProbabilisticDamageModel() = delete;
  ProbabilisticDamageModel(const ProbabilisticDamageModel&) = delete;
  ProbabilisticDamageModel& operator=(const ProbabilisticDamageModel&) = delete;

  //............................................................................
  // Restart methods.
  virtual std::string label() const override { return "ProbabilisticDamageModel"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStrainAlgorithm mStrainAlgorithm;
  bool mDamageInCompression;
  double mkWeibull, mmWeibull, mVolumeMultiplier, mVmin, mVmax, mCriticalDamageThreshold;
  size_t mSeed, mMinFlawsPerNode;
  Field<Dimension, int> mNumFlaws, mMask;
  Field<Dimension, Scalar> mMinFlaw, mMaxFlaw, mInitialVolume, mYoungsModulus, mLongitudinalSoundSpeed, mDdamageDt;
  Field<Dimension, SymTensor> mStrain, mEffectiveStrain;

  // In RZ we need an extra tensor component
  std::unique_ptr<Field<Dimension, Scalar>> mDamageTT, mStrainTT, mEffectiveStrainTT;
};

}

#include "ProbabilisticDamageModelInline.hh"

#endif

