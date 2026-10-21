//---------------------------------Spheral++----------------------------------//
// TensorDamageModel -- Base class for the tensor damage physics models.
// This class does not know how to seed the flaw distribution -- that is 
// required of descendant classes.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Thu Sep 29 15:42:05 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorDamageModel_hh__
#define __Spheral_TensorDamageModel_hh__

#include "DamageModel.hh"

#include <vector>

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  template<typename Dimension> class SolidNodeList;
  template<typename Dimension> class DataBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
  class FileIO;
}

namespace Spheral {

template<typename Dimension>
class TensorDamageModel: 
    public DamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;
  using FlawStorageType = Field<Dimension, std::vector<double> >;

  // Constructors, destructor.
  TensorDamageModel(SolidNodeList<Dimension>& nodeList,
                    const TensorStrainAlgorithm strainAlgorithm,
                    const DamageCouplingAlgorithm damageCouplingAlgorithm,
                    const TableKernel<Dimension>& W,
                    const double crackGrowthMultiplier,
                    const double criticalDamageThreshold,
                    const bool damageInCompression,
                    const FlawStorageType& flaws);
  TensorDamageModel(SolidNodeList<Dimension>& nodeList,
                    const TensorStrainAlgorithm strainAlgorithm,
                    const DamageCouplingAlgorithm damageCouplingAlgorithm,
                    const TableKernel<Dimension>& W,
                    const double crackGrowthMultiplier,
                    const double criticalDamageThreshold,
                    const bool damageInCompression);
  virtual ~TensorDamageModel();

  //...........................................................................
  // Provide the required physics package interface.
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

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  //...........................................................................
  // Optional method to cull the set of flaws to the single weakest one on
  // each point.
  void cullToWeakestFlaws();

  // Get the set of flaw activation energies for the given node index.
  const std::vector<double> flawsForNode(const size_t index) const;

  // Compute a Field with the sum of the activation energies per node.
  Field<Dimension, Scalar> sumActivationEnergiesPerNode() const;

  // Compute a Field with the number of flaws per node.
  Field<Dimension, Scalar> numFlawsPerNode() const;

  // Provide access to the state fields we maintain.
  const Field<Dimension, Scalar>& youngsModulus() const;
  const Field<Dimension, Scalar>& longitudinalSoundSpeed() const;
  const Field<Dimension, SymTensor>& strain() const;
  const Field<Dimension, SymTensor>& effectiveStrain() const;
  const Field<Dimension, Scalar>& DdamageDt() const;
  const FlawStorageType& flaws() const;
  FlawStorageType& flaws();

  void flaws(const FlawStorageType& x);

  // The algorithms to update the strain.
  TensorStrainAlgorithm strainAlgorithm() const;

  // Flag to determine if damage in compression is allowed.
  bool damageInCompression() const;
  void damageInCompression(bool x);

  // The critical damage threshold for not setting the time step.
  double criticalDamageThreshold() const; 
  void criticalDamageThreshold(double x);

  //**************************************************************************
  // Restart methods.
  virtual std::string label() const override { return "TensorDamageModel"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //**************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  FlawStorageType mFlaws;
  Field<Dimension, Scalar> mYoungsModulus;
  Field<Dimension, Scalar> mLongitudinalSoundSpeed;
  Field<Dimension, SymTensor> mStrain;
  Field<Dimension, SymTensor> mEffectiveStrain;
  Field<Dimension, Scalar> mDdamageDt;

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStrainAlgorithm mStrainAlgorithm;
  double mCriticalDamageThreshold, mCriticalNodesPerSmoothingScale;
  bool mDamageInCompression;

  // No default constructor, copying or assignment.
  TensorDamageModel();
  TensorDamageModel(const TensorDamageModel&);
  TensorDamageModel& operator=(const TensorDamageModel&);
};

}

#include "TensorDamageModelInline.hh"

#endif

