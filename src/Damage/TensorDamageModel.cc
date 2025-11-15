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
#include "FileIO/FileIO.hh"
#include "TensorDamageModel.hh"
#include "TensorStrainPolicy.hh"
#include "TensorDamagePolicy.hh"
#include "YoungsModulusPolicy.hh"
#include "LongitudinalSoundSpeedPolicy.hh"
#include "DamageGradientPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/updateStateFields.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"

#include <string>
#include <vector>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamageModel<Dimension>::
TensorDamageModel(SolidNodeList<Dimension>& nodeList,
                  const TensorStrainAlgorithm strainAlgorithm,
                  const DamageCouplingAlgorithm damageCouplingAlgorithm,
                  const TableKernel<Dimension>& W,
                  const double crackGrowthMultiplier,
                  const double criticalDamageThreshold,
                  const bool damageInCompression,
                  const FlawStorageType& flaws):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mFlaws(SolidFieldNames::flaws, flaws),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList),
  mDdamageDt(TensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrainAlgorithm(strainAlgorithm),
  mCriticalDamageThreshold(criticalDamageThreshold),
  mDamageInCompression(damageInCompression) {
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamageModel<Dimension>::
TensorDamageModel(SolidNodeList<Dimension>& nodeList,
                  const TensorStrainAlgorithm strainAlgorithm,
                  const DamageCouplingAlgorithm damageCouplingAlgorithm,
                  const TableKernel<Dimension>& W,
                  const double crackGrowthMultiplier,
                  const double criticalDamageThreshold,
                  const bool damageInCompression):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mFlaws(SolidFieldNames::flaws, nodeList),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList),
  mDdamageDt(TensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrainAlgorithm(strainAlgorithm),
  mCriticalDamageThreshold(criticalDamageThreshold),
  mDamageInCompression(damageInCompression) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamageModel<Dimension>::
~TensorDamageModel() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the tensor damage field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Set the scalar magnitude of the damage evolution.
  const auto* nodeListPtr = &(this->nodeList());
  auto&       DDDt = derivs.field(state.buildFieldKey(TensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeListPtr->name()), 0.0);
  this->computeScalarDDDt(dataBase,
                          state,
                          time,
                          dt,
                          DDDt);
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename TensorDamageModel<Dimension>::TimeStepType
TensorDamageModel<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {

  // // Look at how quickly we're trying to change the damage.
  // double dt = DBL_MAX;
  // const Field<Dimension, SymTensor>& damage = this->nodeList().damage();
  // const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  // const size_t nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &(this->nodeList())));
  // for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
  //      iItr != connectivityMap.end(nodeListi);
  //      ++iItr) {
  //   const int i = *iItr;
  //   const double D0 = damage(i).Trace() / Dimension::nDim;
  //   dt = min(dt, 0.8*max(D0, 1.0 - D0)/
  //            std::sqrt(mDdamageDt(i)*mDdamageDt(i) + 1.0e-20));
  // }
  // return TimeStepType(dt, "Rate of damage change");

  return TimeStepType(1.0e100, "Rate of damage change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Register Youngs modulus and the longitudinal sound speed.
  auto& nodes = this->nodeList();
  state.enroll(mYoungsModulus, std::make_shared<YoungsModulusPolicy<Dimension>>(nodes));
  state.enroll(mLongitudinalSoundSpeed, std::make_shared<LongitudinalSoundSpeedPolicy<Dimension>>(nodes));

  // Register the strain and effective strain.
  state.enroll(mStrain);
  state.enroll(mEffectiveStrain, std::make_shared<TensorStrainPolicy<Dimension>>(mStrainAlgorithm));

  // Register the damage.
  // Note we are overriding the default no-op policy for the damage
  // as originally registered by the SolidSPHHydroBase class.
  auto& damage = nodes.damage();
  state.enroll(damage, std::make_shared<TensorDamagePolicy<Dimension>>(*this));

  // Mask out nodes beyond the critical damage threshold from setting the timestep.
  auto maskKey = state.buildFieldKey(HydroFieldNames::timeStepMask, nodes.name());
  auto& mask = state.field(maskKey, 0);
  const auto n = nodes.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    if (damage(i).Trace() > mCriticalDamageThreshold) mask(i) = 0;
  }
}

//------------------------------------------------------------------------------
// Register the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDdamageDt);
}

//------------------------------------------------------------------------------
// Apply the boundary conditions to the ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  const auto nodeListName = this->nodeList().name();
  const auto DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  auto& D = state.field(DKey, SymTensor::zero());

  // Apply ghost boundaries to the damage.
  for (auto boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyGhostBoundary(D);
  }
}

//------------------------------------------------------------------------------
// Enforce boundary conditions for the physics specific fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  const auto nodeListName = this->nodeList().name();
  const auto DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  auto& D = state.field(DKey, SymTensor::zero());

  // Enforce!
  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceBoundary(D);
  }
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Set the moduli.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(SolidFieldNames::bulkModulus, state, derivs);
  updateStateFields(SolidFieldNames::shearModulus, state, derivs);
  updateStateFields(SolidFieldNames::yieldStrength, state, derivs);
  updateStateFields(SolidFieldNames::YoungsModulus, state, derivs);
  updateStateFields(SolidFieldNames::longitudinalSoundSpeed, state, derivs);
}

//------------------------------------------------------------------------------
// Cull the flaws on each node to the single weakest one.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
cullToWeakestFlaws() {
  const auto n = mFlaws.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    auto& flaws = mFlaws[i];
    if (flaws.size() > 0) {
      const auto maxVal = *max_element(flaws.begin(), flaws.end());
      flaws = vector<double>(maxVal);
    }
  }
}

//------------------------------------------------------------------------------
// Compute a Field with the sum of the activation energies per node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
TensorDamageModel<Dimension>::
sumActivationEnergiesPerNode() const {
  auto& nodeList = this->nodeList();
  const auto n = mFlaws.numInternalElements();
  Field<Dimension, Scalar> result("Sum activation energies", nodeList);
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto& flaws = mFlaws(i);
    for (const auto fij: flaws) {
      result(i) += fij;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute a Field with the number of flaws per node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
TensorDamageModel<Dimension>::
numFlawsPerNode() const {
  auto& nodeList = this->nodeList();
  const auto n = mFlaws.numInternalElements();
  Field<Dimension, Scalar> result("num flaws", nodeList);
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    result(i) = flawsForNode(i).size();
  }
  return result;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mFlaws, pathName + "/flaws");
  file.write(mYoungsModulus, pathName + "/YoungsModulus");
  file.write(mLongitudinalSoundSpeed, pathName + "/LongitudinalSoundSpeed");
  file.write(mStrain, pathName + "/strain");
  file.write(mEffectiveStrain, pathName + "/effectiveStrain");
  file.write(mDdamageDt, pathName + "/DdamageDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DamageModel<Dimension>::restoreState(file, pathName);
  file.read(mFlaws, pathName + "/flaws");
  file.read(mYoungsModulus, pathName + "/YoungsModulus");
  file.read(mLongitudinalSoundSpeed, pathName + "/LongitudinalSoundSpeed");
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mDdamageDt, pathName + "/DdamageDt");
}

}

