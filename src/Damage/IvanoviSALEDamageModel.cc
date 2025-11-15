//---------------------------------Spheral++----------------------------------//
// IvanoviSALEDamageModel
//
// The Ivanov damage model, hopefully close to how it's implemented in iSALE.
// This damage model is most appropriate for rocky materials.
//
// Refs:
// 
// Collins, G. S., Melosh, H. J., & Ivanov, B. A. (2004). Modeling damage and deformation in impact simulations.
//   Meteoritics & Planetary Science, 39(2), 217–231. http://doi.wiley.com/10.1111/j.1945-5100.2004.tb00337.x
//
// Raducan, S. D., Davison, T. M., Luther, R., & Collins, G. S. (2019). The role of asteroid strength, porosity and
//   internal friction in impact momentum transfer. Icarus. https://doi.org/10.1016/J.ICARUS.2019.03.040
//
// Lundborg, N. (1967). The strength-size relation of granite. International Journal of Rock Mechanics and
//   Mining Sciences & Geomechanics Abstracts, 4(3):269–272.f
//
// Created by JMO, Sat Jun 26 11:35:44 PDT 2021
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "IvanoviSALEDamageModel.hh"
#include "TensorStrainPolicy.hh"
#include "IvanoviSALEDamagePolicy.hh"
#include "YoungsModulusPolicy.hh"
#include "LongitudinalSoundSpeedPolicy.hh"
#include "DamageGradientPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/updateStateFields.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/Neighbor.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/uniform_random.hh"

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
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
IvanoviSALEDamageModel<Dimension>::
IvanoviSALEDamageModel(SolidNodeList<Dimension>& nodeList,
                       const TableKernel<Dimension>& W,
                       const double minPlasticFailure,
                       const double plasticFailurePressureSlope,
                       const double plasticFailurePressureOffset,
                       const double tensileFailureStress,
                       const double crackGrowthMultiplier,
                       const DamageCouplingAlgorithm damageCouplingAlgorithm,
                       const double criticalDamageThreshold,
                       const Field<Dimension, int>& mask):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mEpsPfb(minPlasticFailure),
  mB(plasticFailurePressureSlope),
  mPc(plasticFailurePressureOffset),
  mTensileFailureStress(tensileFailureStress),
  mCriticalDamageThreshold(criticalDamageThreshold),
  mMask(mask),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mDdamageDt(IvanoviSALEDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
IvanoviSALEDamageModel<Dimension>::
~IvanoviSALEDamageModel() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
//
// In this model we compute the scalar damage derivative assuming unresolved
// crack growth for every point. However, that is not applied in the tensor
// damage update policy unless the flaws are actually activated.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Set the scalar magnitude of the damage evolution.
  const auto* nodeListPtr = &(this->nodeList());
  auto&       DDDt = derivs.field(state.buildFieldKey(IvanoviSALEDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeListPtr->name()), 0.0);
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
typename IvanoviSALEDamageModel<Dimension>::TimeStepType
IvanoviSALEDamageModel<Dimension>::
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
IvanoviSALEDamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Register Youngs modulus and the longitudinal sound speed.
  auto& nodes = this->nodeList();
  state.enroll(mYoungsModulus, std::make_shared<YoungsModulusPolicy<Dimension>>(nodes));
  state.enroll(mLongitudinalSoundSpeed, std::make_shared<LongitudinalSoundSpeedPolicy<Dimension>>(nodes));

  // Register the strain and effective strain.
  state.enroll(mStrain);
  state.enroll(mEffectiveStrain, std::make_shared<TensorStrainPolicy<Dimension>>(TensorStrainAlgorithm::PseudoPlasticStrain));

  // Register the damage and state it requires.
  // Note we are overriding the default no-op policy for the damage
  // as originally registered by the SolidSPHHydroBase class.
  auto& damage = this->nodeList().damage();
  state.enroll(damage, std::make_shared<IvanoviSALEDamagePolicy<Dimension>>(mEpsPfb,
                                                                            mB,
                                                                            mPc,
                                                                            mTensileFailureStress));
 
  // Mask out nodes beyond the critical damage threshold from setting the timestep.
  auto maskKey = state.buildFieldKey(HydroFieldNames::timeStepMask, this->nodeList().name());
  auto& mask = state.field(maskKey, 0);
  const auto nlocal = this->nodeList().numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    if (damage(i).Trace() > mCriticalDamageThreshold) mask(i) = 0;
  }
}

//------------------------------------------------------------------------------
// Register the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDdamageDt);
}

//------------------------------------------------------------------------------
// Apply the boundary conditions to the ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  auto& D = state.field(DKey, SymTensor::zero());

  // Apply ghost boundaries to the damage.
  for (auto boundaryItr = this->boundaryBegin();
       boundaryItr < this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyGhostBoundary(D);
  }
}

//------------------------------------------------------------------------------
// Enforce boundary conditions for the physics specific fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  auto& D = state.field(DKey, SymTensor::zero());

  // Enforce!
  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr < this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceBoundary(D);
  }
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
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
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mYoungsModulus, pathName + "/YoungsModulus");
  file.write(mLongitudinalSoundSpeed, pathName + "/LongitudinalSoundSpeed");
  file.write(mStrain, pathName + "/strain");
  file.write(mEffectiveStrain, pathName + "/effectiveStrain");
  file.write(mDdamageDt, pathName + "/DdamageDt");
  file.write(mMask, pathName + "/mask");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DamageModel<Dimension>::restoreState(file, pathName);
  file.read(mYoungsModulus, pathName + "/YoungsModulus");
  file.read(mLongitudinalSoundSpeed, pathName + "/LongitudinalSoundSpeed");
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mDdamageDt, pathName + "/DdamageDt");
  file.read(mMask, pathName + "/mask");
}

}

