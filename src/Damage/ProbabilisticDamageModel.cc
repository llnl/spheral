//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamageModel
// A damage model based on Weibull statistics that uses volume based
// probabilities per node to decide when damage starts to accrue.  Should
// generate similar results to the classic Benz-Asphaug (Grady-Kipp) model
// without generating explicit flaws.  Also appropriate for use with varying
// resolution materials.
//
// Created by JMO, Tue Apr 13 15:58:08 PDT 2021
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "ProbabilisticDamageModel.hh"
#include "TensorStrainPolicy.hh"
#include "ProbabilisticDamagePolicy.hh"
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
#include "Distributed/allReduce.hh"
#include "Utilities/uniform_random.hh"

#include <boost/functional/hash.hpp>  // hash_combine

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
ProbabilisticDamageModel<Dimension>::
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
                         const Field<Dimension, int>& mask):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mStrainAlgorithm(strainAlgorithm),
  mDamageInCompression(damageInCompression),
  mkWeibull(kWeibull),
  mmWeibull(mWeibull),
  mVolumeMultiplier(volumeMultiplier),
  mVmin(std::numeric_limits<double>::max()),
  mVmax(std::numeric_limits<double>::min()),
  mCriticalDamageThreshold(criticalDamageThreshold),
  mSeed(seed),
  mMinFlawsPerNode(minFlawsPerNode),
  mNumFlaws(SolidFieldNames::numFlaws, nodeList),
  mMask(mask),
  mMinFlaw(SolidFieldNames::minFlaw, nodeList),
  mMaxFlaw(SolidFieldNames::maxFlaw, nodeList),
  mInitialVolume(SolidFieldNames::initialVolume, nodeList),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mDdamageDt(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamageModel<Dimension>::
~ProbabilisticDamageModel() {
}

//------------------------------------------------------------------------------
// initializeProblemStartupDependencies
//
// After all initial state has been initialize (node positions, masses, etc),
// but before we try to run any physics cycles.  This is when we initialize a
// lot of our state for the damage model.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // How many points are actually being damaged?
  // We have to be careful to use an unsigned size_t here due to overflow
  // problems with large numbers of points.
  size_t nused_local = 0u;
  for (auto i = 0u; i < mMask.numInternalElements(); ++i) {
    if (mMask[i] == 1) ++nused_local;
  }
  const size_t nused_global = allReduce(nused_local, SPHERAL_OP_SUM);

  // Compute the Morton-ordering for hashing with the global seed to seed each
  // point-wise random number generator.
  typedef KeyTraits::Key Key;
  const FieldList<Dimension, Key> keyList = mortonOrderIndices(dataBase);
  const auto& nodes = this->nodeList();
  CHECK(keyList.fieldForNodeList(nodes) < keyList.end());
  const Field<Dimension, Key>& keys = **(keyList.fieldForNodeList(nodes));
  
  // Compute the initial volumes and random seeds for each node.  We hash the
  // morton index of each point with the global seed value to create unique
  // but reproducible seeds for each points random number generator.
  auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodes.name()); };
  const auto& mass = state.field(buildKey(HydroFieldNames::mass), 0.0);
  const auto& rho = (state.registered(buildKey(SolidFieldNames::porositySolidDensity)) ?
                     state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                     state.field(buildKey(HydroFieldNames::massDensity), 0.0));
  const auto  nlocal = nodes.numInternalNodes();
  vector<uniform_random> randomGenerators(nlocal);
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    if (mMask(i) == 1) {
      CHECK(mass(i) > 0.0 and rho(i) > 0.0);
      mInitialVolume(i) = mass(i)/rho(i) * mVolumeMultiplier;
      mVmin = std::min(mVmin, mInitialVolume(i));
      mVmax = std::max(mVmax, mInitialVolume(i));
#ifdef __APPLE__
      std::size_t seedi = mSeed;
#else
      Key seedi = mSeed;
#endif
      boost::hash_combine(seedi, keys(i));
      randomGenerators[i].seed(seedi);      // starting out generating in [0,1)
      randomGenerators[i]();                // Recommended to discard first value in sequence
    }
  }
  mVmin = allReduce(mVmin, SPHERAL_OP_MIN);
  mVmax = allReduce(mVmax, SPHERAL_OP_MAX);

  // Generate min/max ranges of flaws for each point.
  const auto mInv = 1.0/mmWeibull;
  size_t minNumFlaws = std::numeric_limits<size_t>::max();
  size_t maxNumFlaws = 0u;
  size_t totalNumFlaws = 0u;
  auto epsMin = std::numeric_limits<double>::max();
  auto epsMax = std::numeric_limits<double>::min();
  auto numFlawsRatio = 0.0;
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    if (mMask(i) == 1) {
      const auto Nflaws = size_t(mMinFlawsPerNode*mInitialVolume(i)/mVmin + 0.5);  // Target number of flaws on this point
      CHECK(Nflaws >= mMinFlawsPerNode);
      const auto Ai = pow(Nflaws/(mkWeibull*mInitialVolume(i)), mInv);
      mMinFlaw(i) = 1.0;
      mMaxFlaw(i) = 0.0;
      auto ntries = 0u;
      while (mMinFlaw(i) > mMaxFlaw(i) and ntries++ < 100u) {
        mMinFlaw(i) = Ai*pow(1.0 - pow(randomGenerators[i](), 1.0/Nflaws), mInv);
        mMaxFlaw(i) = Ai*pow(randomGenerators[i](), 1.0/(mmWeibull*Nflaws));
      }
      CHECK(mMaxFlaw(i) > mMinFlaw(i));
      mNumFlaws(i) = std::max(1, 1 + int(mInitialVolume(i)*mkWeibull*(pow(mMaxFlaw(i), mmWeibull) - pow(mMinFlaw(i), mmWeibull))));

      // Gather statistics
#pragma omp critical
      {
        minNumFlaws = min(minNumFlaws, size_t(mNumFlaws(i)));
        maxNumFlaws = max(maxNumFlaws, size_t(mNumFlaws(i)));
        totalNumFlaws += size_t(mNumFlaws(i));
        epsMin = std::min(epsMin, mMinFlaw(i));
        epsMax = std::max(epsMax, mMaxFlaw(i));
        numFlawsRatio += double(mNumFlaws(i))/Nflaws;
      }
    }
  }

  // Some diagnostic output.
  if (nused_global > 0) {
    minNumFlaws = allReduce(minNumFlaws, SPHERAL_OP_MIN);
    maxNumFlaws = allReduce(maxNumFlaws, SPHERAL_OP_MAX);
    totalNumFlaws = allReduce(totalNumFlaws, SPHERAL_OP_SUM);
    epsMin = allReduce(epsMin, SPHERAL_OP_MIN);
    epsMax = allReduce(epsMax, SPHERAL_OP_MAX);
    numFlawsRatio = allReduce(numFlawsRatio, SPHERAL_OP_SUM)/nused_global;
    if (Process::getRank() == 0) {
      cerr << "ProbabilisticDamageModel for " << nodes.name() << ":" << endl
           << " Min, max, max/min volumes: " << mVmin << " " << mVmax << " " << mVmax*safeInv(mVmin) << endl
           << "    Min num flaws per node: " << minNumFlaws << endl
           << "    Max num flaws per node: " << maxNumFlaws << endl
           << "    Total num flaws       : " << totalNumFlaws << endl
           << "    Avg flaws per node    : " << totalNumFlaws/nused_global << endl
           << "    Min flaw strain       : " << epsMin << endl
           << "    Max flaw strain       : " << epsMax << endl
           << "    Avg Neff/Nflaws       : " << numFlawsRatio << endl;
    }
  }

  // Set the moduli.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(SolidFieldNames::bulkModulus, state, derivs);
  updateStateFields(SolidFieldNames::shearModulus, state, derivs);
  updateStateFields(SolidFieldNames::yieldStrength, state, derivs);
  updateStateFields(SolidFieldNames::YoungsModulus, state, derivs);
  updateStateFields(SolidFieldNames::longitudinalSoundSpeed, state, derivs);
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
ProbabilisticDamageModel<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Set the scalar magnitude of the damage evolution.
  const auto* nodeListPtr = &(this->nodeList());
  auto&       DDDt = derivs.field(state.buildFieldKey(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeListPtr->name()), 0.0);
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
typename ProbabilisticDamageModel<Dimension>::TimeStepType
ProbabilisticDamageModel<Dimension>::
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
ProbabilisticDamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Register Youngs modulus and the longitudinal sound speed.
  auto& nodes = this->nodeList();
  state.enroll(mYoungsModulus, std::make_shared<YoungsModulusPolicy<Dimension>>(nodes));
  state.enroll(mLongitudinalSoundSpeed, std::make_shared<LongitudinalSoundSpeedPolicy<Dimension>>(nodes));

  // Register the strain and effective strain.
  state.enroll(mStrain);
  state.enroll(mEffectiveStrain, std::make_shared<TensorStrainPolicy<Dimension>>(mStrainAlgorithm));

  // Register the damage and state it requires.
  // Note we are overriding the default no-op policy for the damage
  // as originally registered by the SolidSPHHydroBase class.
  auto& damage = nodes.damage();
  state.enroll(damage, std::make_shared<ProbabilisticDamagePolicy<Dimension>>(mDamageInCompression,
                                                                              mkWeibull,
                                                                              mmWeibull));
  state.enroll(mNumFlaws);
  state.enroll(mMinFlaw);
  state.enroll(mMaxFlaw);
  state.enroll(mInitialVolume);

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
ProbabilisticDamageModel<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDdamageDt);
}

//------------------------------------------------------------------------------
// Apply the boundary conditions to the ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
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
ProbabilisticDamageModel<Dimension>::
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
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mNumFlaws, pathName + "/numFlaws");
  file.write(mMinFlaw, pathName + "/minFlaw");
  file.write(mMaxFlaw, pathName + "/maxFlaw");
  file.write(mYoungsModulus, pathName + "/YoungsModulus");
  file.write(mLongitudinalSoundSpeed, pathName + "/LongitudinalSoundSpeed");
  file.write(mDdamageDt, pathName + "/DdamageDt");
  file.write(mStrain, pathName + "/strain");
  file.write(mEffectiveStrain, pathName + "/effectiveStrain");
  file.write(mMask, pathName + "/mask");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DamageModel<Dimension>::restoreState(file, pathName);
  file.read(mNumFlaws, pathName + "/numFlaws");
  file.read(mMinFlaw, pathName + "/minFlaw");
  file.read(mMaxFlaw, pathName + "/maxFlaw");
  file.read(mYoungsModulus, pathName + "/YoungsModulus");
  file.read(mLongitudinalSoundSpeed, pathName + "/LongitudinalSoundSpeed");
  file.read(mDdamageDt, pathName + "/DdamageDt");
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mMask, pathName + "/mask");
}

}

