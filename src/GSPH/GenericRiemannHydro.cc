//---------------------------------Spheral++----------------------------------//
// GenericRiemannHydro -- pure virtual class for hydros using a Riemann
//                        solver
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"
#include "Physics/Physics.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/PureReplaceState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/updateStateFields.hh"

#include "Hydro/computeSPHVolume.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"

#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/GenericRiemannHydro.hh"
#include "GSPH/initializeGradients.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <limits.h>
#include <fstream>
#include <sstream>

using std::vector;
using std::string;
using std::pair;
using std::to_string;
using std::make_pair;
using std::make_shared;


namespace {
// Provide to_string for Spheral Vector
template<typename Vector>
std::string
vec_to_string(const Vector& vec) {
  std::ostringstream oss;
  oss << vec << std::endl;
  return oss.str();
}
}

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GenericRiemannHydro<Dimension>::
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
                    const Vector& xmax):
  Physics<Dimension>(),
  mRestart(registerWithRestart(*this)),
  mRiemannSolver(riemannSolver),
  mKernel(W),
  mGradientType(gradType),
  mDensityUpdate(densityUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mCorrectVelocityGradient(correctVelocityGradient),
  mUseVelocityMagnitudeForDt(useVelocityMagnitudeForDt),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mSpecificThermalEnergyDiffusionCoefficient(epsDiffusionCoeff),
  mCfl(cfl),
  mxmin(xmin),
  mxmax(xmax),
  mTimeStepMask(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),                       // move up one layer
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),   // move up one layer
  mDrhoDx(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mRiemannDpDx(FieldStorageType::CopyFields),
  mRiemannDvDx(FieldStorageType::CopyFields),
  mNewRiemannDpDx(FieldStorageType::CopyFields),
  mNewRiemannDvDx(FieldStorageType::CopyFields),
  mPairAccelerationsPtr(),
  mPairDepsDtPtr() {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
  mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::XSPHDeltaV);
  mM = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::M_SPHCorrection);
  mDxDt = dataBase.newFluidFieldList(Vector::zero(), IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::hydroAcceleration);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDrhoDx = dataBase.newFluidFieldList(Vector::zero(),GSPHFieldNames::densityGradient);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::velocityGradient);
  mRiemannDpDx = dataBase.newFluidFieldList(Vector::zero(),GSPHFieldNames::RiemannPressureGradient);
  mRiemannDvDx = dataBase.newFluidFieldList(Tensor::zero(),GSPHFieldNames::RiemannVelocityGradient);
  mNewRiemannDpDx = dataBase.newFluidFieldList(Vector::zero(),ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannPressureGradient);
  mNewRiemannDvDx = dataBase.newFluidFieldList(Tensor::zero(),ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannVelocityGradient);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  const auto& connectivityMap = dataBase.connectivityMap();
  const auto mass = dataBase.fluidMass();
  const auto massDensity = dataBase.fluidMassDensity();
  const auto position = dataBase.fluidPosition();
  const auto H = dataBase.fluidHfield();
  auto velocity = dataBase.fluidVelocity();
  
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(HydroFieldNames::soundSpeed, state, derivs);

  computeSPHVolume(mass,massDensity,mVolume);

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
        boundItr != this->boundaryEnd();
        ++boundItr){
    (*boundItr)->applyFieldListGhostBoundary(mVolume);
    (*boundItr)->applyFieldListGhostBoundary(velocity);
    (*boundItr)->applyFieldListGhostBoundary(mPressure);
  }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  initializeGradients(connectivityMap,
                      this->kernel(),
                      position,
                      H,
                      mVolume,
                      mPressure,
                      velocity,
                      mM,
                      mRiemannDpDx,
                      mRiemannDvDx);
 
}


//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mRiemannDpDx, Vector::zero(), GSPHFieldNames::RiemannPressureGradient, false);
  dataBase.resizeFluidFieldList(mRiemannDvDx, Tensor::zero(), GSPHFieldNames::RiemannVelocityGradient, false);

  auto mass = dataBase.fluidMass();
  auto massDensity = dataBase.fluidMassDensity();
  auto position = dataBase.fluidPosition();
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  auto velocity = dataBase.fluidVelocity();

  auto positionPolicy = make_policy<IncrementState<Dimension, Vector>>();
  auto pressurePolicy = make_policy<PressurePolicy<Dimension>>();
  auto csPolicy = make_policy<SoundSpeedPolicy<Dimension>>();
  auto pressureGradientPolicy = make_policy<PureReplaceState<Dimension,Vector>>();
  auto velocityGradientPolicy = make_policy<PureReplaceState<Dimension,Tensor>>();
  auto velocityPolicy = make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::position,HydroFieldNames::specificThermalEnergy},true);

  // normal state variables
  state.enroll(mTimeStepMask);
  state.enroll(mVolume);
  state.enroll(mass);
  state.enroll(massDensity);
  state.enroll(position, positionPolicy);
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);
  state.enroll(velocity, velocityPolicy);

  if (mRiemannSolver.linearReconstruction()){
    state.enroll(mRiemannDpDx, pressureGradientPolicy);
    state.enroll(mRiemannDvDx, velocityGradientPolicy);
  }else{
    state.enroll(mRiemannDpDx);
    state.enroll(mRiemannDvDx);
  }

  // conditional for energy method
  if (mCompatibleEnergyEvolution) {
    auto thermalEnergyPolicy = make_policy<CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>>(dataBase);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }else if (mEvolveTotalEnergy) {
    auto thermalEnergyPolicy = make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>();
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  } else {
    auto thermalEnergyPolicy = make_policy<IncrementState<Dimension, Scalar>>();
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }
  
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  
  // Create the scratch fields.
  dataBase.resizeFluidFieldList(mNewRiemannDpDx, Vector::zero(), ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannPressureGradient, false);
  dataBase.resizeFluidFieldList(mNewRiemannDvDx, Tensor::zero(), ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannVelocityGradient, false);
  dataBase.resizeFluidFieldList(mNormalization, 0.0, HydroFieldNames::normalization, false);
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero(), HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero(), HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero(), HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mDrhoDx, Vector::zero(), GSPHFieldNames::densityGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero(), HydroFieldNames::M_SPHCorrection, false);
  
  // Check if someone already registered DxDt.
  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero(), IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }
  // Check that no-one else is trying to control the hydro vote for DvDt.
  CHECK(not derivs.registered(mDvDt));
  derivs.enroll(mDrhoDx);
  derivs.enroll(mNewRiemannDpDx);
  derivs.enroll(mNewRiemannDvDx);
  derivs.enroll(mDvDt);
  derivs.enroll(mNormalization);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mXSPHDeltaV);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mM);

  if (mCompatibleEnergyEvolution) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    mPairDepsDtPtr = std::make_unique<PairWorkType>(connectivityMap);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
    derivs.enroll(HydroFieldNames::pairWork, *mPairDepsDtPtr);
  }
}

//------------------------------------------------------------------------------
// GSPH Time step
//------------------------------------------------------------------------------
template<typename Dimension>
typename GenericRiemannHydro<Dimension>::TimeStepType
GenericRiemannHydro<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar /*currentTime*/) const{
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  // Get some useful fluid variables from the DataBase.
  const auto  mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto  cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
                                                         this->requireOverlapConnectivity(),
                                                         this->requireIntersectionConnectivity());
  const auto  numNodeLists = connectivityMap.nodeLists().size();

  // Initialize the return value to some impossibly high value.
  auto minDt = make_pair(std::numeric_limits<double>::max(), string());

  // Loop over every fluid node.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& fluidNodeList = **(dataBase.fluidNodeListBegin() + nodeListi);
    const auto nPerh = fluidNodeList.nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);

    // Walk all the nodes in this FluidNodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
    // #pragma omp parallel for reduction(MINPAIR:minDt)
#pragma omp parallel
    {
      auto minDt_local = minDt;
#pragma omp for
      for (auto k = 0u; k < ni; ++k) {
        const auto i = connectivityMap.ithNode(nodeListi, k);

        // If this node is masked, don't worry about it.
        if (mask(nodeListi, i) == 1) {

          // Get this nodes minimum characteristic smoothing scale.
          CHECK2(H(nodeListi, i).Determinant() >  0.0,
                 "Bad H tensor : " << H(nodeListi, i) << " : " << fluidNodeList.name() << " " << i << " " << fluidNodeList.firstGhostNode());
          const auto& Hi = H(nodeListi, i);
          const Scalar nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;
          // const Scalar nodeScalei = 1.0/Dimension::rootnu(Hi.Determinant());
          //     const Scalar nodeScalei = nodeExtent(nodeListi, i).minElement()/kernelExtent;

          // Sound speed limit.
          const auto csi = cs(nodeListi, i);
          const auto csDt = nodeScalei/(csi + tiny);
          if (csDt < minDt_local.first) {
            minDt_local = make_pair(csDt, ("Sound speed limit: dt = " + to_string(csDt) + "\n" +
                                           "                   cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                           "            nodeScale = " + to_string(nodeScalei) + "\n" +
                                           "             material = " + fluidNodeList.name() + "\n" +
                                           "(nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                           "           @ position = " + vec_to_string(position(nodeListi, i))));
          }

          // Velocity divergence limit.
          const auto divVelocity = DvDx(nodeListi, i).Trace();
          const auto divvDt = 1.0/(std::abs(divVelocity) + tiny);
          if (divvDt < minDt_local.first) {
            minDt_local = make_pair(divvDt, ("Velocity divergence limit: dt = " + to_string(divvDt) + "\n" +
                                             "                 div velocity = " + to_string(divVelocity) + "\n" +
                                             "                     material = " + fluidNodeList.name() + "\n" +
                                             "        (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                             "                   @ position = " + vec_to_string(position(nodeListi, i))));
          }

          // Maximum velocity difference limit.
          //const auto& xi = position(nodeListi, i);
          const auto& vi = velocity(nodeListi, i);
          const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
          for (auto nodeListj = 0u; nodeListj != numNodeLists; ++nodeListj) {
            const auto& connectivity = fullConnectivity[nodeListj];
            for (auto jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) {
              const auto  j = *jItr;
              //const auto& xj = position(nodeListj, j);
              const auto& vj = velocity(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              //const auto  rhoj = rho(nodeListj, j);
              const Scalar nodeScalej = 1.0/Hj.eigenValues().maxElement()/nPerh;
              // const Scalar vij = std::abs((vj - vi).dot((xi - xj).unitVector()));
              const auto  vij = vi - vj;
              const auto  dtVelDiff = std::min(nodeScalei, nodeScalej)*safeInvVar(vij.magnitude(), 1e-30);
              if (dtVelDiff < minDt_local.first) {
                minDt_local = make_pair(dtVelDiff, ("Pairwise velocity difference limit: dt = " + to_string(dtVelDiff) + "\n" + 
                                                    "                              material = " + fluidNodeList.name() + "\n" +
                                                    "                  (nodeListi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                                    "                  (nodeListj, i, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(Process::getRank()) + ")\n" +
                                                    "                   @ pos(nodeListi, i) = " + vec_to_string(position(nodeListi, i)) + "\n" +
                                                    "                   @ pos(nodeListj, j) = " + vec_to_string(position(nodeListj, j)) + "\n" +
                                                    "                                   vij = " + to_string(vij.magnitude()) + "\n" +
                                                    "                            nodeScalei = " + to_string(nodeScalei) + "\n" +
                                                    "                            nodeScalej = " + to_string(nodeScalej)));
              }

              
            }
          }


          // Total acceleration limit.
          const auto vmagi = vi.magnitude();
          const auto dtAcc = 0.1*std::max(nodeScalei/(vmagi + tiny), vmagi/(DvDt(nodeListi, i).magnitude() + tiny));
          if (dtAcc < minDt_local.first) {
            minDt_local = make_pair(dtAcc, ("Total acceleration limit: dt = " + to_string(dtAcc) + "\n" + 
                                            "              |acceleration| = " + to_string(DvDt(nodeListi, i).magnitude()) + "\n" +
                                            "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                    material = " + fluidNodeList.name() + "\n" +
                                            "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                            "                  @ position = " + vec_to_string(position(nodeListi, i))));
          }

          // If requested, limit against the absolute velocity.
          if (useVelocityMagnitudeForDt()) {
            const auto velDt = nodeScalei/(velocity(nodeListi, i).magnitude() + 1.0e-10);
            if (velDt < minDt_local.first) {
              minDt_local = make_pair(velDt, ("Velocity magnitude limit: dt = " + to_string(velDt) + "\n" +
                                              "                        |vi| = " + to_string(velocity(nodeListi, i).magnitude()) + "\n" +
                                              "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                              "                    material = " + fluidNodeList.name() + "\n" +
                                              "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                              "                  @ position = " + vec_to_string(position(nodeListi, i))));
            }
          }
        }
      }

#pragma omp critical
      if (minDt_local.first < minDt.first) minDt = minDt_local;
    }
  }
  // Scale by the cfl safety factor.
  minDt.first *= cfl();
  return minDt;

}
//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// GenericRiemannHydro<Dimension>::
// preStepInitialize(const DataBase<Dimension>& dataBase, 
//                   State<Dimension>& state,
//                   StateDerivatives<Dimension>& /*derivs*/) {
//   TIME_GSPHpreStepInitialize.start();

//   TIME_GSPHpreStepInitialize.stop();
// }

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GenericRiemannHydro<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  return false;
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto DpDx = state.fields(GSPHFieldNames::RiemannPressureGradient,Vector::zero()); 
  auto DvDx = state.fields(GSPHFieldNames::RiemannVelocityGradient,Tensor::zero()); 

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(DpDx);
    (*boundaryItr)->applyFieldListGhostBoundary(DvDx);

  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto DpDx = state.fields(GSPHFieldNames::RiemannPressureGradient,Vector::zero()); 
  auto DvDx = state.fields(GSPHFieldNames::RiemannVelocityGradient,Tensor::zero()); 

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(volume);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(DpDx);
    (*boundaryItr)->enforceFieldListBoundary(DvDx);
  }
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mVolume, pathName + "/volume");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");

  file.write(mNormalization, pathName + "/normalization");
  file.write(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  // time derivs
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");

  // spatial derivs
  file.write(mM, pathName + "/M");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mRiemannDvDx, pathName + "/riemannDvDx");
  file.write(mRiemannDpDx, pathName + "/riemannDpDx");
  file.write(mNewRiemannDvDx, pathName + "/newRiemannDvDx");
  file.write(mNewRiemannDpDx, pathName + "/newRiemannDpDx");
  

}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mVolume, pathName + "/volume");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");

  file.read(mNormalization, pathName + "/normalization");
  file.read(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  // time derivs
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");

  // spatial derivs
  file.read(mM, pathName + "/M");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mRiemannDvDx, pathName + "/riemannDvDx");
  file.read(mRiemannDpDx, pathName + "/riemannDpDx");
  file.read(mNewRiemannDvDx, pathName + "/newRiemannDvDx");
  file.read(mNewRiemannDpDx, pathName + "/newRiemannDpDx");

}

}
