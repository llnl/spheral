//---------------------------------Spheral++----------------------------------//
// SPH -- The classic SPH/ASPH hydrodynamic packages for Spheral++.
// 
// Created by JMO, Thu Nov 21 16:36:40 PST 2024
//----------------------------------------------------------------------------//
#include "SPH/SPHBase.hh"

#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "correctSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/updateStateFields.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosityView.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/range.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
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
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHBase<Dimension>::
SPHBase(DataBase<Dimension>& dataBase,
        ArtificialViscosity<Dimension>& Q,
        const TableKernel<Dimension>& W,
        const TableKernel<Dimension>& WPi,
        const double cfl,
        const bool useVelocityMagnitudeForDt,
        const bool compatibleEnergyEvolution,
        const bool evolveTotalEnergy,
        const bool gradhCorrection,
        const bool XSPH,
        const bool correctVelocityGradient,
        const bool sumMassDensityOverAllNodeLists,
        const MassDensityType densityUpdate,
        const double epsTensile,
        const double nTensile,
        const Vector& xmin,
        const Vector& xmax):
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mKernel(W),
  mPiKernel(WPi),
  mDensityUpdate(densityUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mGradhCorrection(gradhCorrection),
  mXSPH(XSPH),
  mCorrectVelocityGradient(correctVelocityGradient),
  mSumMassDensityOverAllNodeLists(sumMassDensityOverAllNodeLists),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mxmin(xmin),
  mxmax(xmax),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mOmegaGradh(FieldStorageType::CopyFields),
  mEntropy(FieldStorageType::CopyFields),
  mMassDensityCorrection(FieldStorageType::CopyFields),
  mMassDensitySum(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mGradRho(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mLocalM(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mOmegaGradh = dataBase.newFluidFieldList(1.0, HydroFieldNames::omegaGradh);
  mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  mMassDensityCorrection = dataBase.newFluidFieldList(0.0, HydroFieldNames::massDensityCorrection);
  mMassDensitySum = dataBase.newFluidFieldList(0.0, ReplaceState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity);
  mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
  mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mGradRho = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::massDensityGradient);
  mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
  mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SPHBaseInitializeStartup");

  // Set the moduli.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(HydroFieldNames::soundSpeed, state, derivs);

  // If we're using the grad-h corrections they need to be initialized.  Our
  // postStateUpdate already does this.
  this->postStateUpdate(0.0, 1.0, dataBase, state, derivs);

  // dataBase.fluidEntropy(mEntropy);

  // // In some cases we need the volume per node as well.
  // const bool updateVolume = (this->densityUpdate() == VoronoiCellDensity or
  //                            this->densityUpdate() == SumVoronoiCellDensity);
  // if (updateVolume) {
  //   Mesh<Dimension> mesh;
  //   NodeList<Dimension> voidNodes("void", 0, 0);
  //   vector<const NodeList<Dimension>*> nodeLists(dataBase.nodeListBegin(), dataBase.nodeListEnd());
  //   nodeLists.push_back(&voidNodes);
  //   sort(nodeLists.begin(), nodeLists.end(), typename NodeListRegistrar<Dimension>::NodeListComparator());
  //   Vector xmin, xmax;
  //   boundingBox(dataBase.fluidPosition(), xmin, xmax);
  //   generateMesh<Dimension, 
  //                           typename vector<const NodeList<Dimension>*>::iterator,
  //                           ConstBoundaryIterator>(nodeLists.begin(), nodeLists.end(),
  //                                                  this->boundaryBegin(), this->boundaryEnd(),
  //                                                  xmin, xmax, 
  //                                                  true,           // generateVoid
  //                                                  false,          // generateParallelConnectivity
  //                                                  false,          // remove boundary zones
  //                                                  2.0,            // voidThreshold
  //                                                  mesh,
  //                                                  voidNodes);

  //   mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  //   for (unsigned nodeListi = 0; nodeListi != dataBase.numFluidNodeLists(); ++nodeListi) {
  //     const unsigned n = mVolume[nodeListi]->numInternalElements();
  //     const unsigned offset = mesh.offset(nodeListi);
  //     for (unsigned i = 0; i != n; ++i) {
  //       mVolume(nodeListi, i) = mesh.zone(offset + i).volume();
  //     }
  //   }
  // }
  TIME_END("SPHBaseInitializeStartup");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("SPHBaseRegister");

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mOmegaGradh, 1.0, HydroFieldNames::omegaGradh);

  // We may need the volume per node as well.
  const bool updateVolume = (this->densityUpdate() == MassDensityType::VoronoiCellDensity or
                             this->densityUpdate() == MassDensityType::SumVoronoiCellDensity);
  if (updateVolume) {
    dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  }

  // Now register away.
  // Mass.
  auto mass = dataBase.fluidMass();
  state.enroll(mass);

  // Mass density
  auto massDensity = dataBase.fluidMassDensity();
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    state.enroll(*massDensity[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>(fluidNodeListPtr->rhoMin(),
                                                                                                fluidNodeListPtr->rhoMax()));
  }

  // Volume.
  if (updateVolume) state.enroll(mVolume);

  // Register the position update.
  // We make this dependent on specific thermal energy to cover cases like the RZ specialization
  auto position = dataBase.fluidPosition();
  state.enroll(position, make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::specificThermalEnergy}));

  // Register the velocity
  // We make this dependent on the thermal energy in case we're using the compatible energy update
  auto velocity = dataBase.fluidVelocity();
  state.enroll(velocity, make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::position,
                                                                         HydroFieldNames::specificThermalEnergy},
                                                                        true));  // Use all DvDt sources (wildcard)

  // Register the time step mask, initialized to 1 so that everything defaults to being
  // checked.
  state.enroll(mTimeStepMask);

  // Register the pressure and sound speed.
  state.enroll(mPressure, make_policy<PressurePolicy<Dimension>>());
  state.enroll(mSoundSpeed, make_policy<SoundSpeedPolicy<Dimension>>());

  // Register the grad h correction terms
  // We deliberately make this non-dynamic here.  These corrections are computed
  // during SPHBase::initialize, not as part of our usual state update.
  state.enroll(mOmegaGradh);
  TIME_END("SPHBaseRegister");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SPHBaseRegisterDerivs");

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mMassDensityCorrection, 0.0, HydroFieldNames::massDensityCorrection, false);
  dataBase.resizeFluidFieldList(mMassDensitySum, 0.0, ReplaceState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mNormalization, 0.0, HydroFieldNames::normalization, false);
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mGradRho, Vector::zero, HydroFieldNames::massDensityGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero, HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection, false);
  derivs.enroll(mMassDensityCorrection);
  derivs.enroll(mMassDensitySum);
  derivs.enroll(mNormalization);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mXSPHDeltaV);

  // Check if someone already registered DxDt.
  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }

  // Check that no-one else is trying to control the hydro vote for DvDt.
  CHECK(not derivs.registered(mDvDt));
  derivs.enroll(mDvDt);

  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enroll(mGradRho);
  derivs.enroll(mM);
  derivs.enroll(mLocalM);
  TIME_END("SPHBaseRegisterDerivs");
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SPHBasePreStepInitialize");

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  switch(densityUpdate()) {

  case MassDensityType::IntegrateDensity:
    break;

  case MassDensityType::RigorousSumDensity:
  case MassDensityType::CorrectedSumDensity:
    {
      const auto& connectivityMap = state.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
      if (densityUpdate() == MassDensityType::CorrectedSumDensity) {
        correctSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
        for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensity);
        for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
      }
    }
    break;

  case MassDensityType::SumDensity:
    {
      auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      const auto massDensitySum = derivs.fields(ReplaceState<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                HydroFieldNames::massDensity, 0.0);
      massDensity.assignFields(massDensitySum);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
    }
    break;

  case MassDensityType::HybridSumDensity:
    {
      auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      const auto massDensitySum = derivs.fields(ReplaceState<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                HydroFieldNames::massDensity, 0.0);
      const auto normalization = state.fields(HydroFieldNames::normalization, 0.0);
      const auto numNodeLists = normalization.size();
      for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
        const auto n = normalization[nodeListi]->numInternalElements();
        for (auto i = 0u; i < n; ++i) {
          if (normalization(nodeListi, i) > 0.95) massDensity(nodeListi, i) = massDensitySum(nodeListi, i);
        }
      }
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
    }
    break;

  case MassDensityType::VoronoiCellDensity:
    {
      this->updateVolume(state, false);
      const auto mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto volume = state.fields(HydroFieldNames::volume, 0.0);
      auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      massDensity = mass / volume;
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
    }
    break;

  case MassDensityType::SumVoronoiCellDensity:
    {
      this->updateVolume(state, true);
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  volume = state.fields(HydroFieldNames::volume, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeSumVoronoiCellMassDensity(connectivityMap, this->kernel(), position, mass, volume, H, massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensity);
      for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
    }
    break;

  default:
    break;
  }

  // // This form looks for points that are too close based on specific volume.
  // if (mfilter > 0.0) {
  //   const TableKernel<Dimension>& W = this->kernel();
  //   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  //   const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  //   const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  //   const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //   const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   const unsigned numNodeLists = mass.size();
  //   const Scalar W0 = W.kernelValue(0.0, 1.0);
  //   FieldList<Dimension, Vector> deltar = dataBase.newFluidFieldList(Vector::zero, "delta position");
  //   FieldList<Dimension, Scalar> deltav = dataBase.newFluidFieldList(0.0, "delta velocity");
  //   FieldList<Dimension, Scalar> weightsum = dataBase.newFluidFieldList(0.0, "delta velocity weight sum");
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const Scalar nPerh = position[nodeListi]->nodeListPtr()->nodesPerSmoothingScale();
  //     for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
  //          iItr != connectivityMap.end(nodeListi);
  //          ++iItr) {
  //       const int i = *iItr;
  //       const Vector& ri = position(nodeListi, i);
  //       const Vector& vi = velocity(nodeListi, i);
  //       const Scalar mi = mass(nodeListi, i);
  //       const Scalar rhoi = massDensity(nodeListi, i);
  //       const SymTensor& Hi = H(nodeListi, i);
  //       const SymTensor Hinvi = Hi.Inverse();
  //       const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
  //       for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
  //         for (typename vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
  //              jItr != fullConnectivity[nodeListj].end();
  //              ++jItr) {
  //           const unsigned j = *jItr;
  //           const Vector& rj = position(nodeListj, j);
  //           const Vector& vj = velocity(nodeListj, j);
  //           const Scalar mj = mass(nodeListj, j);
  //           const Scalar rhoj = massDensity(nodeListj, j);
  //           const SymTensor& Hj = H(nodeListj, j);
  //           const SymTensor Hinvj = Hj.Inverse();
  //           const Vector rji = rj - ri;
  //           const Vector rjihat = rji.unitVector();
  //           const Vector etai = Hi*rji;
  //           const Vector etaj = Hj*rji;
  //           const Scalar etaMagi = etai.magnitude();
  //           const Scalar etaMagj = etaj.magnitude();
  //           const Vector delta = 0.5*(max(0.0, 1.0/nPerh - etaMagi)*Hinvi + max(0.0, 1.0/nPerh - etaMagj)*Hinvj)*rjihat;
  //           const Scalar weight = 0.5*(W.kernelValue(etaMagi, 1.0) + W.kernelValue(etaMagj, 1.0))/W0 * (vj - vi).magnitude();
  //           deltar(nodeListi, i) -= weight*delta;
  //           weightsum(nodeListi, i) += weight;
  //           deltav(nodeListi, i) += weight*(vj - vi).magnitude();
  //         }
  //       }
  //     }
  //   }

  //   // Apply the filtering.
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       // const Scalar hi = 1.0/(H(nodeListi, i).eigenValues().maxElement());
  //       // const Scalar mag0 = DvDx(nodeListi, i).eigenValues().maxAbsElement()*hi*dt;
  //       // const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       const Scalar mag0 = deltav(nodeListi, i)*safeInv(weightsum(nodeListi, i))*dt;
  //       if (mag0 > 0.0) {
  //         const Scalar deltamag = deltar(nodeListi, i).magnitude();
  //         // const Scalar effmag = mfilter*deltamag;
  //         const Scalar effmag = mfilter*min(mag0, deltamag);
  //         position(nodeListi, i) += effmag*deltar(nodeListi, i).unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }
  TIME_END("SPHBasePreStepInitialize");
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("SPHBaseFinalizeDerivs");

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    auto DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
      boundaryPtr->applyFieldListGhostBoundary(accelerations);
      boundaryPtr->applyFieldListGhostBoundary(DepsDt);
    }
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();
  }
  TIME_END("SPHBaseFinalizeDerivs");
}

//------------------------------------------------------------------------------
// Post-state update work
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SPHBase<Dimension>::
postStateUpdate(const typename Dimension::Scalar time,
                const typename Dimension::Scalar dt,
                const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SPHBasePostStateUpdate");

  // Compute the grad h corrrections if needed.  We have to do this in the post-state update
  // because we need boundary conditions applied to position and H first.
  bool result = false;
  if (mGradhCorrection) {
    const auto& WT = this->kernel();
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
    auto        omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
    computeSPHOmegaGradhCorrection(connectivityMap, WT, position, H, omega);
    result = true;
  }

  TIME_END("SPHBasePostStateUpdate");
  return result;
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("SPHBaseGhostBounds");

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> omega = state.fields(HydroFieldNames::omegaGradh, 0.0);

  // FieldList<Dimension, Scalar> volume;
  // const bool updateVolume = (this->densityUpdate() == MassDensityType::VoronoiCellDensity or
  //                            this->densityUpdate() == MassDensityType::SumVoronoiCellDensity);
  // if (updateVolume) {
  //   CHECK(state.fieldNameRegistered(HydroFieldNames::volume));
  //   volume = state.fields(HydroFieldNames::volume, 0.0);
  // }

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(mass);
    boundaryPtr->applyFieldListGhostBoundary(massDensity);
    boundaryPtr->applyFieldListGhostBoundary(specificThermalEnergy);
    boundaryPtr->applyFieldListGhostBoundary(velocity);
    boundaryPtr->applyFieldListGhostBoundary(pressure);
    boundaryPtr->applyFieldListGhostBoundary(soundSpeed);
    boundaryPtr->applyFieldListGhostBoundary(omega);
    // if (updateVolume) boundaryPtr->applyFieldListGhostBoundary(volume);
  }
  TIME_END("SPHBaseGhostBounds");
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("SPHBaseEnforceBounds");

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> omega = state.fields(HydroFieldNames::omegaGradh, 0.0);

  // FieldList<Dimension, Scalar> volume;
  // const bool updateVolume = (this->densityUpdate() == MassDensityType::VoronoiCellDensity or
  //                            this->densityUpdate() == MassDensityType::SumVoronoiCellDensity);
  // if (updateVolume) {
  //   CHECK(state.fieldNameRegistered(HydroFieldNames::volume));
  //   volume = state.fields(HydroFieldNames::volume, 0.0);
  // }

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->enforceFieldListBoundary(mass);
    boundaryPtr->enforceFieldListBoundary(massDensity);
    boundaryPtr->enforceFieldListBoundary(specificThermalEnergy);
    boundaryPtr->enforceFieldListBoundary(velocity);
    boundaryPtr->enforceFieldListBoundary(pressure);
    boundaryPtr->enforceFieldListBoundary(soundSpeed);
    boundaryPtr->enforceFieldListBoundary(omega);
    // if (updateVolume) boundaryPtr->enforceFieldListBoundary(volume);
  }
  TIME_END("SPHBaseEnforceBounds");
}

//------------------------------------------------------------------------------
// Update the volume field in the State.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
updateVolume(State<Dimension>& state,
             const bool boundaries) const {
  TIME_BEGIN("SPHBaseUpdateVol");

  // Pre-conditions.
  REQUIRE(state.fieldNameRegistered(HydroFieldNames::position));
  REQUIRE(state.fieldNameRegistered(HydroFieldNames::volume));
  REQUIRE(state.meshRegistered());

  // Find the global bounding box.
  Vector xmin, xmax;
  const FieldList<Dimension, Vector> positions = state.fields(HydroFieldNames::position, Vector::zero);
  globalBoundingBox<Dimension>(positions, xmin, xmax, 
                               false);     // ghost points

  // Puff things up a bit.
  const Vector delta = 0.1*(xmax - xmin);
  xmin -= delta;
  xmax += delta;

  // Create the mesh.
  Mesh<Dimension>& mesh = state.mesh();
  mesh.clear();
  NodeList<Dimension> voidNodes("void", 0, 0);
  vector<const NodeList<Dimension>*> nodeLists(positions.nodeListPtrs().begin(),
                                               positions.nodeListPtrs().end());
  nodeLists.push_back(&voidNodes);
  generateMesh<Dimension, 
                          typename vector<const NodeList<Dimension>*>::iterator,
                          ConstBoundaryIterator>
    (nodeLists.begin(), nodeLists.end(),
     this->boundaryBegin(),
     this->boundaryEnd(),
     xmin, xmax,
     true,                             // meshGhostNodes
     false,                            // generateVoid
     false,                            // generateParallelConnectivity
     false,                            // removeBoundaryZones
     2.0,                              // voidThreshold
     mesh,
     voidNodes);

  // Extract the volume.
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);

  // Now walk the NodeLists and set the volume.
  const unsigned numNodeLists = volume.size();
  unsigned nodeListi, i, offset, numInternal;
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    offset = mesh.offset(nodeListi);
    numInternal = volume[nodeListi]->numInternalElements();
    for (i = 0; i != numInternal; ++i) {
      volume(nodeListi, i) = mesh.zone(i + offset).volume();
    }
    fill(volume[nodeListi]->begin() + numInternal,
         volume[nodeListi]->end(),
         1.0e-10);
  }

  // Optionally fill in the boundary values for the volume.
  if (boundaries) {
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->applyFieldListGhostBoundary(volume);
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();
  }

  // That's it.
  TIME_END("SPHBaseUpdateVol");
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mVolume, pathName + "/volume");
  // file.write(mEntropy, pathName + "/entropy");
  file.write(mMassDensitySum, pathName + "/massDensitySum");
  file.write(mNormalization, pathName + "/normalization");
  file.write(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.write(mOmegaGradh, pathName + "/omegaGradh");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mGradRho, pathName + "/gradRho");
  file.write(mMassDensityCorrection, pathName + "/massDensityCorrection");
  file.write(mM, pathName + "/M");
  file.write(mLocalM, pathName + "/localM");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mVolume, pathName + "/volume");
  // file.read(mEntropy, pathName + "/entropy");
  file.read(mMassDensitySum, pathName + "/massDensitySum");
  file.read(mNormalization, pathName + "/normalization");
  file.read(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.read(mOmegaGradh, pathName + "/omegaGradh");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mGradRho, pathName + "/gradRho");
  file.read(mMassDensityCorrection, pathName + "/massDensityCorrection");
  file.read(mM, pathName + "/M");
  file.read(mLocalM, pathName + "/localM");

  // // For backwards compatibility on change 3597 -- drop in the near future.
  // for (auto DvDtPtr: mDvDt) DvDtPtr->name(HydroFieldNames::hydroAcceleration);
}

}
