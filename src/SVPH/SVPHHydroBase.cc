//---------------------------------Spheral++----------------------------------//
// SVPHHydroBase -- The SVPH hydrodynamic package for Spheral++.
//
// Created by JMO, Sun Jul 28 20:57:01 PDT 2013
//----------------------------------------------------------------------------//
#include "SVPH/SVPHHydroBase.hh"
#include "SVPH/computeSVPHCorrections.hh"
#include "SVPH/SVPHCorrectionsPolicy.hh"
#include "SVPH/SVPHFieldNames.hh"
#include "SPH/computeSumVoronoiCellMassDensity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosityHandle.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/range.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "FileIO/FileIO.hh"
#include "Mesh/Mesh.hh"

#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

namespace Spheral {


//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHHydroBase<Dimension>::
SVPHHydroBase(const TableKernel<Dimension>& W,
              ArtificialViscosityHandle<Dimension>& Q,
              const double cfl,
              const bool useVelocityMagnitudeForDt,
              const bool compatibleEnergyEvolution,
              const bool XSVPH,
              const bool linearConsistent,
              const MassDensityType densityUpdate,
              const Scalar fcentroidal,
              const Vector& xmin,
              const Vector& xmax):
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mKernel(W),
  mDensityUpdate(densityUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mXSVPH(XSVPH),
  mLinearConsistent(linearConsistent),
  mfcentroidal(fcentroidal),
  mXmin(xmin),
  mXmax(xmax),
  mMeshPtr(MeshPtr(new Mesh<Dimension>())),
  mA(FieldStorageType::Copy),
  mB(FieldStorageType::Copy),
  mGradB(FieldStorageType::Copy),
  mTimeStepMask(FieldStorageType::Copy),
  mPressure(FieldStorageType::Copy),
  mSoundSpeed(FieldStorageType::Copy),
  mVolume(FieldStorageType::Copy),
  mMassDensitySum(FieldStorageType::Copy),
  mXSVPHDeltaV(FieldStorageType::Copy),
  mDxDt(FieldStorageType::Copy),
  mDvDt(FieldStorageType::Copy),
  mDmassDensityDt(FieldStorageType::Copy),
  mDspecificThermalEnergyDt(FieldStorageType::Copy),
  mDvDx(FieldStorageType::Copy),
  mInternalDvDx(FieldStorageType::Copy),
  mPairAccelerations(FieldList<FieldStorageType::Copy),
  mRestart(registerWithRestart(*this)) {
  // Delegate range checking to our assignment methods.
  this->fcentroidal(mfcentroidal);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHHydroBase<Dimension>::
~SVPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Create storage for the pressure and sound speed.
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);

//------------------------------------------------------------------------------
// Second stage of problem start up: initialize values
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Set the moduli.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(HydroFieldNames::soundSpeed, state, derivs);

  // Construct the mesh and volumes.
  NodeList<Dimension> voidNodes("internal void", 0, 0);
  vector<const NodeList<Dimension>*> nodeLists(dataBase.nodeListBegin(), dataBase.nodeListEnd());
  nodeLists.push_back(&voidNodes);
  // std::sort(nodeLists.begin(), nodeLists.end(), typename NodeListRegistrar<Dimension>::NodeListComparator());
  generateMesh<Dimension,
                          typename vector<const NodeList<Dimension>*>::iterator,
                          ConstBoundaryIterator>
    (nodeLists.begin(),
     nodeLists.end(),
     this->boundaryBegin(),
     this->boundaryEnd(),
     mXmin,
     mXmax,
     true,              // mesh ghost nodes
     false,             // generateVoid
     true,              // generateParallelConnectivity
     false,             // removeBoundaryZones
     2.0,               // voidThreshold
     *mMeshPtr,
     voidNodes);
  for (unsigned nodeListi = 0; nodeListi != dataBase.numFluidNodeLists(); ++nodeListi) {
    const unsigned n = mVolume[nodeListi]->numInternalElements();
    const unsigned offset = mMeshPtr->offset(nodeListi);
    for (unsigned i = 0; i != n; ++i) {
      mVolume(nodeListi, i) = mMeshPtr->zone(offset + i).volume();
    }
  }

  // Compute the SVPH normalization and corrections.
  computeSVPHCorrections<Dimension>(dataBase.connectivityMap(),
                                    this->kernel(),
                                    mVolume, dataBase.fluidPosition(), dataBase.fluidHfield(),
                                    mA, mB, mGradB);

}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mA, 0.0, SVPHFieldNames::A_SVPH);
  dataBase.resizeFluidFieldList(mB, Vector::zero(), SVPHFieldNames::B_SVPH);
  dataBase.resizeFluidFieldList(mGradB, Tensor::zero(), SVPHFieldNames::gradB_SVPH);
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);

  // Now register away.
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {

    // Mass.
    state.enroll(fluidNodeListPtr->mass());

    // Mass density.
    if (densityUpdate() == IntegrateDensity) {
      state.enroll(fluidNodeListPtr->massDensity(), make_policy<IncrementBoundedState<Dimension, Scalar>>(fluidNodeListPtr->rhoMin(),
                                                                                                          fluidNodeListPtr->rhoMax()));P
    } else {
      state.enroll(fluidNodeListPtr->massDensity(), make_policy<ReplaceBoundedState<Dimension, Scalar>>(fluidNodeListPtr->rhoMin(),
                                                                                                        fluidNodeListPtr->rhoMax()));
    }

    // Mesh and volume.
    state.enrollMesh(mMeshPtr);
    state.enroll(HydroFieldNames::mesh, make_policy<MeshPolicy<Dimension>>(*this, mXmin, mXmax));
    state.enroll(*mVolume[nodeListi], make_policy<VolumePolicy<Dimension>>());

    // SVPH corrections.
    // All of these corrections are computed in the same method/policy, so we register
    // the A field with the update policy and the others just come along for the ride.
    state.enroll(*mA[nodeListi], make_policy<SVPHCorrectionsPolicy<Dimension>>(dataBase, this->kernel()));
    state.enroll(*mB[nodeListi]);
    state.enroll(*mGradB[nodeListi]);

    // Register the position update.
    state.enroll(fluidNodeListPtr->positions(), make_policy<IncrementState<Dimension, Vector>>());

    // Are we using the compatible energy evolution scheme?
    if (compatibleEnergyEvolution()) {
      state.enroll(fluidNodeListPtr->specificThermalEnergy(), make_policy<NonSymmetricSpecificThermalEnergyPolicy<Dimension>>(dataBase));
      state.enroll(fluidNodeListPtr->velocity(), make_policy<IncrementState<Dimension, Vector>>(HydroFieldNames::position,
                                                                                                HydroFieldNames::specificThermalEnergy));
    } else {
      state.enroll(fluidNodeListPtr->specificThermalEnergy(), make_policy<IncrementState<Dimension, Scalar>>());
      state.enroll(fluidNodeListPtr->velocity(), make_policy<IncrementState<Dimension, Vector>>());
    }

    // Register the time step mask, initialized to 1 so that everything defaults to being
    // checked.
    state.enroll(*mTimeStepMask[nodeListi]);

    // Register the pressure and sound speed.
    state.enroll(*mPressure[nodeListi], make_policy<PressurePolicy<Dimension>>());
    state.enroll(*mSoundSpeed[nodeListi], make_policy<SoundSpeedPolicy<Dimension>>());
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  const auto DxDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const auto DvDtName = HydroFieldNames::hydroAcceleration;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero(), ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMassDensitySum, 0.0, ReplaceState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassFirstMoment, Vector::zero(), HydroFieldNames::massFirstMoment, false);
  dataBase.resizeFluidFieldList(mMassSecondMomentEta, SymTensor::zero(), HydroFieldNames::massSecondMomentEta, false);
  dataBase.resizeFluidFieldList(mMassSecondMomentLab, SymTensor::zero(), HydroFieldNames::massSecondMomentLab, false);
  dataBase.resizeFluidFieldList(mXSVPHDeltaV, Vector::zero(), HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero(), IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero(), HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero(), IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero(), HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero(), HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mPairAccelerations, vector<Vector>(), HydroFieldNames::pairAccelerations, false);

  for (auto [i, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    derivs.enroll(*mHideal[i]);
    derivs.enroll(*mMassDensitySum[i]);
    derivs.enroll(*mWeightedNeighborSum[i]);
    derivs.enroll(*mMassFirstMoment[i]);
    derivs.enroll(*mMassSecondMomentEta[i]);
    derivs.enroll(*mMassSecondMomentLab[i]);
    derivs.enroll(*mXSVPHDeltaV[i]);

    // These two (the position and velocity updates) may be registered
    // by other physics packages as well, so we need to be careful
    // not to duplicate if so.
    const auto DxDtKey = State<Dimension>::buildFieldKey(DxDtName, fluidNodeListPtr->name());
    const auto DvDtKey = State<Dimension>::buildFieldKey(DvDtName, fluidNodeListPtr->name());
    if (not derivs.registered(DxDtKey)) derivs.enroll(*mDxDt[i]);
    if (not derivs.registered(DvDtKey)) derivs.enroll(*mDvDt[i]);

    derivs.enroll(*mDmassDensityDt[i]);
    derivs.enroll(*mDspecificThermalEnergyDt[i]);
    derivs.enroll(*mDHDt[i]);
    derivs.enroll(*mDvDx[i]);
    derivs.enroll(*mInternalDvDx[i]);
    derivs.enroll(*mPairAccelerations[i]);
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Depending on the type of the ArtificialViscosity, dispatch the call to
  // the secondDerivativesLoop
  auto& Qhandle = this->artificialViscosity();
  if (Qhandle.QPiTypeIndex() == std::type_index(typeid(Scalar))) {
      const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Scalar>&>(Qhandle);
      this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  } else {
    CHECK(Qhandle.QPiTypeIndex() == std::type_index(typeid(Tensor)));
    const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Tensor>&>(Qhandle);
    this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  }
}
  
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename QType>
void
SVPHHydroBase<Dimension>::
evaluateDerivativesImpl(const typename Dimension::Scalar time,
                        const typename Dimension::Scalar dt,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivatives,
                        const QType& Q) const {

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero());
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Scalar> A = state.fields(SVPHFieldNames::A_SVPH, 0.0);
  const auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  const auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero());
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero());
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Vector> XSVPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSVPHDeltaV.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Start our big loop over all FluidNodeLists.
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    const int firstGhostNodei = fluidNodeListPtr->firstGhostNode();
    CONTRACT_VAR(firstGhostNodei);

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = fluidNodeListPtr->work();

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const Scalar& Vi = volume(nodeListi, i);
      const Scalar& Ai = A(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Vi > 0.0);
      CHECK(Ai > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSVPHDeltaVi = XSVPHDeltaV(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#if defined __INTEL_COMPILER
#pragma vector always
#endif
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const Vector& rj = position(nodeListj, j);
              const Scalar& mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const Scalar& epsj = specificThermalEnergy(nodeListj, j);
              const Scalar& Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar& Vj = volume(nodeListj, j);
              const Scalar& Aj = A(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Vj > 0.0);
              CHECK(Aj > 0.0);

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Vector& XSVPHDeltaVj = XSVPHDeltaV(nodeListj, j);

              // Node displacement.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // Symmetrized kernel weight and gradient.
              const Vector Hetai = Hi*etai.unitVector();
              Scalar Wi, gWi;
              const std::pair<double, double> WWi = W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
              const Vector gradWi = gWi*Hetai;

              const Vector Hetaj = Hj*etaj.unitVector();
              Scalar Wj, gWj;
              W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
              const Vector gradWj = gWj*Hetaj;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
              }

              // Mass density evolution.
              const Vector vij = vi - vj;
              DrhoDti += Vj*rhoi*vij.dot(gradWj);
              DrhoDtj += Vi*rhoj*vij.dot(gradWi);

              // // Compute the pair-wise artificial viscosity.
              // const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
              //                                           ri, etai, vi, rhoi, ci, Hi,
              //                                           rj, etaj, vj, rhoj, cj, Hj);
              // const Vector Qacci = 0.5*mj*(QPiij.first *gradWi);
              // const Vector Qaccj = 0.5*mi*(QPiij.second*gradWj);
              // const Scalar workQi = vij.dot(Qacci);
              // const Scalar workQj = vij.dot(Qaccj);
              // const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              // const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              // maxViscousPressurei = max(maxViscousPressurei, Qi);
              // maxViscousPressurej = max(maxViscousPressurej, Qj);

              // Compute the pair-wise artificial viscosity.
              QPiType QPiij, QPiji;
              Scalar Qi, Qj;
              Q.QPiij(QPiij, QPiji, Qi, Qj,
                      nodeListi, i, nodeListj, j,
                      ri, Hi, etai, vi, rhoi, ci,  
                      rj, Hj, etaj, vj, rhoj, cj,
                      fClQ, fCqQ, DvDxQ); 
              const Vector Qacci = Ai*Vj*(rhoi*rhoi*QPiij - rhoj*rhoj*QPiji)/rhoi * gradWj;
              const Vector Qaccj = Aj*Vi*(rhoi*rhoi*QPiij - rhoj*rhoj*QPiji)/rhoj * gradWi;
              // const Vector Qacci = -rhoj*QPiij.second*Ai*Vj * gradWj;
              // const Vector Qaccj =  rhoi*QPiij.first *Aj*Vi * gradWi;
              const Scalar workQi = Ai*Vj*rhoi*QPiij*vij.dot(gradWj);
              const Scalar workQj = Aj*Vi*rhoj*QPiji*vij.dot(gradWi);
              // const Scalar workQi = -mi/(mi + mj)*(vi.dot(Qacci) + vj.dot(Qaccj));
              // const Scalar workQj = -mj/(mi + mj)*(vi.dot(Qacci) + vj.dot(Qaccj));
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const Vector deltaDvDti = Ai*Vj*(Pi - Pj)/rhoi*gradWj + Qacci; // - Qacci - Qaccj;
              const Vector deltaDvDtj = Aj*Vi*(Pi - Pj)/rhoj*gradWi + Qaccj; // + Qacci + Qaccj;
              // const Vector aij = (Pi - Pj)*Ai*Vj/rhoi * gradWj + Qacci;
              // const Vector aji = (Pi - Pj)*Aj*Vi/rhoj * gradWi + Qaccj;
              // const Vector Fc = mi*aij + mj*aji;
              // const Vector da = Fc/(mi + mj);
              // const Vector deltaDvDti = aij - da;
              // const Vector deltaDvDtj = aji - da;
              // CHECK2(fuzzyEqual(Fc.dot(mi*deltaDvDti + mj*deltaDvDtj), Fc.magnitude2(), 1.0e-10),
              //        "Pair-wise forces should sum to same central force:  "
              //        << Fc << " "
              //        << (mi*deltaDvDti + mj*deltaDvDtj));
              // CHECK2(fuzzyEqual(-mi*mj*deltaDvDti.dot(deltaDvDtj), mi*mi*deltaDvDti.magnitude2(), 1.0e-10),
              //        "Pair-wise forces should be equal and opposite:  "
              //        << mi*deltaDvDti << " "
              //        << mj*deltaDvDtj);
              DvDti += deltaDvDti;
              DvDtj += deltaDvDtj;

              // Specific thermal energy evolution.
              DepsDti += Ai*Vj*Pi/rhoi*vij.dot(gradWj) + workQi;
              DepsDtj += Aj*Vi*Pj/rhoj*vij.dot(gradWi) + workQj;
              if (mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(deltaDvDti);
                pairAccelerationsj.push_back(deltaDvDtj);
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = Vj*vij.dyad(gradWi);
              const Tensor deltaDvDxj = Vi*vij.dyad(gradWj);
              DvDxi -= deltaDvDxi;
              DvDxj -= deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= deltaDvDxi;
                localDvDxj -= deltaDvDxj;
              }

              // Estimate of delta v (for XSVPH).
              if (mXSVPH and (nodeListi == nodeListj)) {
                const double fXSVPH = max(0.0, min(1.0, abs(vij.dot(rij)*safeInv(vij.magnitude()*rij.magnitude()))));
                CHECK(fXSVPH >= 0.0 and fXSVPH <= 1.0);
                XSVPHDeltaVi -= fXSVPH*Vj*Wj*vij;
                XSVPHDeltaVj += fXSVPH*Vi*Wj*vij;
              }

            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK2(not mCompatibleEnergyEvolution or 
             (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
             (pairAccelerationsi.size() == numNeighborsi),
             "Bad sizing for pair accelerations!  "
             << i << " "
             << firstGhostNodei << " "
             << pairAccelerationsi.size() << " "
             << numNeighborsi);

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Finish the density sum.
      rhoSumi += mi*W0*Hdeti;
      rhoSumi *= Ai;

      // Finish the continuity equation.
      DrhoDti *= Ai;

      // Finish the gradient of the velocity.
      DvDxi *= Ai;
      localDvDxi *= Ai;

      // Determine the position evolution, based on whether we're doing XSVPH or not.
      if (mXSVPH) {
        DxDti = vi + Ai*XSVPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // // Apply any centroidal filtering.
      // DxDti = (1.0 - mfcentroidal)*DxDti + mfcentroidal*(zonei.position() - ri)/dt;

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());

      // Now add the pairwise time for each neighbor we computed here.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          Field<Dimension, Scalar>& workFieldj = nodeLists[nodeListj]->work();
#if defined __INTEL_COMPILER
#pragma vector always
#endif
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              workFieldj(j) += deltaTimePair;
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    FieldList<Dimension, Vector> accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (densityUpdate() == RigorousSumDensity or
      densityUpdate() == SumVoronoiCellDensity) {
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero());
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    computeSumVoronoiCellMassDensity(connectivityMap, this->kernel(), position, mass, volume, H, massDensity);
  } else if (densityUpdate() == SumDensity) {
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    FieldList<Dimension, Scalar> massDensitySum = derivs.fields(ReplaceState<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                                HydroFieldNames::massDensity, 0.0);
    massDensity.assignFields(massDensitySum);
  } else if (densityUpdate() == VoronoiCellDensity) {
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    massDensity = mass / volume;
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> A = state.fields(SVPHFieldNames::A_SVPH, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    (*boundaryItr)->applyFieldListGhostBoundary(A);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> A = state.fields(SVPHFieldNames::A_SVPH, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    (*boundaryItr)->applyFieldListGhostBoundary(A);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mVolume, pathName + "/volume");
  file.write(mMassDensitySum, pathName + "/massDensitySum");
  file.write(mXSVPHDeltaV, pathName + "/XSVPHDeltaV");

  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mVolume, pathName + "/volume");
  file.read(mMassDensitySum, pathName + "/massDensitySum");
  file.read(mXSVPHDeltaV, pathName + "/XSVPHDeltaV");

  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
}

}
