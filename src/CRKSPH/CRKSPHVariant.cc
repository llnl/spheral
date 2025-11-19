//---------------------------------Spheral++----------------------------------//
// CRKSPHVariant -- A development variant of CRKSPH for experimentation.
//
// Created by JMO, Thu Oct 12 14:24:43 PDT 2017
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "CRKSPH/CRKSPHHydroBase.hh"
#include "RK/computeVoronoiVolume.hh"
#include "computeHullVolumes.hh"
#include "computeCRKSPHSumVolume.hh"
#include "computeHVolumes.hh"
#include "editMultimaterialSurfaceTopology.hh"
#include "Utilities/SurfaceNodeCoupling.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/correctSPHSumMassDensity.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeCRKSPHMoments.hh"
#include "detectSurface.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHIntegral.hh"
#include "gradientCRKSPH.hh"
#include "centerOfMass.hh"
#include "volumeSpacing.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "DataBase/updateStateFields.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "ContinuityVolumePolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

#include "CRKSPHVariant.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHVariant<Dimension>::
CRKSPHVariant(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
              ArtificialViscosity<Dimension>& Q,
              const TableKernel<Dimension>& W,
              const TableKernel<Dimension>& WPi,
              const double filter,
              const double cfl,
              const bool useVelocityMagnitudeForDt,
              const bool compatibleEnergyEvolution,
              const bool evolveTotalEnergy,
              const bool XSPH,
              const MassDensityType densityUpdate,
              const HEvolutionType HUpdate,
              const RKOrder correctionOrder,
              const RKVolumeType volumeType,
              const double epsTensile,
              const double nTensile,
              const bool limitMultimaterialTopology):
  CRKSPHHydroBase<Dimension>(smoothingScaleMethod,
                             Q,
                             W,
                             WPi,
                             filter,
                             cfl,
                             useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution,
                             evolveTotalEnergy,
                             XSPH,
                             densityUpdate,
                             HUpdate,
                             correctionOrder,
                             volumeType,
                             epsTensile,
                             nTensile,
                             limitMultimaterialTopology) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHVariant<Dimension>::
~CRKSPHVariant() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVariant<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  cout << "VARIENT INIT PROBLEM START" << endl;

  // Create storage for our internal state.
  this->mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  this->mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  this->mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  this->mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  this->mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  this->mHideal = dataBase.newFluidFieldList(SymTensor::zero(), ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  this->mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  this->mEffViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::effectiveViscousPressure);
  this->mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  this->mViscousWork = dataBase.newFluidFieldList(0.0, HydroFieldNames::viscousWork);
  this->mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  this->mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero(), HydroFieldNames::massSecondMoment);
  this->mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::XSPHDeltaV);
  this->mDxDt = dataBase.newFluidFieldList(Vector::zero(), IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  this->mDvDt = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::hydroAcceleration);
  this->mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  this->mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  this->mDHDt = dataBase.newFluidFieldList(SymTensor::zero(), IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H);
  this->mDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::velocityGradient);
  this->mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::internalVelocityGradient);
  this->mPairAccelerations.clear();
  this->mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero(), "delta centroid");

  this->mA = dataBase.newFluidFieldList(0.0,                        HydroFieldNames::A_CRKSPH);
  this->mB = dataBase.newFluidFieldList(Vector::zero(),               HydroFieldNames::B_CRKSPH);
  this->mGradA = dataBase.newFluidFieldList(Vector::zero(),           HydroFieldNames::gradA_CRKSPH);
  this->mGradB = dataBase.newFluidFieldList(Tensor::zero(),           HydroFieldNames::gradB_CRKSPH);
  this->mM0 = dataBase.newFluidFieldList(0.0,                       HydroFieldNames::m0_CRKSPH);
  this->mM1 = dataBase.newFluidFieldList(Vector::zero(),              HydroFieldNames::m1_CRKSPH);
  this->mM2 = dataBase.newFluidFieldList(SymTensor::zero(),           HydroFieldNames::m2_CRKSPH);
  this->mGradm0 = dataBase.newFluidFieldList(Vector::zero(),          HydroFieldNames::gradM0_CRKSPH);
  this->mGradm1 = dataBase.newFluidFieldList(Tensor::zero(),          HydroFieldNames::gradM1_CRKSPH);
  this->mGradm2 = dataBase.newFluidFieldList(ThirdRankTensor::zero(), HydroFieldNames::gradM2_CRKSPH);
  if (this->mCorrectionOrder == RKOrder::QuadraticOrder) {
    this->mC = dataBase.newFluidFieldList(Tensor::zero(),                HydroFieldNames::C_CRKSPH);
    this->mGradC = dataBase.newFluidFieldList(ThirdRankTensor::zero(),   HydroFieldNames::gradC_CRKSPH);
    this->mM3 = dataBase.newFluidFieldList(ThirdRankTensor::zero(),      HydroFieldNames::m3_CRKSPH);
    this->mM4 = dataBase.newFluidFieldList(FourthRankTensor::zero(),     HydroFieldNames::m4_CRKSPH);
    this->mGradm3 = dataBase.newFluidFieldList(FourthRankTensor::zero(), HydroFieldNames::gradM3_CRKSPH);
    this->mGradm4 = dataBase.newFluidFieldList(FifthRankTensor::zero(),  HydroFieldNames::gradM4_CRKSPH);
  }

  this->mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  this->mEtaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::etaVoidPoints);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVariant<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  cout << "VARIENT INIT PROBLEM START DEPEND" << endl;

  // We need volumes in order to prepare the surface detection.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();

  // Compute the volumes for real.
  if (this->mVolumeType == RKVolumeType::CRKMassOverDensity) {
    this->mVolume.assignFields(mass/massDensity);
  } else if (this->mVolumeType == RKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, mass, H, this->mVolume);
  } else if (this->mVolumeType == RKVolumeType::CRKVoronoiVolume) {
    this->mVolume.assignFields(mass/massDensity);
    FieldList<Dimension, typename Dimension::FacetedVolume> cells;
    FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;
    const FieldList<Dimension, typename Dimension::SymTensor> damage = dataBase.solidEffectiveDamage();
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         vector<typename Dimension::FacetedVolume>(),               // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),      // no holes
                         vector<Boundary<Dimension>*>(),                            // no boundaries
                         FieldList<Dimension, typename Dimension::Scalar>(),        // no weights
                         this->mSurfacePoint, this->mVolume, this->mDeltaCentroid, this->mEtaVoidPoints,    // return values
                         cells, cellFaceFlags);                                     // no return cells
  } else if (this->mVolumeType == RKVolumeType::CRKHullVolume) {
    computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, this->mVolume);
  } else if (this->mVolumeType == RKVolumeType::HVolume) {
    const Scalar nPerh = this->mVolume.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, this->mVolume);
  } else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(this->mVolume);
    if (this->mVolumeType == RKVolumeType::CRKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(this->mVolume);
      (*boundItr)->applyFieldListGhostBoundary(this->mSurfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(this->mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // if (mVolumeType == RKVolumeType::CRKVoronoiVolume) {
  //   // flagSurfaceNeighbors(mSurfacePoint, connectivityMap);
  //   // mVolume = computeShepardsInterpolation(mVolume,
  //   //                                        connectivityMap,
  //   //                                        W,
  //   //                                        position,
  //   //                                        H,
  //   //                                        mVolume);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mVolume);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // }

  // Compute the corrections.
  const NodeCoupling couple;
  computeCRKSPHMoments(connectivityMap, W, this->mVolume, position, H, this->correctionOrder(), couple, this->mM0, this->mM1, this->mM2, this->mM3, this->mM4, this->mGradm0, this->mGradm1, this->mGradm2, this->mGradm3, this->mGradm4);
  computeCRKSPHCorrections(this->mM0, this->mM1, this->mM2, this->mM3, this->mM4, this->mGradm0, this->mGradm1, this->mGradm2, this->mGradm3, this->mGradm4, H, this->surfacePoint(), this->correctionOrder(), this->mA, this->mB, this->mC, this->mGradA, this->mGradB, this->mGradC);

  // This breaks domain independence, so we'll try being inconsistent on the first step.
  // // We need to initialize the velocity gradient if we're using the CRKSPH artificial viscosity.
  // const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  // mDvDx.assignFields(gradientCRKSPH(velocity, position, mVolume, H, mA, mB, mC, mGradA, mGradB, mGradC, connectivityMap, correctionOrder(), W, NodeCoupling()));

  // Initialize the pressure, sound speed, and entropy.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(HydroFieldNames::soundSpeed, state, derivs);
  updateStateFields(HydroFieldNames::entropy, state, derivs);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVariant<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  cout << "VARIENT INIT PROBLEM START" << endl;

  // Create storage for our internal state.
  this->mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  this->mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  this->mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  this->mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  this->mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  this->mHideal = dataBase.newFluidFieldList(SymTensor::zero(), ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  this->mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  this->mEffViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::effectiveViscousPressure);
  this->mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  this->mViscousWork = dataBase.newFluidFieldList(0.0, HydroFieldNames::viscousWork);
  this->mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  this->mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero(), HydroFieldNames::massSecondMoment);
  this->mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::XSPHDeltaV);
  this->mDxDt = dataBase.newFluidFieldList(Vector::zero(), IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  this->mDvDt = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::hydroAcceleration);
  this->mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  this->mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  this->mDHDt = dataBase.newFluidFieldList(SymTensor::zero(), IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H);
  this->mDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::velocityGradient);
  this->mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::internalVelocityGradient);
  this->mPairAccelerations.clear();
  this->mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero(), "delta centroid");

  this->mA = dataBase.newFluidFieldList(0.0,                        HydroFieldNames::A_CRKSPH);
  this->mB = dataBase.newFluidFieldList(Vector::zero(),               HydroFieldNames::B_CRKSPH);
  this->mGradA = dataBase.newFluidFieldList(Vector::zero(),           HydroFieldNames::gradA_CRKSPH);
  this->mGradB = dataBase.newFluidFieldList(Tensor::zero(),           HydroFieldNames::gradB_CRKSPH);
  this->mM0 = dataBase.newFluidFieldList(0.0,                       HydroFieldNames::m0_CRKSPH);
  this->mM1 = dataBase.newFluidFieldList(Vector::zero(),              HydroFieldNames::m1_CRKSPH);
  this->mM2 = dataBase.newFluidFieldList(SymTensor::zero(),           HydroFieldNames::m2_CRKSPH);
  this->mGradm0 = dataBase.newFluidFieldList(Vector::zero(),          HydroFieldNames::gradM0_CRKSPH);
  this->mGradm1 = dataBase.newFluidFieldList(Tensor::zero(),          HydroFieldNames::gradM1_CRKSPH);
  this->mGradm2 = dataBase.newFluidFieldList(ThirdRankTensor::zero(), HydroFieldNames::gradM2_CRKSPH);
  if (this->mCorrectionOrder == RKOrder::QuadraticOrder) {
    this->mC = dataBase.newFluidFieldList(Tensor::zero(),                HydroFieldNames::C_CRKSPH);
    this->mGradC = dataBase.newFluidFieldList(ThirdRankTensor::zero(),   HydroFieldNames::gradC_CRKSPH);
    this->mM3 = dataBase.newFluidFieldList(ThirdRankTensor::zero(),      HydroFieldNames::m3_CRKSPH);
    this->mM4 = dataBase.newFluidFieldList(FourthRankTensor::zero(),     HydroFieldNames::m4_CRKSPH);
    this->mGradm3 = dataBase.newFluidFieldList(FourthRankTensor::zero(), HydroFieldNames::gradM3_CRKSPH);
    this->mGradm4 = dataBase.newFluidFieldList(FifthRankTensor::zero(),  HydroFieldNames::gradM4_CRKSPH);
  }

  // We need volumes in order to prepare the surface detection.
  this->mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  this->mEtaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::etaVoidPoints);
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();

  // Compute the volumes for real.
  if (this->mVolumeType == RKVolumeType::CRKMassOverDensity) {
    this->mVolume.assignFields(mass/massDensity);
  } else if (this->mVolumeType == RKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, mass, H, this->mVolume);
  } else if (this->mVolumeType == RKVolumeType::CRKVoronoiVolume) {
    this->mVolume.assignFields(mass/massDensity);
    FieldList<Dimension, typename Dimension::FacetedVolume> cells;
    FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;
    const FieldList<Dimension, typename Dimension::SymTensor> damage = dataBase.solidEffectiveDamage();
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         vector<typename Dimension::FacetedVolume>(),               // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),      // no holes
                         vector<Boundary<Dimension>*>(),                            // no boundaries
                         FieldList<Dimension, typename Dimension::Scalar>(),        // no weights
                         this->mSurfacePoint, this->mVolume, this->mDeltaCentroid, this->mEtaVoidPoints,    // return values
                         cells, cellFaceFlags);                                     // no return cells
  } else if (this->mVolumeType == RKVolumeType::CRKHullVolume) {
    computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, this->mVolume);
  } else if (this->mVolumeType == RKVolumeType::HVolume) {
    const Scalar nPerh = this->mVolume.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, this->mVolume);
  } else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(this->mVolume);
    if (this->mVolumeType == RKVolumeType::CRKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(this->mVolume);
      (*boundItr)->applyFieldListGhostBoundary(this->mSurfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(this->mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // if (mVolumeType == RKVolumeType::CRKVoronoiVolume) {
  //   // flagSurfaceNeighbors(mSurfacePoint, connectivityMap);
  //   // mVolume = computeShepardsInterpolation(mVolume,
  //   //                                        connectivityMap,
  //   //                                        W,
  //   //                                        position,
  //   //                                        H,
  //   //                                        mVolume);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mVolume);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // }

  // Compute the corrections.
  const NodeCoupling couple;
  computeCRKSPHMoments(connectivityMap, W, this->mVolume, position, H, this->correctionOrder(), couple, this->mM0, this->mM1, this->mM2, this->mM3, this->mM4, this->mGradm0, this->mGradm1, this->mGradm2, this->mGradm3, this->mGradm4);
  computeCRKSPHCorrections(this->mM0, this->mM1, this->mM2, this->mM3, this->mM4, this->mGradm0, this->mGradm1, this->mGradm2, this->mGradm3, this->mGradm4, H, this->surfacePoint(), this->correctionOrder(), this->mA, this->mB, this->mC, this->mGradA, this->mGradB, this->mGradC);

  // This breaks domain independence, so we'll try being inconsistent on the first step.
  // // We need to initialize the velocity gradient if we're using the CRKSPH artificial viscosity.
  // const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  // mDvDx.assignFields(gradientCRKSPH(velocity, position, mVolume, H, mA, mB, mC, mGradA, mGradB, mGradC, connectivityMap, correctionOrder(), W, NodeCoupling()));

  // Initialize the pressure, sound speed, and entropy.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(HydroFieldNames::soundSpeed, state, derivs);
  updateStateFields(HydroFieldNames::entropy, state, derivs);
}

//------------------------------------------------------------------------------
// Initialize the hydro before evaluating derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CRKSPHVariant<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // Compute the kernel correction fields.
  cout << "VARIENT INITIALIZE" << endl;
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero());
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const FieldList<Dimension, int> surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero());
  FieldList<Dimension, Tensor> C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero());
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero());
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero());
  FieldList<Dimension, ThirdRankTensor> gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero());
  FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CRKSPH, 0.0);
  FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CRKSPH, Vector::zero());
  FieldList<Dimension, SymTensor> m2 = state.fields(HydroFieldNames::m2_CRKSPH, SymTensor::zero());
  FieldList<Dimension, ThirdRankTensor> m3 = state.fields(HydroFieldNames::m3_CRKSPH, ThirdRankTensor::zero());
  FieldList<Dimension, FourthRankTensor> m4 = state.fields(HydroFieldNames::m4_CRKSPH, FourthRankTensor::zero());
  FieldList<Dimension, Vector> gradm0 = state.fields(HydroFieldNames::gradM0_CRKSPH, Vector::zero());
  FieldList<Dimension, Tensor> gradm1 = state.fields(HydroFieldNames::gradM1_CRKSPH, Tensor::zero());
  FieldList<Dimension, ThirdRankTensor> gradm2 = state.fields(HydroFieldNames::gradM2_CRKSPH, ThirdRankTensor::zero());
  FieldList<Dimension, FourthRankTensor> gradm3 = state.fields(HydroFieldNames::gradM3_CRKSPH, FourthRankTensor::zero());
  FieldList<Dimension, FifthRankTensor> gradm4 = state.fields(HydroFieldNames::gradM4_CRKSPH, FifthRankTensor::zero());

  // Change CRKSPH weights here if need be!
  const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  const NodeCoupling couple;
  computeCRKSPHMoments(connectivityMap, W, vol, position, H, this->correctionOrder(), couple, m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, surfacePoint, this->correctionOrder(), A, B, C, gradA, gradB, gradC);

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(A);
    (*boundItr)->applyFieldListGhostBoundary(B);
    (*boundItr)->applyFieldListGhostBoundary(C);
    (*boundItr)->applyFieldListGhostBoundary(gradA);
    (*boundItr)->applyFieldListGhostBoundary(gradB);
    (*boundItr)->applyFieldListGhostBoundary(gradC);
  }

  // Get the artificial viscosity and initialize it.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();
  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               W);
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVariant<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();
  cout << "VARIENT EVALUATE" << endl;

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  order = this->correctionOrder();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero(), true);
  const auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero(), true);
  const auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero());
  const auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero(), true);
  const auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero(), true);
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);
  CHECK(B.size() == numNodeLists or order == RKOrder::ZerothOrder);
  CHECK(C.size() == numNodeLists or order != RKOrder::QuadraticOrder);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists or order == RKOrder::ZerothOrder);
  CHECK(gradC.size() == numNodeLists or order != RKOrder::QuadraticOrder);
  CHECK(surfacePoint.size() == numNodeLists);

  // Derivative FieldLists.
  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero());
  auto DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto DepsDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto DHDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero());
  auto Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero());
  auto maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero());
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (this->mCompatibleEnergyEvolution) {
    auto nodeListi = 0;
    for (auto itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        const size_t n = connectivityMap.numNeighborsForNode(*itr, i);
        pairAccelerations(nodeListi, i).reserve(n);
      }
    }
  }

  // Some scratch variables.
  Scalar Ai, Aj;
  Vector gradAi, gradAj, forceij, forceji;
  Vector Bi = Vector::zero(), Bj = Vector::zero();
  Tensor Ci = Tensor::zero(), Cj = Tensor::zero();
  Tensor gradBi = Tensor::zero(), gradBj = Tensor::zero();
  ThirdRankTensor gradCi = ThirdRankTensor::zero(), gradCj = ThirdRankTensor::zero();
  Scalar gWi, gWj, Wi, Wj, gW0i, gW0j, W0i, W0j;
  Vector gradWi, gradWj, gradW0i, gradW0j;
  Vector deltagrad;

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const auto& nodeList = **itr;
    const auto firstGhostNodei = nodeList.firstGhostNode();
    const auto hmin = nodeList.hmin();
    const auto hmax = nodeList.hmax();
    const auto hminratio = nodeList.hminratio();
    const auto maxNumNeighbors = nodeList.maxNumNeighbors();
    const auto nPerh = nodeList.nodesPerSmoothingScale();

    // Get the work field for this NodeList.
    auto& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto rhoi = massDensity(nodeListi, i);
      const auto epsi = specificThermalEnergy(nodeListi, i);
      const auto Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto ci = soundSpeed(nodeListi, i);
      Ai = A(nodeListi, i);
      gradAi = gradA(nodeListi, i);
      if (order != RKOrder::ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
      }
      if (order == RKOrder::QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
      }
      const auto Hdeti = Hi.Determinant();
      const auto weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);
      // CHECK2(Ai > 0.0, i << " " << Ai);
      CHECK2(Hdeti > 0.0, i << " " << Hdeti);
      CHECK2(weighti > 0.0, i << " " << weighti);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& viscousWorki = viscousWork(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#if defined __INTEL_COMPILER
#pragma vector always
#endif
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const auto& rj = position(nodeListj, j);
              const auto  mj = mass(nodeListj, j);
              const auto& vj = velocity(nodeListj, j);
              const auto  rhoj = massDensity(nodeListj, j);
              const auto  epsj = specificThermalEnergy(nodeListj, j);
              const auto  Pj = pressure(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              const auto  cj = soundSpeed(nodeListj, j);
              Aj = A(nodeListj, j);
              gradAj = gradA(nodeListj, j);
              if (order != RKOrder::ZerothOrder) {
                Bj = B(nodeListj, j);
                gradBj = gradB(nodeListj, j);
              }
              if (order == RKOrder::QuadraticOrder) {
                Cj = C(nodeListj, j);
                gradCj = gradC(nodeListj, j);
              }
              const auto Hdetj = Hj.Determinant();
              const auto weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Aj > 0.0 or j >= firstGhostNodej);
              CHECK(Hdetj > 0.0);
              CHECK(weightj > 0.0);

              auto& DxDtj = DxDt(nodeListj, j);
              auto& DrhoDtj = DrhoDt(nodeListj, j);
              auto& DvDtj = DvDt(nodeListj, j);
              auto& DepsDtj = DepsDt(nodeListj, j);
              auto& DvDxj = DvDx(nodeListj, j);
              auto& localDvDxj = localDvDx(nodeListj, j);
              auto& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              auto& effViscousPressurej = effViscousPressure(nodeListj, j);
              auto& viscousWorkj = viscousWork(nodeListj, j);
              auto& pairAccelerationsj = pairAccelerations(nodeListj, j);
              auto& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              auto& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              auto& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Find the effective weights of i->j and j->i.
              // const auto wi = 2.0*weighti*weightj/(weighti + weightj);
              const auto wij = 0.5*(weighti + weightj);

              // Node displacement.
              const auto rij = ri - rj;
              const auto etai = Hi*rij;
              const auto etaj = Hj*rij;
              const auto etaMagi = etai.magnitude();
              const auto etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);
              const auto vij = vi - vj;

              // Symmetrized kernel weight and gradient.
              CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, order,  rij,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi);
              CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, order, -rij, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj);
              deltagrad = gradWj - gradWi;
              const auto gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);
              const auto gradWSPHj = (Hj*etaj.unitVector())*W.gradValue(etaj.magnitude(), Hdetj);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const auto fweightij = nodeListi == nodeListj ? 1.0 : mj*rhoi/(mi*rhoj);
              const auto rij2 = rij.magnitude2();
              const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
              weightedNeighborSumi +=     fweightij*std::abs(gWi);
              weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
              massSecondMomenti +=     fweightij*gradWSPHi.magnitude2()*thpt;
              massSecondMomentj += 1.0/fweightij*gradWSPHj.magnitude2()*thpt;

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const auto QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                        ri, etai, vi, rhoi, ci, Hi,
                                        rj, etaj, vj, rhoj, cj, Hj);
              const auto Qaccij = (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second).dot(deltagrad);
              // const auto workQij = 0.5*(vij.dot(Qaccij));
              const auto workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);                // CRK
              const auto workQj = rhoi*rhoi*QPiij.first .dot(vij).dot(deltagrad);                // CRK
              // const auto workQVi =  vij.dot((rhoj*rhoj*QPiij.second).dot(gradWj));               //RK V and RK I Work
              // const auto workQVj =  vij.dot((rhoi*rhoi*QPiij.first).dot(gradWi));                //RK V and RK I Work
              const auto Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const auto Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                                 // We need tighter timestep controls on the Q with CRK
              maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
              effViscousPressurei += wij * Qi * Wj;
              effViscousPressurej += wij * Qj * Wi;
              viscousWorki += 0.5*wij*wij/mi*workQi;
              viscousWorkj += 0.5*wij*wij/mj*workQj;

              // Velocity gradient.
              DvDxi -= wij*vij.dyad(gradWj);
              DvDxj += wij*vij.dyad(gradWi);
              if (nodeListi == nodeListj) {
                localDvDxi -= wij*vij.dyad(gradWj);
                localDvDxj += wij*vij.dyad(gradWi);
              }

              // // Mass density gradient.
              // gradRhoi += wij*(rhoj - rhoi)*gradWj;
              // gradRhoj += wij*(rhoi - rhoj)*gradWi;

              // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
              // Momentum
              forceij = 0.5*wij*wij*((Pi + Pj)*deltagrad + Qaccij);                    // Type III CRK interpoint force.
              forceji = 0.5*wij*wij*((Pi + Pj)*deltagrad + Qaccij);                    // Type III CRK interpoint force.
              DvDti -= forceij/mi;
              DvDtj += forceji/mj; 
              if (this->mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(-forceij/mi);
                pairAccelerationsj.push_back( forceji/mj);
              }

              // Energy
              DepsDti += 0.5*wij*wij*(Pj*vij.dot(deltagrad) + workQi)/mi;              // CRK
              DepsDtj += 0.5*wij*wij*(Pi*vij.dot(deltagrad) + workQj)/mj;              // CRK

              // Estimate of delta v (for XSPH).
              if (this->mXSPH and (nodeListi == nodeListj)) {
                XSPHDeltaVi -= wij*Wj*vij;
                XSPHDeltaVj += wij*Wi*vij;
              }
                
            }
          }
        }
      }
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not this->mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // // For a surface point, add the RK thermal energy evolution.
      // // DepsDti -= Pi/rhoi*DvDxi.Trace();
      // if (surfacePoint(nodeListi, i) > 1) DepsDti -= Pi/rhoi*DvDxi.Trace();

      // Time evolution of the mass density.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (this->mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (this->mXSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = this->mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                                   ri,
                                                                   DvDxi,
                                                                   hmin,
                                                                   hmax,
                                                                   hminratio,
                                                                   nPerh);
      Hideali = this->mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                              ri,
                                                              weightedNeighborSumi,
                                                              massSecondMomenti,
                                                              W,
                                                              hmin,
                                                              hmax,
                                                              hminratio,
                                                              nPerh,
                                                              connectivityMap,
                                                              nodeListi,
                                                              i);
    }
  }
}

}

