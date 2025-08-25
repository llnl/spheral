//---------------------------------Spheral++----------------------------------//
// SolidFSISPH -- SolidSPHHydro modified to better handle 
//                multimaterial material problems with interfaces
//                and large density ratios. 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "Physics/GenericHydro.hh"

#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/SolidEquationOfState.hh" 

#include "Hydro/computeSPHVolume.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"

#include "Strength/SolidFieldNames.hh"
#include "Strength/DeviatoricStressPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/PlasticStrainPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/PureReplaceState.hh"
#include "DataBase/updateStateFields.hh"
#include "DataBase/ReplaceWithRatioPolicy.hh"

#include "ArtificialViscosity/ArtificialViscosityView.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/Timer.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/range.hh"

#include "FSISPH/SolidFSISPH.hh"
#include "FSISPH/FSIFieldNames.hh"
#include "FSISPH/computeFSISPHSumMassDensity.hh"
#include "FSISPH/computeHWeightedFSISPHSumMassDensity.hh"
#include "FSISPH/computeInterfacePressureCorrectedSumMassDensity.hh"
#include "FSISPH/SlideSurface.hh"

#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

using std::make_shared;
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


inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor::zero;
  }
}

inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

inline
Dim<3>::SymTensor
tensileStressCorrection(const Dim<3>::SymTensor& sigma) {
  const EigenStruct<3> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  const double lambdaz = eigen.eigenValues.z();
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,                              0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0,                              0.0,                              (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}


//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidFSISPH<Dimension>::
SolidFSISPH(DataBase<Dimension>& dataBase,
            ArtificialViscosity<Dimension>& Q,
            SlideSurface<Dimension>& slides,
            const TableKernel<Dimension>& W,
            const double cfl,
            const double surfaceForceCoefficient,
            const double densityStabilizationCoefficient,
            const double specificThermalEnergyDiffusionCoefficient,
            const double xsphCoefficient,
            const InterfaceMethod interfaceMethod,
            const KernelAveragingMethod kernelAveragingMethod,
            const std::vector<int> sumDensityNodeLists,
            const bool useVelocityMagnitudeForDt,
            const bool compatibleEnergyEvolution,
            const bool evolveTotalEnergy,
            const bool linearCorrectGradients,
            const bool decoupleDamagedMaterial,
            const double interfacePmin,
            const double interfaceNeighborAngleThreshold,
            const FSIMassDensityMethod densityUpdate,
            const double epsTensile,
            const double nTensile,
            const Vector& xmin,
            const Vector& xmax):
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mKernel(W),
  mSlideSurface(slides),
  mDensityUpdate(densityUpdate),
  mInterfaceMethod(interfaceMethod),
  mKernelAveragingMethod(kernelAveragingMethod),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mLinearCorrectGradients(linearCorrectGradients),
  mDecoupleDamagedMaterial(decoupleDamagedMaterial),
  mApplySelectDensitySum(std::any_of(sumDensityNodeLists.begin(), sumDensityNodeLists.end(), [](int x) { return x == 1; })),
  mSumDensityNodeLists(sumDensityNodeLists),
  mSurfaceForceCoefficient(surfaceForceCoefficient),
  mDensityStabilizationCoefficient(densityStabilizationCoefficient),
  mSpecificThermalEnergyDiffusionCoefficient(specificThermalEnergyDiffusionCoefficient),
  mXSPHCoefficient(xsphCoefficient),
  mInterfacePmin(interfacePmin),
  mInterfaceNeighborAngleThreshold(interfaceNeighborAngleThreshold),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mxmin(xmin),
  mxmax(xmax),
  mPairAccelerationsPtr(),
  mPairDepsDtPtr(),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mDamagedPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mBulkModulus(FieldStorageType::CopyFields),
  mShearModulus(FieldStorageType::CopyFields),
  mYieldStrength(FieldStorageType::CopyFields),
  mPlasticStrain0(FieldStorageType::CopyFields),
  //mInverseEquivalentDeviatoricStress(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDdeviatoricStressDt(FieldStorageType::CopyFields),
  mDPDx(FieldStorageType::CopyFields),
  mDepsDx(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mLocalM(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mInterfaceFlags(FieldStorageType::CopyFields),
  mInterfaceAreaVectors(FieldStorageType::CopyFields),
  mInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceSmoothness(FieldStorageType::CopyFields),
  mNewInterfaceFlags(FieldStorageType::CopyFields),
  mNewInterfaceAreaVectors(FieldStorageType::CopyFields),
  mNewInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceSmoothnessNormalization(FieldStorageType::CopyFields),
  mInterfaceFraction(FieldStorageType::CopyFields),
  mNewInterfaceSmoothness(FieldStorageType::CopyFields),
  mInterfaceAngles(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mDamagedPressure = dataBase.newFluidFieldList(0.0, SolidFieldNames::damagedPressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
  mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
  //mInverseEquivalentDeviatoricStress = dataBase.newFluidFieldList(0.0, FSIFieldNames::inverseEquivalentDeviatoricStress);
  mVolume = dataBase.newFluidFieldList(0.0,HydroFieldNames::volume);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
  mDPDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::pressureGradient);
  mDepsDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::specificThermalEnergyGradient);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
  mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection);
  mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
  mInterfaceFlags = dataBase.newFluidFieldList(int(0),  FSIFieldNames::interfaceFlags);
  mInterfaceAreaVectors = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceAreaVectors);
  mInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceNormals);
  mInterfaceSmoothness = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceSmoothness);
  mNewInterfaceFlags = dataBase.newFluidFieldList(int(0), PureReplaceState<Dimension,int>::prefix() + FSIFieldNames::interfaceFlags);
  mNewInterfaceAreaVectors = dataBase.newFluidFieldList(Vector::one, PureReplaceState<Dimension,Vector>::prefix() + FSIFieldNames::interfaceAreaVectors);
  mNewInterfaceNormals = dataBase.newFluidFieldList(Vector::one, PureReplaceState<Dimension,Vector>::prefix() + FSIFieldNames::interfaceNormals);
  mInterfaceSmoothnessNormalization = dataBase.newFluidFieldList(0.0, FSIFieldNames::interfaceSmoothnessNormalization);
  mInterfaceFraction = dataBase.newFluidFieldList(0.0, FSIFieldNames::interfaceFraction);
  mNewInterfaceSmoothness = dataBase.newFluidFieldList(0.0, PureReplaceState<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness);
  mInterfaceAngles = dataBase.newFluidFieldList(0.0, FSIFieldNames::interfaceAngles);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Set the moduli.
  updateStateFields(HydroFieldNames::pressure, state, derivs);
  updateStateFields(HydroFieldNames::soundSpeed, state, derivs);
  updateStateFields(SolidFieldNames::bulkModulus, state, derivs);
  updateStateFields(SolidFieldNames::shearModulus, state, derivs);
  updateStateFields(SolidFieldNames::yieldStrength, state, derivs);

  const auto& mass = dataBase.fluidMass();
  const auto& massDensity = dataBase.fluidMassDensity();
  computeSPHVolume(mass,massDensity,mVolume);
}

//------------------------------------------------------------------------------
// Register states
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("SolidFSISPHregisterState");

  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "FSISPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  FieldList<Dimension, Vector> position = dataBase.solidPosition();
  FieldList<Dimension, Scalar> mass = dataBase.solidMass();
  FieldList<Dimension, Scalar> massDensity = dataBase.solidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.solidHfield();
  FieldList<Dimension, Vector> velocity = dataBase.solidVelocity();
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.solidSpecificThermalEnergy();
  FieldList<Dimension, SymTensor> deviatoricStress = dataBase.solidDeviatoricStress();
  FieldList<Dimension, Scalar> plasticStrain = dataBase.solidPlasticStrain();
  FieldList<Dimension, SymTensor> damage = dataBase.solidDamage();
  FieldList<Dimension, int> fragIDs = dataBase.solidFragmentIDs();
  FieldList<Dimension, int> pTypes = dataBase.solidParticleTypes();
  
  for (auto [nodeListi, nodeListPtr]: enumerate(dataBase.solidNodeListBegin(), dataBase.solidNodeListEnd())) {
    *mPlasticStrain0[nodeListi] = nodeListPtr->plasticStrain();
    (*mPlasticStrain0[nodeListi]).name(SolidFieldNames::plasticStrain + "0");
  }

  // Create the local storage.
  dataBase.resizeSolidFieldList(mTimeStepMask, int(1), HydroFieldNames::timeStepMask);
  dataBase.resizeSolidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  dataBase.resizeSolidFieldList(mPressure, 0.0, HydroFieldNames::pressure, false);
  dataBase.resizeSolidFieldList(mSoundSpeed, 0.0, HydroFieldNames::soundSpeed, false);
  dataBase.resizeSolidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeSolidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeSolidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeSolidFieldList(mPlasticStrain0, 0.0, SolidFieldNames::plasticStrain + "0", false);
  dataBase.resizeSolidFieldList(mDamagedPressure, 0.0, SolidFieldNames::damagedPressure, false);
  dataBase.resizeSolidFieldList(mInterfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals, false);
  dataBase.resizeSolidFieldList(mInterfaceFraction, 0.0, FSIFieldNames::interfaceFraction, false); 
  dataBase.resizeSolidFieldList(mInterfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness, false);
  //dataBase.resizeSolidFieldList(mInverseEquivalentDeviatoricStress, 0.0, FSIFieldNames::inverseEquivalentDeviatoricStress, false);
  dataBase.resizeFluidFieldList(mInterfaceFlags, int(0), FSIFieldNames::interfaceFlags,false);
  dataBase.resizeFluidFieldList(mInterfaceAreaVectors, Vector::zero, FSIFieldNames::interfaceAreaVectors,false);
  dataBase.resizeFluidFieldList(mInterfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mInterfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness,false);

  // Register the deviatoric stress and plastic strain to be evolved.
  auto positionPolicy = make_policy<IncrementState<Dimension, Vector>>();
  auto velocityPolicy = make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::position,HydroFieldNames::specificThermalEnergy},true);
  auto Ppolicy = make_policy<PressurePolicy<Dimension>>();
  auto csPolicy = make_policy<SoundSpeedPolicy<Dimension>>();
  auto plasticStrainPolicy = make_policy<PlasticStrainPolicy<Dimension>>();
  auto bulkModulusPolicy = make_policy<BulkModulusPolicy<Dimension>>();
  auto shearModulusPolicy = make_policy<ShearModulusPolicy<Dimension>>();
  auto yieldStrengthPolicy = make_policy<YieldStrengthPolicy<Dimension>>();
  auto volumePolicy = make_policy<ReplaceWithRatioPolicy<Dimension,Scalar>>({HydroFieldNames::massDensity},
                                                                            HydroFieldNames::mass, HydroFieldNames::massDensity);
  auto interfaceFlagsPolicy = make_policy<PureReplaceState<Dimension,int>>();
  auto interfaceAreaVectorsPolicy = make_policy<PureReplaceState<Dimension,Vector>>();
  auto interfaceNormalsPolicy = make_policy<PureReplaceState<Dimension,Vector>>();
  auto interfaceSmoothnessPolicy = make_policy<PureReplaceState<Dimension,Scalar>>();
  
  auto deviatoricStressPolicy = make_policy<DeviatoricStressPolicy<Dimension>>();
  state.enroll(deviatoricStress, deviatoricStressPolicy);

  if(this->compatibleEnergyEvolution()){
    auto  thermalEnergyPolicy = make_policy<CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>>(dataBase);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }else if (mEvolveTotalEnergy) {
    auto  thermalEnergyPolicy = make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>();
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  } else {
    auto  thermalEnergyPolicy = make_policy<IncrementState<Dimension, Scalar>>();
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }

  for (auto [nodeListi, nodeListPtr]: enumerate(dataBase.solidNodeListBegin(), dataBase.solidNodeListEnd())) {
    state.enroll(*massDensity[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>(nodeListPtr->rhoMin(),
                                                                                                nodeListPtr->rhoMax()));
  }

  state.enroll(position,             positionPolicy);
  state.enroll(velocity,             velocityPolicy);
  state.enroll(mVolume,              volumePolicy);
  state.enroll(mSoundSpeed,          csPolicy);
  state.enroll(mPressure,            Ppolicy);
  state.enroll(mDamagedPressure);                              // Updated by PressurePolicy
  state.enroll(mBulkModulus,         bulkModulusPolicy);
  state.enroll(mShearModulus,        shearModulusPolicy);
  state.enroll(mYieldStrength,       yieldStrengthPolicy);
  state.enroll(plasticStrain,        plasticStrainPolicy);

  state.enroll(mInterfaceFlags,       interfaceFlagsPolicy);
  state.enroll(mInterfaceAreaVectors, interfaceAreaVectorsPolicy); 
  state.enroll(mInterfaceNormals,     interfaceNormalsPolicy); 
  state.enroll(mInterfaceSmoothness,  interfaceSmoothnessPolicy); 

  state.enroll(mTimeStepMask);
  state.enroll(mass);
  //state.enroll(mInverseEquivalentDeviatoricStress);
  state.enroll(damage);
  state.enroll(fragIDs);
  state.enroll(pTypes);
  state.enroll(mPlasticStrain0);

  TIME_END("SolidFSISPHregisterState");
}

//------------------------------------------------------------------------------
// Register Derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
registerDerivatives(DataBase<Dimension>&  dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidFSISPHregisterDerivs");

  FieldList<Dimension, Scalar> plasticStrainRate = dataBase.solidPlasticStrainRate();

  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress, false);
  dataBase.resizeFluidFieldList(mDPDx, Vector::zero, FSIFieldNames::pressureGradient, false);
  dataBase.resizeFluidFieldList(mDepsDx, Vector::zero, FSIFieldNames::specificThermalEnergyGradient, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero, HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mNormalization, 0.0, HydroFieldNames::normalization, false);
  dataBase.resizeFluidFieldList(mNewInterfaceFlags, int(0),  PureReplaceState<Dimension,int>::prefix() + FSIFieldNames::interfaceFlags,false);
  dataBase.resizeFluidFieldList(mNewInterfaceAreaVectors, Vector::zero,  PureReplaceState<Dimension,Vector>::prefix() + FSIFieldNames::interfaceAreaVectors,false);
  dataBase.resizeFluidFieldList(mNewInterfaceNormals, Vector::zero,  PureReplaceState<Dimension,Vector>::prefix() + FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mInterfaceSmoothnessNormalization, 0.0, FSIFieldNames::interfaceSmoothnessNormalization,false); 
  dataBase.resizeFluidFieldList(mInterfaceFraction, 0.0, FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mNewInterfaceSmoothness, 0.0,  PureReplaceState<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness,false);
  dataBase.resizeFluidFieldList(mInterfaceAngles, 0.0,  FSIFieldNames::interfaceAngles,false);

  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }

  CHECK(not derivs.registered(mDvDt));

  if (mCompatibleEnergyEvolution) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    mPairDepsDtPtr = std::make_unique<PairWorkType>(connectivityMap);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
    derivs.enroll(HydroFieldNames::pairWork, *mPairDepsDtPtr);
  }

  derivs.enroll(plasticStrainRate);
  derivs.enroll(mXSPHDeltaV);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mDvDt);
  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDdeviatoricStressDt);
  derivs.enroll(mDPDx);
  derivs.enroll(mDepsDx);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enroll(mM);
  derivs.enroll(mLocalM);
  derivs.enroll(mNormalization);
  derivs.enroll(mNewInterfaceFlags);
  derivs.enroll(mNewInterfaceAreaVectors);
  derivs.enroll(mNewInterfaceNormals);
  derivs.enroll(mInterfaceSmoothnessNormalization);
  derivs.enroll(mInterfaceFraction);
  derivs.enroll(mNewInterfaceSmoothness);
  derivs.enroll(mInterfaceAngles);

  TIME_END("SolidFSISPHregisterDerivs");
}

//------------------------------------------------------------------------------
// FSI specialized density summmation
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("SolidFSISPHpreStepInitialize");
  if (mApplySelectDensitySum){
  switch(this->densityUpdate()){
   case FSIMassDensityMethod::FSISumMassDensity:
    {
      const auto& W = this->kernel();
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
            auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeFSISPHSumMassDensity(connectivityMap, W, mSumDensityNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
     break;
    }
    case FSIMassDensityMethod::PressureCorrectSumMassDensity:
    {
      const auto& W = this->kernel();
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
            auto  pressure = state.fields(HydroFieldNames::pressure, 0.0);
            auto  soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
            auto  volume = state.fields(HydroFieldNames::volume, 0.0);
            auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeInterfacePressureCorrectedSumMassDensity(connectivityMap, 
                                                      W, 
                                                      mSumDensityNodeLists, 
                                                      position, 
                                                      mass, 
                                                      H,
                                                      volume,
                                                      pressure,
                                                      soundSpeed,
                                                      massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      break;
    }
    case FSIMassDensityMethod::HWeightedSumMassDensity:
    {
      const auto& W = this->kernel();
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
            auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeHWeightedFSISPHSumMassDensity(connectivityMap, W, mSumDensityNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      break;
    }
    default:
      break;
  }
  }
  TIME_END("SolidFSISPHpreStepInitialize");
}

//------------------------------------------------------------------------------
// For compatible energy we need to apply the bc conditions to acceleration
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
finalizeDerivatives(const Scalar /*time*/, 
                    const Scalar  /*dt*/,
                    const DataBase<Dimension>&  /*dataBase*/, 
                    const State<Dimension>& /*state*/,
                          StateDerivatives<Dimension>&  derivs) const {                 
  if (this->compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
           (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    }
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }

} // finalize


//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> damagedPressure = state.fields(SolidFieldNames::damagedPressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  //FieldList<Dimension, Scalar> invSqrtJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);
  FieldList<Dimension, int> interfaceFlags = state.fields(FSIFieldNames::interfaceFlags, int(0));
  FieldList<Dimension, Vector> interfaceAreaVectors = state.fields(FSIFieldNames::interfaceAreaVectors, Vector::zero);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(damagedPressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
    (*boundaryItr)->applyFieldListGhostBoundary(fragIDs);
    (*boundaryItr)->applyFieldListGhostBoundary(pTypes);
    //(*boundaryItr)->applyFieldListGhostBoundary(invSqrtJ2);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceFlags);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceAreaVectors);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceNormals);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceSmoothness);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> damagedPressure = state.fields(SolidFieldNames::damagedPressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  //FieldList<Dimension, Scalar> invSqrtJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);
  FieldList<Dimension, int> interfaceFlags = state.fields(FSIFieldNames::interfaceFlags, int(0));
  FieldList<Dimension, Vector> interfaceAreaVectors = state.fields(FSIFieldNames::interfaceAreaVectors, Vector::zero);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(damagedPressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
    (*boundaryItr)->enforceFieldListBoundary(fragIDs);
    (*boundaryItr)->enforceFieldListBoundary(pTypes);
    //(*boundaryItr)->enforceFieldListBoundary(invSqrtJ2);
    (*boundaryItr)->enforceFieldListBoundary(interfaceFlags);
    (*boundaryItr)->enforceFieldListBoundary(interfaceAreaVectors);
    (*boundaryItr)->enforceFieldListBoundary(interfaceNormals);
    (*boundaryItr)->enforceFieldListBoundary(interfaceSmoothness);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mDamagedPressure, pathName + "/damagedPressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mVolume, pathName + "/volume");
  file.write(mBulkModulus, pathName + "/bulkModulus");
  file.write(mShearModulus, pathName + "/shearModulus");
  file.write(mYieldStrength, pathName + "/yieldStrength");
  file.write(mPlasticStrain0, pathName + "/plasticStrain0");
  //file.write(mInverseEquivalentDeviatoricStress, pathName + "/inverseEquivalentDeviatoricStress");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.write(mDPDx, pathName + "/DpDx");
  file.write(mDepsDx, pathName + "/DepsDx");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mM, pathName + "/M");
  file.write(mLocalM, pathName + "/localM");
  file.write(mNormalization, pathName + "/normalization");
  file.write(mInterfaceFlags, pathName + "/interfaceFlags");
  file.write(mInterfaceAreaVectors, pathName + "/interfaceAreaVectors");
  file.write(mInterfaceNormals, pathName + "/interfaceNormals");
  file.write(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.write(mNewInterfaceFlags, pathName + "/newInterfaceFlags");
  file.write(mNewInterfaceAreaVectors, pathName + "/newInterfaceAreaVectors");
  file.write(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.write(mInterfaceSmoothnessNormalization, pathName + "/interfaceSmoothnessNormalization");
  file.write(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
  file.write(mInterfaceFraction, pathName + "/interfaceFraction");
  file.write(mInterfaceAngles, pathName + "/interfaceAngles");
  
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mDamagedPressure, pathName + "/damagedPressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mVolume, pathName + "/volume");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");
  file.read(mPlasticStrain0, pathName + "/plasticStrain0");
  //file.read(mInverseEquivalentDeviatoricStress, pathName + "/inverseEquivalentDeviatoricStress");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mDPDx, pathName + "/DpDx");
  file.read(mDepsDx, pathName + "/DepsDx");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mM, pathName + "/M");
  file.read(mLocalM, pathName + "/localM");
  file.read(mNormalization, pathName + "/normalization");
  file.read(mInterfaceFlags, pathName + "/interfaceFlags");
  file.read(mInterfaceAreaVectors, pathName + "/interfaceAreaVectors");
  file.read(mInterfaceNormals, pathName + "/interfaceNormals");
  file.read(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.read(mNewInterfaceFlags, pathName + "/newInterfaceFlags");
  file.read(mNewInterfaceAreaVectors, pathName + "/newInterfaceAreaVectors");
  file.read(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.read(mInterfaceSmoothnessNormalization, pathName + "/interfaceSmoothnessNormalization");
  file.read(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
  file.read(mInterfaceFraction, pathName + "/interfaceFraction");
  file.read(mInterfaceAngles, pathName + "/interfaceAngles");

  // For backwards compatibility on change 3597 -- drop in the near future.
  for (auto DvDtPtr: mDvDt) DvDtPtr->name(HydroFieldNames::hydroAcceleration);

}

//------------------------------------------------------------------------------
// method for limited linear reconstruction between nodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPH<Dimension>::
linearReconstruction(const typename Dimension::Vector& ri,
                     const typename Dimension::Vector& rj,
                     const typename Dimension::Scalar& yi,
                     const typename Dimension::Scalar& yj,
                     const typename Dimension::Vector& DyDxi,
                     const typename Dimension::Vector& DyDxj,
                           typename Dimension::Scalar& ytildei,
                           typename Dimension::Scalar& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);
  const auto Dyj = 0.5*DyDxj.dot(rij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0) * std::max(tiny,abs(Dy0)));
  const auto xi = Dyi * denom;
  const auto xj = Dyj * denom;

  // limiter function - vanleer 1979
  const auto phii = ( xi > 0.0 ?  min(4.0*xi/((1.0 + xi)*(1.0 + xi)),1.0) : 0.0 );
  const auto phij = ( xj > 0.0 ?  min(4.0*xj/((1.0 + xj)*(1.0 + xj)),1.0) : 0.0 );        
  const auto phi = 0.5*(phii+phij);

  // linear constructed inteface values
  ytildei = yi - phi * Dyi;
  ytildej = yj + phi * Dyj;
}

} // Spheral namespace



