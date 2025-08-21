//---------------------------------Spheral++----------------------------------//
// SolidFSISPHRZ
//
// RZ version of SolidFSISPH
// This RZ version is a naive area-weighting implementation.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
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
#include "Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "RK/ContinuityVolumePolicyRZ.hh"

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

#include "Geometry/GeometryRegistrar.hh"

#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/Timer.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/range.hh"

#include "FSISPH/SolidFSISPHRZ.hh"
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

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SolidFSISPHRZ::
SolidFSISPHRZ(DataBase<Dimension>& dataBase,
            ArtificialViscosityHandle<Dimension>& Q,
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
  SolidFSISPH<Dimension>(dataBase,
	      Q,
	      slides,
	      W,
	      cfl,
	      surfaceForceCoefficient,
	      densityStabilizationCoefficient,
	      specificThermalEnergyDiffusionCoefficient,
	      xsphCoefficient,
	      interfaceMethod,
	      kernelAveragingMethod,
	      sumDensityNodeLists,
	      useVelocityMagnitudeForDt,
	      compatibleEnergyEvolution,
	      evolveTotalEnergy,
	      linearCorrectGradients,
	      decoupleDamagedMaterial,
	      interfacePmin,
	      interfaceNeighborAngleThreshold,
	      densityUpdate,
	      epsTensile,
	      nTensile,
	      xmin,
	      xmax),
  mPairAccelerationsPtr(),
  mSelfAccelerations(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
// RZ area-weighting at startup for consistency.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  auto mass = dataBase.fluidMass();
  const auto pos = dataBase.fluidPosition();
  const unsigned numNodeLists = mass.numFields();
  // Divide mass by 2*pi*r
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }
  // Call base method
  SolidFSISPH<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);
  // Multiply mass back by 2*pi*r
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Register states
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  // The base class does most of it.
  SolidFSISPH<Dimension>::registerState(dataBase, state);

  // XXX TODO: add the following in, but its currently causing errors when run in (compatibleEnergy) mode
  // RuntimeError: Verification failed: StateBase ERROR: failed to return type for key pair-wise accelerations

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy.
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    state.enroll(specificThermalEnergy, make_policy<RZNonSymmetricSpecificThermalEnergyPolicy>(dataBase));

  } else if (evolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
registerDerivatives(DataBase<Dimension>&  dataBase,
                    StateDerivatives<Dimension>& derivs) {
  SolidFSISPH<Dimension>::registerDerivatives(dataBase, derivs);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    dataBase.resizeFluidFieldList(mSelfAccelerations, Vector::zero, HydroFieldNames::selfAccelerations, false);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
    derivs.enroll(mSelfAccelerations);
  }
}

//------------------------------------------------------------------------------
// FSI specialized density summmation
// RZ area-weighting: divide mass by 2*pi*r before, multiply back after base call.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
preStepInitialize(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Convert the mass to mass per unit length first.
  // Map: x->z, y->r (r = cylindrical radius)
  if (mApplySelectDensitySum){
    auto mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    const auto numNodeLists = mass.numFields();
    // Divide mass by 2*pi*r
    for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = mass[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const auto circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
        mass(nodeListi, i) /= circi;
      }
    }
  }

  // Base class does most of the work.
  SolidFSISPH<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Now convert back to true masses. (multiply by 2*pi*r)
  if (mApplySelectDensitySum){
    const auto position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
    auto       mass = state.fields(HydroFieldNames::mass, 0.0);
    auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const auto numNodeLists = massDensity.numFields();
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = mass[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const auto& xi = position(nodeListi, i);
        const auto  circi = 2.0*M_PI*abs(xi.y());
        mass(nodeListi, i) *= circi;
      }
    }
  }
}

//------------------------------------------------------------------------------
// For compatible energy we need to apply the bc conditions to acceleration
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
finalizeDerivatives(const Scalar time, 
                    const Scalar  dt,
                    const DataBase<Dimension>&  dataBase, 
                    const State<Dimension>& state,
		    StateDerivatives<Dimension>&  derivs) const {
  SolidFSISPH<Dimension>::finalizeDerivatives(time, dt, dataBase, state, derivs);
} // finalize


//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary SolidFSISPH BCs.
  SolidFSISPH<Dimension>::applyGhostBoundaries(state, derivs);
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Scale back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary SolidFSISPH BCs.
  SolidFSISPH<Dimension>::enforceBoundaries(state, derivs);

  // Scale back to mass.
  // We also ensure no point approaches the z-axis too closely.
  FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    //const Scalar nPerh = mass[nodeListi]->nodeList().nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      Vector& posi = pos(nodeListi, i);
      const Scalar circi = 2.0*M_PI*abs(posi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
dumpState(FileIO& file, const string& pathName) const {
  SolidFSISPH<Dimension>::dumpState(file, pathName);
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
restoreState(const FileIO& file, const string& pathName) {
  SolidFSISPH<Dimension>::restoreState(file, pathName);
}

//------------------------------------------------------------------------------
// method for limited linear reconstruction between nodes
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
linearReconstruction(const typename Dimension::Vector& ri,
                     const typename Dimension::Vector& rj,
                     const typename Dimension::Scalar& yi,
                     const typename Dimension::Scalar& yj,
                     const typename Dimension::Vector& DyDxi,
                     const typename Dimension::Vector& DyDxj,
                           typename Dimension::Scalar& ytildei,
                           typename Dimension::Scalar& ytildej) const {
  // XXX TODO: do I need 2*pi*r factor conversions here too?
  SolidFSISPH<Dimension>::linearReconstruction(ri, rj, yi, yj,
					    DyDxi, DyDxj, ytildei, ytildej);
}

//------------------------------------------------------------------------------
// Dispatch evaluateDerivatives based on type of Q
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
evaluateDerivatives(const Dimension::Scalar time,
                    const Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  this->firstDerivativesLoop(time,dt,dataBase,state,derivatives);

  // Depending on the type of the ArtificialViscosity, dispatch the call to
  // the secondDerivativesLoop
  auto& Qhandle = this->artificialViscosity();
  if (Qhandle.QPiTypeIndex() == std::type_index(typeid(Scalar))) {
      const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Scalar>&>(Qhandle);
      this->secondDerivativesLoop(time,dt,dataBase,state,derivatives,Q);
  } else {
    CHECK(Qhandle.QPiTypeIndex() == std::type_index(typeid(Tensor)));
    const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Tensor>&>(Qhandle);
    this->secondDerivativesLoop(time,dt,dataBase,state,derivatives,Q);
  }

  //this->setH(time,dt,dataBase,state,derivatves)
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename QType>
void
SolidFSISPHRZ::
secondDerivativesLoop(const Dimension::Scalar time,
                      const Dimension::Scalar dt,
                      const DataBase<Dimension>& dataBase,
                      const State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const QType& Q) const { 

  // XXX TODO: add R factors

  using QPiType = typename QType::ReturnType;

  // Get the SlideSurfaces.
  auto& slides = this->slideSurface();

  // The kernels and such.
  const auto& W = this->kernel();

  // huge amount of tinies
  const auto tiny = std::numeric_limits<double>::epsilon();
  const auto tinyScalarDamage = 1.0e-2;
  const auto tinyNonDimensional = 1.0e-9;

  // map to get the new interface flag from the neighbor details
  std::vector<std::vector<int>> interfaceFlagMap(2, vector<int>(6,0));
  interfaceFlagMap[0][2] = 1;
  interfaceFlagMap[0][4] = 1;
  interfaceFlagMap[1][0] = 3;
  interfaceFlagMap[1][1] = 3;
  interfaceFlagMap[1][2] = 3;
  interfaceFlagMap[1][3] = 3;
  interfaceFlagMap[1][4] = 3;

  // constants
  const auto W0 = W(0.0, 1.0);
  const auto interfaceNeighborAngleThreshold = this->interfaceNeighborAngleThreshold();
  const auto interfacePmin = this->interfacePmin();
  const auto decoupleDamagedMaterial = this->decoupleDamagedMaterial();
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto totalEnergy = this->evolveTotalEnergy();
  const auto epsDiffusionCoeff = this->specificThermalEnergyDiffusionCoefficient();
  const auto rhoStabilizeCoeff = this->densityStabilizationCoefficient();
  const auto surfaceForceCoeff = this->surfaceForceCoefficient();
  const auto xsphCoeff = this->xsphCoefficient();
  const auto XSPH = xsphCoeff > tiny;
  const auto diffuseEnergy = epsDiffusionCoeff>tiny and compatibleEnergy;
  const auto stabilizeDensity = rhoStabilizeCoeff>tiny;
  const auto alwaysAverageKernels = (mKernelAveragingMethod==KernelAveragingMethod::AlwaysAverageKernels);
  const auto averageInterfaceKernels = (mKernelAveragingMethod==KernelAveragingMethod::AverageInterfaceKernels);
  const auto constructHLLC = (mInterfaceMethod == InterfaceMethod::HLLCInterface);
  const auto activateConstruction = !(mInterfaceMethod == InterfaceMethod::NoInterface);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  numPairs = pairs.size();
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

  // Get the state and derivative FieldLists.
  const auto interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  const auto interfaceFlags = state.fields(FSIFieldNames::interfaceFlags, int(0));
  const auto interfaceAreaVectors = state.fields(FSIFieldNames::interfaceAreaVectors, Vector::zero);
  const auto interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto damagedPressure = state.fields(SolidFieldNames::damagedPressure, 0.0);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  const auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  const auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);

  //const auto yield = state.fields(SolidFieldNames::yieldStrength, 0.0);
  //const auto invJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);

  CHECK(interfaceFlags.size() == numNodeLists);
  CHECK(interfaceAreaVectors.size() == numNodeLists);
  CHECK(interfaceNormals.size() == numNodeLists);
  CHECK(interfaceSmoothness.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(damagedPressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(K.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);
  //CHECK(yield.size() == numNodeLists);
  //CHECK(invJ2.size() == numNodeLists);

  // Derivative FieldLists.
  const auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  localM = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  DepsDx = derivs.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  const auto  DPDx = derivs.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  newInterfaceNormals = derivs.fields(PureReplaceState<Dimension, Vector>::prefix() + FSIFieldNames::interfaceNormals, Vector::zero);
  auto  newInterfaceFlags = derivs.fields(PureReplaceState<Dimension, int>::prefix() + FSIFieldNames::interfaceFlags, int(0));
  auto  newInterfaceAreaVectors = derivs.fields(PureReplaceState<Dimension, Vector>::prefix() + FSIFieldNames::interfaceAreaVectors, Vector::zero);
  auto  interfaceSmoothnessNormalization = derivs.fields(FSIFieldNames::interfaceSmoothnessNormalization, 0.0);
  auto  interfaceFraction = derivs.fields(FSIFieldNames::interfaceFraction, 0.0);
  auto  newInterfaceSmoothness = derivs.fields(PureReplaceState<Dimension, Scalar>::prefix() + FSIFieldNames::interfaceSmoothness, 0.0);
  auto  interfaceAngles = derivs.fields(FSIFieldNames::interfaceAngles, 0.0);
  auto  normalization = derivs.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  DSDt = derivs.fields(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto* pairAccelerationsPtr = (compatibleEnergy ?
                                &derivs.template get<PairAccelerationsType>(HydroFieldNames::pairAccelerations) :
                                nullptr);
  auto* pairDepsDtPtr = (compatibleEnergy ?
                         &derivs.template get<PairWorkType>(HydroFieldNames::pairWork) :
                         nullptr);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(DepsDx.size() == numNodeLists);
  CHECK(DPDx.size() == numNodeLists);
  CHECK(newInterfaceFlags.size() == numNodeLists);
  CHECK(newInterfaceAreaVectors.size() == numNodeLists);
  CHECK(newInterfaceNormals.size() == numNodeLists);
  CHECK(interfaceSmoothnessNormalization.size() == numNodeLists);
  CHECK(interfaceFraction.size() == numNodeLists);
  CHECK(newInterfaceSmoothness.size() == numNodeLists);
  CHECK(interfaceAngles.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);
  CHECK(not compatibleEnergy or pairAccelerationsPtr->size() == numPairs);
  CHECK(not compatibleEnergy or pairDepsDtPtr->size() == numPairs);

  //this->computeMCorrection(time,dt,dataBase,state,derivs);

// Now we calculate  the hydro deriviatives
// Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private  scratch variables.
    unsigned i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj, PLineari, PLinearj, epsLineari, epsLinearj, Qi, Qj;
    QPiType QPiij, QPiji;
    SymTensor sigmai, sigmaj;
    Vector sigmarhoi, sigmarhoj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto newInterfaceFlags_thread = newInterfaceFlags.threadCopy(threadStack, ThreadReduction::MAX);
    auto newInterfaceAreaVectors_thread = newInterfaceAreaVectors.threadCopy(threadStack);
    auto newInterfaceNormals_thread = newInterfaceNormals.threadCopy(threadStack);
    auto newInterfaceSmoothness_thread = newInterfaceSmoothness.threadCopy(threadStack);
    auto interfaceSmoothnessNormalization_thread = interfaceSmoothnessNormalization.threadCopy(threadStack);
    auto interfaceFraction_thread = interfaceFraction.threadCopy(threadStack);
    auto interfaceAngles_thread = interfaceAngles.threadCopy(threadStack, ThreadReduction::MAX);
    auto normalization_thread = normalization.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);
    auto DSDt_thread = DSDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    
#pragma omp for
    for (auto kk = 0u; kk < numPairs; ++kk) {
      
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& interfaceFlagsi = interfaceFlags(nodeListi,i);
      const auto& interfaceAreaVectorsi = interfaceAreaVectors(nodeListi,i);
      const auto& interfaceNormalsi = interfaceNormals(nodeListi,i);
      const auto& interfaceSmoothnessi = interfaceSmoothness(nodeListi,i);
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi,i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Pdi = damagedPressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& pTypei = pTypes(nodeListi, i);
      const auto  fragIDi = fragIDs(nodeListi, i);
      //const auto  Yi = yield(nodeListi, i);
      //const auto  invJ2i = invJ2(nodeListi, i);
      const auto  voli = mi/rhoi;
      const auto  mui = max(mu(nodeListi,i),tiny);
      const auto  Ki = max(tiny,K(nodeListi,i))+4.0/3.0*mui;
      const auto  Hdeti = Hi.Determinant();
      epsLineari = epsi;
      PLineari = Pdi;

      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DepsDxi = DepsDx(nodeListi, i);
      const auto& DPDxi = DPDx(nodeListi, i);
      const auto& Mi = M(nodeListi, i);
      auto& normi = normalization_thread(nodeListi,i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& newInterfaceFlagsi = newInterfaceFlags_thread(nodeListi,i);
      auto& newInterfaceAreaVectorsi = newInterfaceAreaVectors_thread(nodeListi,i);
      auto& newInterfaceNormalsi = newInterfaceNormals_thread(nodeListi,i);
      auto& newInterfaceSmoothnessi = newInterfaceSmoothness_thread(nodeListi,i);
      auto& interfaceFractioni = interfaceFraction_thread(nodeListi,i);
      auto& interfaceSmoothnessNormalizationi = interfaceSmoothnessNormalization_thread(nodeListi,i);
      auto& minNeighborAnglei = interfaceAngles_thread(nodeListi,i);
      
      // Get the state for node j
      const auto& interfaceFlagsj = interfaceFlags(nodeListj,j);
      const auto& interfaceAreaVectorsj = interfaceAreaVectors(nodeListj,j);
      const auto& interfaceNormalsj = interfaceNormals(nodeListj,j);
      const auto& interfaceSmoothnessj = interfaceSmoothness(nodeListj,j);
      const auto& rj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj,j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Pdj = damagedPressure(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto& pTypej = pTypes(nodeListj, j);
      const auto  fragIDj = fragIDs(nodeListj, j);
      //const auto  Yj = yield(nodeListj, j);
      //const auto  invJ2j = invJ2(nodeListj, j);
      const auto  volj = mj/rhoj;
      const auto  muj = max(mu(nodeListj,j),tiny);
      const auto  Kj = max(tiny,K(nodeListj,j))+4.0/3.0*muj;
      const auto  Hdetj = Hj.Determinant();
      epsLinearj = epsj;
      PLinearj = Pdj;

      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      const auto& DepsDxj = DepsDx(nodeListj, j);
      const auto& DPDxj = DPDx(nodeListj, j);
      const auto& Mj = M(nodeListj,j);
      auto& normj = normalization_thread(nodeListj,j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& newInterfaceFlagsj = newInterfaceFlags_thread(nodeListj,j);
      auto& newInterfaceAreaVectorsj = newInterfaceAreaVectors_thread(nodeListj,j);
      auto& newInterfaceNormalsj = newInterfaceNormals_thread(nodeListj,j);
      auto& newInterfaceSmoothnessj = newInterfaceSmoothness_thread(nodeListj,j);
      auto& interfaceFractionj = interfaceFraction_thread(nodeListj,j);
      auto& interfaceSmoothnessNormalizationj = interfaceSmoothnessNormalization_thread(nodeListj,j);
      auto& minNeighborAnglej = interfaceAngles_thread(nodeListj,j);

      // line of action
      const auto rij = ri - rj;
      const auto rhatij = rij.unitVector();
      
      // decoupling and boolean switches
      //-------------------------------------------------------
      // Flag if this is a contiguous material pair or not.
      const auto sameMatij =  (nodeListi == nodeListj and fragIDi==fragIDj);
      const auto differentMatij = !sameMatij; 
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);
      
      // pairwise damage and nodal damage
      const auto Di = max(0.0, min(1.0, damage(nodeListi, i).dot(rhatij).magnitude()));
      const auto Dj = max(0.0, min(1.0, damage(nodeListj, j).dot(rhatij).magnitude()));
      const auto fDi =  (sameMatij ? (1.0-Di)*(1.0-Di) : 0.0 );
      const auto fDj =  (sameMatij ? (1.0-Dj)*(1.0-Dj) : 0.0 );
      const auto fDij = (sameMatij ? pow(1.0-std::abs(Di-Dj),2.0) : 0.0 );

      // is Pmin being activated? (Pmin -> interface Pmin)
      const auto pLimiti = (sameMatij ? (Pdi-rhoi*ci*ci*tinyNonDimensional) : interfacePmin);
      const auto pLimitj = (sameMatij ? (Pdj-rhoj*cj*cj*tinyNonDimensional) : interfacePmin);
      const auto pminActivei = (Pi < pLimiti);
      const auto pminActivej = (Pj < pLimitj);
      
      // decoupling criteria, we want material interface to be able to separate and if 
      // decoupleDamagedMaterial is active we want damaged material to behave like gravel
      const auto isExpanding = (ri-rj).dot(vi-vj) > 0.0;
      const auto isFullyDamaged = (fDi<tinyScalarDamage) or (fDj<tinyScalarDamage);
      const auto isPastAdhesionThreshold = pminActivei or pminActivej;
      const auto canDecouple = (isFullyDamaged and decoupleDamagedMaterial) or differentMatij;

      const auto decouple = isExpanding  and isPastAdhesionThreshold and canDecouple;

      // do we need to construct our interface velocity?
      const auto constructInterface = (fDij < 1.0-tinyScalarDamage) and activateConstruction;
      const auto negligableShearWave = max(mui,muj) < tinyNonDimensional*min(Ki,Kj);

      // do we reduce our deviatoric stress
      const auto isTensile = (((Si+Sj)-(Pdi+Pdj)*SymTensor::one).dot(rhatij)).dot(rhatij) > 0;
      const auto damageReduceStress = isTensile or differentMatij;
    
      // Kernels
      //--------------------------------------
      const auto Hij = 0.5*(Hi+Hj);
      const auto etaij = Hij*rij;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagij = etaij.magnitude();
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagij >= 0.0);
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      const auto Hetai = Hi*etai.unitVector();
      const auto Hetaj = Hj*etaj.unitVector();
      auto gradWi = gWi*Hetai;
      auto gradWj = gWj*Hetaj;
      auto gradWiMi = gradWi;
      auto gradWjMj = gradWj;

      // average our kernels
      const auto gradWij = 0.5*(gradWi+gradWj);
      const auto Wij = 0.5*(Wi+Wj);

      if(averageKernelij){
        const auto gWij = 0.5*(gWi+gWj);
        Wi = Wij;
        Wj = Wij;
        gWi = gWij;
        gWj = gWij;
        gradWi = gradWij;
        gradWj = gradWij;
      }

      gradWiMi = Mi.Transpose()*gradWi;
      gradWjMj = Mj.Transpose()*gradWj;

      // interface fields
      //----------------------------------------------------------------------------------
      const auto fSij = ( sameMatij ? 1.0 : -1.0);
      const auto interfaceSwitch = std::min(std::min(interfaceFlagsj,interfaceFlagsi),1);
      const int sameMatMask = (sameMatij ? 1 : 0);
      const int diffMatMask = 1-sameMatMask;
      const int controlNodeMaski = interfaceFlagMap[0][interfaceFlagsi];
      const int controlNodeMaskj = interfaceFlagMap[0][interfaceFlagsj];

      const auto AijMij = fSij*(voli*volj)*(gradWjMj+gradWiMi); // surface area vector

      const auto alignment = (interfaceFlagsi == 5 or interfaceFlagsj == 5 ? 
                              0.0:
                              std::max(fSij*interfaceNormalsi.dot(interfaceNormalsj),0.0));

      newInterfaceAreaVectorsi -= AijMij;
      newInterfaceAreaVectorsj += AijMij;

      interfaceFractioni += sameMatMask * volj * Wi;
      interfaceFractionj += sameMatMask * voli * Wj;

      // we use the angle between normal and neighbors to check for interface control particles (type 2 or 4)
      const auto minNeighborMaski = (interfaceFlagsj != 5 and sameMatij ? 1.0 : -1.0); // turns following into no-op
      const auto minNeighborMaskj = (interfaceFlagsi != 5 and sameMatij ? 1.0 : -1.0);
      minNeighborAnglei = std::max(std::min(minNeighborMaski,-interfaceAreaVectorsi.unitVector().dot(rhatij)),minNeighborAnglei);
      minNeighborAnglej = std::max(std::min(minNeighborMaskj, interfaceAreaVectorsj.unitVector().dot(rhatij)),minNeighborAnglej);

      // if the node neighbors a control node, set its value to 1 if initially 0
      newInterfaceFlagsi = std::max(newInterfaceFlagsi, interfaceFlagMap[diffMatMask][interfaceFlagsj]);
      newInterfaceFlagsj = std::max(newInterfaceFlagsj, interfaceFlagMap[diffMatMask][interfaceFlagsi]);

      // only flagged nodes get normals and they're based on their control node neighbors
      newInterfaceNormalsi += interfaceSwitch*controlNodeMaskj*fSij*volj*interfaceAreaVectorsj*Wij;
      newInterfaceNormalsj += interfaceSwitch*controlNodeMaski*fSij*voli*interfaceAreaVectorsi*Wij;

      // smoothness is calculated for all flagged interface/surface nodes
      interfaceSmoothnessNormalizationi += interfaceSwitch*volj*Wij;
      interfaceSmoothnessNormalizationj += interfaceSwitch*voli*Wij;
      newInterfaceSmoothnessi += interfaceSwitch*alignment*volj*Wij;
      newInterfaceSmoothnessj += interfaceSwitch*alignment*voli*Wij;

      if (!decouple){

        // Stress state
        //---------------------------------------------------------------
        const auto rhoij = 0.5*(rhoi+rhoj); 
        const auto cij = 0.5*(ci+cj); 
        const auto vij = vi - vj;

        // raw AV
        Q.QPiij(QPiij, QPiji, Qi, Qj,
                nodeListi, i, nodeListj, j,
                ri, Hij, etaij, vi, rhoij, cij,  
                rj, Hij, etaij, vj, rhoij, cij,
                fClQ, fCqQ, DvDxQ); 

        // slide correction
        if (slides.isSlideSurface(nodeListi,nodeListj)){
          const auto slideCorr = slides.slideCorrection(interfaceSmoothnessi,
                                                        interfaceSmoothnessj,
                                                        interfaceNormalsi,
                                                        interfaceNormalsj,
                                                        vi,
                                                        vj);
          QPiij *= slideCorr;
          QPiji *= slideCorr;
        }

        // save our max pressure from the AV for each node
        const auto conversionFactorQ = rhoi*rhoj/(max(rhoij,tiny)*max(rhoij,tiny));
        Qi *= conversionFactorQ;
        Qj *= conversionFactorQ;
        maxViscousPressurei = max(maxViscousPressurei, Qi);
        maxViscousPressurej = max(maxViscousPressurej, Qj);
        effViscousPressurei += volj * Qi * Wi;
        effViscousPressurej += voli * Qj * Wj;

        // stress tensor
        //{
          // apply yield pairwise 
          //const auto Yij = std::max(0.0,std::min(Yi,Yj));
          //const auto fYieldi = std::min(Yij*invJ2i,1.0);
          //const auto fYieldj = std::min(Yij*invJ2j,1.0);
          //const auto Seffi = (sameMatij ? 1.0 : 0.0 ) * Si;
          //const auto Seffj = (sameMatij ? 1.0 : 0.0 ) * Sj;
          //const auto Seffi = (min(Pi,Pj) < 0.0 or differentMatij ? maxfDij : 1.0)  * Si;
          //const auto Seffj = (min(Pi,Pj) < 0.0 or differentMatij ? maxfDij : 1.0)  * Sj;
          //const auto Peffi = maxfDij*min(Pi,0.0) + max(Pi,0.0);
          //const auto Peffj = maxfDij*min(Pj,0.0) + max(Pj,0.0);
          //const auto Seffi = maxfDij * Si; //(damageReduceStress ? fDij : 1.0) * Si;
          //const auto Seffj = maxfDij * Sj; //(damageReduceStress ? fDij : 1.0) * Sj;
          const auto Peffi = (differentMatij ? max(Pdi,interfacePmin) : Pdi);
          const auto Peffj = (differentMatij ? max(Pdj,interfacePmin) : Pdj);
          const auto Seffi = (damageReduceStress ? fDij : 1.0) * Si;
          const auto Seffj = (damageReduceStress ? fDij : 1.0) * Sj;
          sigmai = Seffi - Peffi * SymTensor::one;
          sigmaj = Seffj - Peffj * SymTensor::one;
        //}

        // Compute the tensile correction to add to the stress as described in 
        // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
        {
          const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
          const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
          const auto Ri = fi*tensileStressCorrection(sigmai);
          const auto Rj = fj*tensileStressCorrection(sigmaj);
          sigmai += Ri;
          sigmaj += Rj;
        }

        // accelerations
        //---------------------------------------------------------------
        const auto rhoirhoj = 1.0/(rhoi*rhoj);
        const auto sf = (sameMatij ? 1.0 : 1.0 + surfaceForceCoeff*abs((rhoi-rhoj)/(rhoi+rhoj+tiny)));
        sigmarhoi = sf*(rhoirhoj*sigmai*gradWiMi - 0.5*QPiij*gradWiMi);
        sigmarhoj = sf*(rhoirhoj*sigmaj*gradWjMj - 0.5*QPiji*gradWjMj);

        if (averageKernelij){
          const auto sigmarhoij = 0.5*(sigmarhoi+sigmarhoj);
          sigmarhoi = sigmarhoij;
          sigmarhoj = sigmarhoij;
        }
      
        const auto deltaDvDt = sigmarhoi + sigmarhoj;

        if (freeParticle) {
          DvDti += mj*deltaDvDt;
          DvDtj -= mi*deltaDvDt;
        } 
      
        // Velocity Gradient
        //-----------------------------------------------------------
        // construct our interface velocity 
        auto vstar = 0.5*(vi+vj);

        linearReconstruction(ri,rj,Pdi,Pdj,DPDxi,DPDxj,PLineari,PLinearj);
        if (constructInterface){
          // components
          const auto ui = vi.dot(rhatij);
          const auto uj = vj.dot(rhatij);
          const auto wi = vi - ui*rhatij;
          const auto wj = vj - uj*rhatij;
          
          // weights weights
          const auto Ci =  (constructHLLC ? std::sqrt(rhoi*Ki)  : Ki  ) + tiny;
          const auto Cj =  (constructHLLC ? std::sqrt(rhoj*Kj)  : Kj  ) + tiny;
          const auto Csi = (constructHLLC ? std::sqrt(rhoi*mui) : mui ) + tiny;
          const auto Csj = (constructHLLC ? std::sqrt(rhoj*muj) : muj ) + tiny;
          const auto CiCjInv = safeInv(Ci+Cj,tiny);
          const auto CsiCsjInv = safeInv(Csi+Csj,tiny);

          // weights
          const auto weightUi = max(0.0, min(1.0, Ci*CiCjInv));
          const auto weightUj = 1.0 - weightUi;
          const auto weightWi = (negligableShearWave ? weightUi : max(0.0, min(1.0, Csi*CsiCsjInv )) );
          const auto weightWj = 1.0 - weightWi;

          // interface velocity
          const auto ustar = weightUi*ui + weightUj*uj + (constructHLLC ? (PLinearj - PLineari)*CiCjInv : 0.0); 
          const auto wstar = weightWi*wi + weightWj*wj;// - (constructHLLC ? (Seffj - Seffi).dot(rhatij)*CsiCsjInv : Vector::zero);
          vstar = fDij * vstar + (1.0-fDij)*(ustar*rhatij + wstar);
        }

        // TODO:
        //------------------------------------------------------------
        // fix pressure/rawPressure convert to damagePressure/pressure
        // get wave speeds in there in a good manner
        //------------------------------------------------------------

        // local velocity gradient for DSDt
        if (sameMatij) {
          localDvDxi -=  2.0*volj*((vi-vstar).dyad(gradWi));
          localDvDxj -=  2.0*voli*((vstar-vj).dyad(gradWj)); 
        }
        
        // diffuse to stabilize things
        if (stabilizeDensity and (ci>tiny and cj>tiny)){
          const auto cFactor = 1.0 + max(min( (vi-vj).dot(rhatij)/max(cij,tiny), 0.0), -1.0);
          const auto effCoeff = (differentMatij ? 1.0 : rhoStabilizeCoeff*cFactor);
          vstar += (constructHLLC ? fDij : 1.0) * effCoeff * rhatij * cij * min(max((PLinearj-PLineari)/(Ki + Kj),-0.25),0.25);
        }

        // global velocity gradient
        DvDxi -= 2.0*volj*(vi-vstar).dyad(gradWiMi);
        DvDxj -= 2.0*voli*(vstar-vj).dyad(gradWjMj);

        // energy conservation
        // ----------------------------------------------------------
        const auto deltaDepsDti = 2.0*sigmarhoi.dot(vi-vstar);
        const auto deltaDepsDtj = 2.0*sigmarhoj.dot(vstar-vj);

        DepsDti -= mj*deltaDepsDti;
        DepsDtj -= mi*deltaDepsDtj;

        if(compatibleEnergy){
          (*pairAccelerationsPtr)[kk] = - deltaDvDt;
          (*pairDepsDtPtr)[kk][0] = - deltaDepsDti; 
          (*pairDepsDtPtr)[kk][1] = - deltaDepsDtj;
        }
        
        // thermal diffusion
        //-----------------------------------------------------------
        if (sameMatij and diffuseEnergy){
          linearReconstruction(ri,rj,epsi,epsj,DepsDxi,DepsDxj,epsLineari,epsLinearj);
          const auto cijEff = max(min(cij + (vi-vj).dot(rhatij), cij),0.0);
          const auto diffusion =  epsDiffusionCoeff*cijEff*(epsLineari-epsLinearj)*etaij.dot(gradWij)/(rhoij*etaMagij*etaMagij+tiny);
          if (compatibleEnergy) {
            (*pairDepsDtPtr)[kk][0] += diffusion; 
            (*pairDepsDtPtr)[kk][1] -= diffusion;
          }
        }

        // normalization 
        //-----------------------------------------------------------
        normi += volj*Wi;
        normj += voli*Wj;

        // XSPH -- we use this to handle tensile instability here
        //-----------------------------------------------------------
        if (sameMatij and XSPH) {
          const auto fxsph = (std::min(Pdi,Pdj) < 0 ? 1 : 0);
          XSPHWeightSumi += fxsph * volj*Wi;
          XSPHWeightSumj += fxsph * voli*Wj;
          XSPHDeltaVi -= fDij*volj*Wi*(vi-vj);
          XSPHDeltaVj -= fDij*voli*Wj*(vj-vi);
        }

      } // if damageDecouple 
    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  } // OpenMP parallel region


  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto& interfaceFlagsi = interfaceFlags(nodeListi,i);
      const auto& interfaceAreaVectorsi = interfaceAreaVectors(nodeListi,i);
      const auto  Hdeti = Hi.Determinant();
      const auto psi = Hdeti*mi/rhoi*W0;
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DvDti = DvDt(nodeListi,i);
      const auto& localMi = localM(nodeListi, i);
      auto& normi = normalization(nodeListi,i);
      auto& DepsDti = DepsDt(nodeListi,i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      auto& newInterfaceNormalsi = newInterfaceNormals(nodeListi,i);
      auto& newInterfaceSmoothnessi = newInterfaceSmoothness(nodeListi,i);
      auto& interfaceSmoothnessNormalizationi = interfaceSmoothnessNormalization(nodeListi,i);
      auto& interfaceFractioni = interfaceFraction(nodeListi,i);
      auto& newInterfaceFlagsi = newInterfaceFlags(nodeListi,i);
      auto& minNeighborAnglei = interfaceAngles(nodeListi,i);

      // finish our normalization
      normi += psi;

      // finish our interface fields.
      if(minNeighborAnglei < interfaceNeighborAngleThreshold)  newInterfaceFlagsi = std::max(2,newInterfaceFlagsi+1);
      if(interfaceFractioni/(1.0-psi)<0.1) newInterfaceFlagsi = 5.0;

      interfaceFractioni += psi;

      if ( interfaceFlagsi > 0 ){
        const double proxWeighti = 100*(1.0 - interfaceFlagsi % 2);
        newInterfaceNormalsi = (newInterfaceNormalsi + proxWeighti * psi * interfaceAreaVectorsi).unitVector();
        newInterfaceSmoothnessi = newInterfaceSmoothnessi/max(interfaceSmoothnessNormalizationi,tiny);
      } else {
        newInterfaceNormalsi = Vector::zero;
      }

      DrhoDti -=  rhoi*DvDxi.Trace();

      if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      DxDti = vi;
      if (XSPH) {
        CHECK(normi >= 0.0);
        const auto invNormi = 1.0/max(normi,tiny);
        XSPHWeightSumi = (interfaceFlagsi == 0 ? XSPHWeightSumi*invNormi : 0);
        DxDti += xsphCoeff*XSPHWeightSumi*XSPHDeltaVi*invNormi;
      }

      localDvDxi = localDvDxi*localMi;

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - deformation.Trace()/3.0*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti += spinCorrection + 2.0*mui*deviatoricDeformation;
      
    } //loop-nodes
  } //loop-nodeLists
} // evaluateDerivatives method


//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
void
SolidFSISPHRZ::
firstDerivativesLoop(const Dimension::Scalar /*time*/,
                     const Dimension::Scalar /*dt*/,
                     const DataBase<Dimension>& dataBase,
                     const State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) const {

  // XXX TODO: add R factors

  // The kernels and such.
  const auto& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  const auto alwaysAverageKernels = (mKernelAveragingMethod==KernelAveragingMethod::AlwaysAverageKernels);
  const auto averageInterfaceKernels = (mKernelAveragingMethod==KernelAveragingMethod::AverageInterfaceKernels);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  numPairs = pairs.size();

  // Get the state and derivative FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto damagedPressure = state.fields(SolidFieldNames::damagedPressure, 0.0);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(damagedPressure.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DepsDx = derivs.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  auto  DPDx = derivs.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  
  CHECK(DepsDx.size() == numNodeLists);
  CHECK(DPDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);

#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto DPDx_thread = DPDx.threadCopy(threadStack);
    auto DepsDx_thread = DepsDx.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < numPairs; ++kk) {

      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& fragIDi = fragIDs(nodeListi,i);
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = damagedPressure(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      // Get the state for node j
      const auto& fragIDj = fragIDs(nodeListj,j);
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = damagedPressure(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DPDxi = DPDx_thread(nodeListi, i);
      auto& DPDxj = DPDx_thread(nodeListj, j);
      auto& DepsDxi = DepsDx_thread(nodeListi, i);
      auto& DepsDxj = DepsDx_thread(nodeListj, j);
      auto& localMi = localM_thread(nodeListi,i);
      auto& localMj = localM_thread(nodeListj,j);
      auto& Mi = M_thread(nodeListi,i);
      auto& Mj = M_thread(nodeListj,j);

      const auto rij = ri - rj;
      const auto Pij = Pi - Pj;
      const auto epsij = epsi - epsj;

      // logic
      //---------------------------------------
      const auto sameMatij = (nodeListi == nodeListj and fragIDi == fragIDj);
      const auto differentMatij = (nodeListi!=nodeListj);
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);

      // Kernels
      //--------------------------------------
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      const auto gWi = W.gradValue(etaMagi, Hdeti);
      const auto gWj = W.gradValue(etaMagj, Hdetj);
      
      const auto Hetai = Hi*etai.unitVector();
      const auto Hetaj = Hj*etaj.unitVector();
      
      auto gradWi = gWi*Hetai;
      auto gradWj = gWj*Hetaj;
      
      //Wi & Wj --> Wij for interface better agreement DrhoDt and DepsDt
      if(averageKernelij){
        const auto gradWij = 0.5*(gradWi+gradWj);
        gradWi = gradWij;
        gradWj = gradWij;
      }

      gradWi *= mj/rhoj;
      gradWj *= mi/rhoi;

      // spatial gradients and correction
      //---------------------------------------------------------------
      const auto deltaRi = rij.dyad(gradWi);
      const auto deltaRj = rij.dyad(gradWj);

      Mi -= deltaRi;
      Mj -= deltaRj;

      DPDxi -= Pij*gradWi;
      DPDxj -= Pij*gradWj;

      if(sameMatij){
        localMi -=  deltaRi;
        localMj -=  deltaRj;
        DepsDxi -= epsij*gradWi;
        DepsDxj -= epsij*gradWj;
      }
    } // loop over pairs
      // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

   
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto& nodeList = mass[nodeListi]->nodeList();
      const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
        auto& Mi = M(nodeListi, i);
        auto& localMi = localM(nodeListi, i);
        auto& DepsDxi = DepsDx(nodeListi, i);
        auto& DPDxi = DPDx(nodeListi, i);

        const auto Mdeti = Mi.Determinant();
        const auto goodM = ( Mdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
        Mi =  (goodM and this->linearCorrectGradients() ? Mi.Inverse() : Tensor::one);
        
        const auto localMdeti = localMi.Determinant();
        const auto goodLocalM = ( localMdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
        localMi =  (goodLocalM and this->linearCorrectGradients() ? localMi.Inverse() : Tensor::one);

        DPDxi = Mi.Transpose()*DPDxi;
        DepsDxi = localMi.Transpose()*DepsDxi;
      } // for each node
    }   // for each nodelist

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(M);
      (*boundaryItr)->applyFieldListGhostBoundary(DPDx);
      (*boundaryItr)->applyFieldListGhostBoundary(DepsDx);
    }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

 } // method 

} // Spheral namespace

