//---------------------------------Spheral++----------------------------------//
// MFV -- This is an Arbitrary Eulerian-Lagrangian extension of the
//        MFV approach of Hopkins 2015. Its got several node-motion
//        approaches which promote more regular particle distributions.
//
//        Each of the ALE options defines the velocity of the nodes 
//        differently. The flux that results from the difference
//        between the nodes velocities and the fluid velocity.
//        The velocities are defined as follows for the 
//        NodeMotionTypes:
//
//        1) Eulerian ---- static Nodes
//        2) Lagrangian -- nodal velocity = fluid velocity. (This is
//                         a spheralized version of MFV so there
//                         is some flux between nodes)
//        3) Fician ------ nodal velocity = fluid velocity + Fician
//                         PST correction
//        4) XSPH -------- nodal velocity = xsph velocity
//
//   Hopkins P.F. (2015) "A New Class of Accurate, Mesh-Free Hydrodynamic 
//   Simulation Methods," MNRAS, 450(1):53-110
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
// TODO:
//   1 backpressure and fician particle shifting
//   2 Eulerian model will still crash on the Noh implosion due to void particles
//   3 Good implementation of Ngb update
//   4 treatment for material interfaces
//---------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/PureReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceWithRatioPolicy.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"

#include "GSPH/MFV.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/computeSumVolume.hh"
#include "GSPH/computeMFMDensity.hh"
#include "GSPH/Policies/MassFluxPolicy.hh"
#include "GSPH/Policies/MFVIncrementSpecificThermalEnergyPolicy.hh"
#include "GSPH/Policies/MFVIncrementVelocityPolicy.hh"
#include "GSPH/Policies/CompatibleMFVSpecificThermalEnergyPolicy.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <sstream>

using std::string;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MFV<Dimension>::
MFV(DataBase<Dimension>& dataBase,
    RiemannSolverBase<Dimension>& riemannSolver,
    const TableKernel<Dimension>& W,
    const Scalar epsDiffusionCoeff,
    const double cfl,
    const bool useVelocityMagnitudeForDt,
    const bool compatibleEnergyEvolution,
    const bool evolveTotalEnergy,
    const bool XSPH,
    const bool correctVelocityGradient,
    const double nodeMotionCoefficient,
    const NodeMotionType nodeMotionType,
    const GradientType gradType,
    const MassDensityType densityUpdate,
    const double epsTensile,
    const double nTensile,
    const Vector& xmin,
    const Vector& xmax):
  GenericRiemannHydro<Dimension>(dataBase,
                                 riemannSolver,
                                 W,
                                 epsDiffusionCoeff,
                                 cfl,
                                 useVelocityMagnitudeForDt,
                                 compatibleEnergyEvolution,
                                 evolveTotalEnergy,
                                 XSPH,
                                 correctVelocityGradient,
                                 gradType,
                                 densityUpdate,
                                 epsTensile,
                                 nTensile,
                                 xmin,
                                 xmax),
  mNodeMotionCoefficient(nodeMotionCoefficient),
  mNodeMotionType(nodeMotionType),
  mNodalVelocity(FieldStorageType::CopyFields),
  mDmassDt(FieldStorageType::CopyFields),
  mDthermalEnergyDt(FieldStorageType::CopyFields),
  mDmomentumDt(FieldStorageType::CopyFields),
  mDvolumeDt(FieldStorageType::CopyFields),
  //mHStretchTensor(FieldStorageType::CopyFields),
  mPairMassFluxPtr() {
    // mNodalVelocity = dataBase.newFluidFieldList(Vector::zero(), GSPHFieldNames::nodalVelocity);
    mDmassDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::mass);
    mDthermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy);
    mDmomentumDt = dataBase.newFluidFieldList(Vector::zero(), IncrementState<Dimension, Vector>::prefix() + GSPHFieldNames::momentum);
    mDvolumeDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::volume);
    //mHStretchTensor = dataBase.newFluidFieldList(SymTensor::zero(), "HStretchTensor");
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  GenericRiemannHydro<Dimension>::initializeProblemStartup(dataBase);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  GenericRiemannHydro<Dimension>::registerState(dataBase,state);

  // dataBase.resizeFluidFieldList(mNodalVelocity, Vector::zero(), GSPHFieldNames::nodalVelocity,false);

  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto position = state.fields(HydroFieldNames::position,Vector::zero());
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);

 // We use the thermal energy to update the specific thermal energy
  state.removePolicy(specificThermalEnergy,false);

  CHECK(position.numFields() == dataBase.numFluidNodeLists());
  CHECK(velocity.numFields() == dataBase.numFluidNodeLists());
  CHECK(volume.numFields() == dataBase.numFluidNodeLists());
  CHECK(mass.numFields() == dataBase.numFluidNodeLists());
  CHECK(specificThermalEnergy.numFields() == dataBase.numFluidNodeLists());

  auto nodeListi = 0u;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr < dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    auto& massi = (*itr)->mass();
    auto  minVolume = massi.min()/(*itr)->rhoMax();
    auto  maxVolume = massi.max()/(*itr)->rhoMin();
    state.enroll(*volume[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>(minVolume,
                                                                                           maxVolume));
  }

  
  state.enroll(massDensity, make_policy<ReplaceWithRatioPolicy<Dimension,Scalar>>({HydroFieldNames::mass,
                                                                                   HydroFieldNames::volume},
                                                                                  HydroFieldNames::mass,
                                                                                  HydroFieldNames::volume));

  state.enroll(mass,  make_policy<MassFluxPolicy<Dimension>>({HydroFieldNames::velocity, 
                                                              HydroFieldNames::specificThermalEnergy}));

  state.enroll(velocity, 
               make_policy<MFVIncrementVelocityPolicy<Dimension>>({HydroFieldNames::specificThermalEnergy}));

  
  if (this->compatibleEnergyEvolution()) {
    auto thermalEnergyPolicy = make_policy<CompatibleMFVSpecificThermalEnergyPolicy<Dimension>>(dataBase);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }else if (this->evolveTotalEnergy()) {
    std::cout <<"evolve total energy not implemented for MFV" << std::endl;
  } else {
    auto thermalEnergyPolicy = make_policy<MFVIncrementSpecificThermalEnergyPolicy<Dimension>>();
    state.enroll(specificThermalEnergy,thermalEnergyPolicy);
  }

}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);
  dataBase.resizeFluidFieldList(mDmassDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::mass, false);
  dataBase.resizeFluidFieldList(mDthermalEnergyDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy, false);
  dataBase.resizeFluidFieldList(mDmomentumDt, Vector::zero(), IncrementState<Dimension, Vector>::prefix() + GSPHFieldNames::momentum, false);
  dataBase.resizeFluidFieldList(mDvolumeDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::volume, false);
  //dataBase.resizeFluidFieldList(mHStretchTensor,SymTensor::zero(), "HStretchTensor", false);
  derivs.enroll(mDmassDt);
  derivs.enroll(mDthermalEnergyDt);
  derivs.enroll(mDmomentumDt);
  derivs.enroll(mDvolumeDt);
  //derivs.enroll(mHStretchTensor);
  const auto& connectivityMap = dataBase.connectivityMap();
  mPairMassFluxPtr = std::make_unique<PairMassFluxType>(connectivityMap);
  derivs.enroll(GSPHFieldNames::pairMassFlux, *mPairMassFluxPtr);
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::preStepInitialize(dataBase,state,derivs);

  if(this->densityUpdate() == MassDensityType::RigorousSumDensity){
    
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
          auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
          auto  volume = state.fields(HydroFieldNames::volume, 0.0);

    computeSumVolume(dataBase.connectivityMap(),this->kernel(),position,H,volume);
    computeMFMDensity(mass,volume,massDensity);
  
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr){
      (*boundaryItr)->applyFieldListGhostBoundary(volume);
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    }
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr < this->boundaryEnd(); 
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MFV<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  return GenericRiemannHydro<Dimension>::initialize(time,dt,dataBase,state,derivs); 
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  // hackish solution and I should be ashamed.
  if (this->compatibleEnergyEvolution()) {
    auto DpDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + GSPHFieldNames::momentum, Vector::zero());
    auto DmDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::mass, 0.0);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr){
      (*boundaryItr)->applyFieldListGhostBoundary(DpDt);
      (*boundaryItr)->applyFieldListGhostBoundary(DmDt);
  }
    
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
MFV<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);

  // auto nodalVelocity = state.fields(GSPHFieldNames::nodalVelocity, Vector::zero());

  // for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //      boundaryItr != this->boundaryEnd();
  //      ++boundaryItr) {
  //   (*boundaryItr)->applyFieldListGhostBoundary(nodalVelocity);
  // }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);

  // auto nodalVelocity = state.fields(GSPHFieldNames::nodalVelocity, Vector::zero());

  // for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //      boundaryItr != this->boundaryEnd();
  //      ++boundaryItr) {
  //   (*boundaryItr)->enforceFieldListBoundary(nodalVelocity);
  // }
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  GenericRiemannHydro<Dimension>::dumpState(file,pathName);
  // file.write(mNodalVelocity, pathName + "/nodalVelocity");
  file.write(mDmassDt, pathName + "/DmassDt");
  file.write(mDthermalEnergyDt, pathName + "/DthermalEnergyDt");
  file.write(mDmomentumDt, pathName + "/DmomentumDt");
  file.write(mDvolumeDt, pathName + "/DvolumeDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  GenericRiemannHydro<Dimension>::restoreState(file,pathName);
  // file.read(mNodalVelocity, pathName + "/nodalVelocity");
  file.read(mDmassDt, pathName + "/DmassDt");
  file.read(mDthermalEnergyDt, pathName + "/DthermalEnergyDt");
  file.read(mDmomentumDt, pathName + "/DmomentumDt");
  file.read(mDvolumeDt, pathName + "/DvolumeDt");
}

}

