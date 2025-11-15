//---------------------------------Spheral++----------------------------------//
// MFM -- spheralized verions of "Meshless Finite Mass" 
//   Hopkins P.F. (2015) "A New Class of Accurate, Mesh-Free Hydrodynamic 
//   Simulation Methods," MNRAS, 450(1):53-110
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceWithRatioPolicy.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"

#include "GSPH/MFM.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/computeSumVolume.hh"
#include "GSPH/computeMFMDensity.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <sstream>

using std::string;
using std::min;
using std::max;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MFM<Dimension>::
MFM(DataBase<Dimension>& dataBase,
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
  mDvolumeDt(FieldStorageType::CopyFields){
  mDvolumeDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::volume);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  GenericRiemannHydro<Dimension>::registerState(dataBase,state);
  
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
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

  auto massDensity = dataBase.fluidMassDensity();
  state.enroll(massDensity, make_policy<ReplaceWithRatioPolicy<Dimension,Scalar>>({HydroFieldNames::volume},
                                                                                  HydroFieldNames::mass,
                                                                                  HydroFieldNames::volume));
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);

  dataBase.resizeFluidFieldList(mDvolumeDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::volume, false);
  derivs.enroll(mDvolumeDt);
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::preStepInitialize(dataBase,state,derivs);

  if(this->densityUpdate() == MassDensityType::RigorousSumDensity){
    // plop into an intialize volume function
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
MFM<Dimension>::
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
MFM<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  GenericRiemannHydro<Dimension>::finalizeDerivatives(time,dt,dataBase,state,derivs);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  GenericRiemannHydro<Dimension>::dumpState(file,pathName);

  file.write(mDvolumeDt, pathName + "/DvolumeDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFM<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  GenericRiemannHydro<Dimension>::restoreState(file,pathName);

  file.read(mDvolumeDt, pathName + "/DvolumeDt");
}

}

