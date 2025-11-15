//---------------------------------Spheral++----------------------------------//
// GSPH -- A Riemann-solver-based implementation of SPH. Compared to 
//         MFM/MFV this approach requires a larger neighbor set. 2.5
//         nodes per kernel extent instead of 2-2.25 for MFM/MFV but
//         does perform better on certain tests (Noh implosion)
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "Hydro/computeSPHVolume.hh"
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
#include "Utilities/timingUtilities.hh"
#include "Utilities/Timer.hh"

#include "GSPH/GSPH.hh"
#include "GSPH/GSPHFieldNames.hh"
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
GSPH<Dimension>::
GSPH(DataBase<Dimension>& dataBase,
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
  mDmassDensityDt(FieldStorageType::CopyFields){

  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("GSPHinitializeStartup");
  GenericRiemannHydro<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);
  TIME_END("GSPHinitializeStartup");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("GSPHregister");

  GenericRiemannHydro<Dimension>::registerState(dataBase,state);

  auto massDensity = dataBase.fluidMassDensity();
  auto nodeListi = 0u;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr < dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    state.enroll(*massDensity[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>((*itr)->rhoMin(),
                                                                                                (*itr)->rhoMax()));
  }

  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  state.enroll(volume, make_policy<ReplaceWithRatioPolicy<Dimension,Scalar>>({HydroFieldNames::massDensity},
                                                                             HydroFieldNames::mass,
                                                                             HydroFieldNames::massDensity));
                                                                             

  TIME_END("GSPHregister");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("GSPHregisterDerivs");

  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);

  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  derivs.enroll(mDmassDensityDt);

  TIME_END("GSPHregisterDerivs");
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("GSPHpreStepInitialize");
  GenericRiemannHydro<Dimension>::preStepInitialize(dataBase,state,derivs);

  if(this->densityUpdate() == MassDensityType::RigorousSumDensity){
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
          auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
          auto  volume = state.fields(HydroFieldNames::volume, 0.0);

    computeSPHSumMassDensity(connectivityMap, this->kernel(), true, position, mass, H, massDensity);
    computeSPHVolume(mass,massDensity,volume);
  
    for (auto boundaryItr = this->boundaryBegin(); 
              boundaryItr < this->boundaryEnd(); 
              ++boundaryItr){
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      (*boundaryItr)->applyFieldListGhostBoundary(volume);
    }
    for (auto boundaryItr = this->boundaryBegin(); 
              boundaryItr < this->boundaryEnd(); 
              ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
  TIME_END("GSPHpreStepInitialize");
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GSPH<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("GSPHinitialize");
  auto result = GenericRiemannHydro<Dimension>::initialize(time,dt,dataBase,state,derivs);
  TIME_END("GSPHinitialize");
  return result;
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("GSPHfinalizeDerivs");
  GenericRiemannHydro<Dimension>::finalizeDerivatives(time,dt,dataBase,state,derivs);
  TIME_END("GSPHfinalizeDerivs");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("GSPHghostBounds");
  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);
  TIME_END("GSPHghostBounds");
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("GSPHenforceBounds");
  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);
  TIME_END("GSPHenforceBounds");
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  GenericRiemannHydro<Dimension>::dumpState(file,pathName);

  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPH<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  GenericRiemannHydro<Dimension>::restoreState(file,pathName);

  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
}

}
