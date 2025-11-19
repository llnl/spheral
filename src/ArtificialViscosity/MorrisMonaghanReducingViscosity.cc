//---------------------------------Spheral++----------------------------------//
// A simple form of the reducing artificial viscosity from Morris & Monaghan.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementBoundedState.hh"

#include "MorrisMonaghanReducingViscosity.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {
    
//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MorrisMonaghanReducingViscosity<Dimension>::
MorrisMonaghanReducingViscosity(const Scalar nhQ,
                                const Scalar nhL,
                                const Scalar aMin,
                                const Scalar aMax,
                                const Scalar negligibleSoundSpeed):
  Physics<Dimension>(),
  mnhQ(nhQ),
  mnhL(nhL),
  maMin(aMin),
  maMax(aMax),
  mNegCs(negligibleSoundSpeed),
  mClMultiplier(FieldStorageType::CopyFields),
  mCqMultiplier(FieldStorageType::CopyFields),
  mDrvAlphaDtQ(FieldStorageType::CopyFields),
  mDrvAlphaDtL(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  mClMultiplier = dataBase.newFluidFieldList(0.0, HydroFieldNames::ArtificialViscousClMultiplier);
  mCqMultiplier = dataBase.newFluidFieldList(0.0, HydroFieldNames::ArtificialViscousCqMultiplier);
  mDrvAlphaDtQ = dataBase.newFluidFieldList(0.0, IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::ArtificialViscousCqMultiplier);
  mDrvAlphaDtL = dataBase.newFluidFieldList(0.0, IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::ArtificialViscousClMultiplier);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
registerState(DataBase<Dimension>& /*dataBase*/,
              State<Dimension>& state) {
  state.enroll(mCqMultiplier, make_policy<IncrementBoundedState<Dimension, Scalar>>(maMin,maMax));
  state.enroll(mClMultiplier, make_policy<IncrementBoundedState<Dimension, Scalar>>(maMin,maMax));
}
    
//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDrvAlphaDtQ);
  derivs.enroll(mDrvAlphaDtL);
}

//------------------------------------------------------------------------------
// Boundaries
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  auto ClMult = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  auto CqMult = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  for (auto* bcPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    bcPtr->applyFieldListGhostBoundary(ClMult);
    bcPtr->applyFieldListGhostBoundary(CqMult);
  }
}
    
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  // Get Qtys for Derivs
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto Hsmooth = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto reducingViscosityMultiplierQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  const auto reducingViscosityMultiplierL = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
    
  // Derivative FieldLists
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto DrvAlphaDtQ = derivs.fields(IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  auto DrvAlphaDtL = derivs.fields(IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
    
  auto numNodeLists = pressure.numFields();
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(reducingViscosityMultiplierQ.size() == numNodeLists);
  CHECK(reducingViscosityMultiplierL.size() == numNodeLists);
  CHECK(DrvAlphaDtQ.size() == numNodeLists);
  CHECK(DrvAlphaDtL.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(Hsmooth.size() == numNodeLists);
        
  // Walk the nodes
  for (auto i = 0u; i < numNodeLists; ++i) {
    const auto numNodes = pressure[i]->numInternalElements();
    for (auto j = 0u; j < numNodes; j++) {
      const Scalar rvQ = reducingViscosityMultiplierQ(i,j);
      const Scalar rvL = reducingViscosityMultiplierL(i,j);
      const Scalar cs = soundSpeed(i,j);
      const Scalar csSafe = (cs*cs + mNegCs*mNegCs)/cs;
      const Scalar pmin = max(mNegCs*mNegCs*massDensity(i,j),abs(pressure(i,j)));
      const Scalar source = max(-DvDx(i,j).Trace(),0.0);
      const Scalar adiabatIndex = max(cs*cs*massDensity(i,j)/pmin,1.0+mNegCs);
      const Scalar decayConstQ = (1.0/mnhQ)*sqrt((adiabatIndex-1.0)/(2.0*adiabatIndex));
      const Scalar decayConstL = (1.0/mnhL)*sqrt((adiabatIndex-1.0)/(2.0*adiabatIndex));
      const Scalar h = 1.0/(Dimension::rootnu(Hsmooth(i,j).Determinant()));
      const Scalar decayTimeQ = h/(decayConstQ*csSafe);
      const Scalar decayTimeL = h/(decayConstL*csSafe);
            
      // safeInverse instead of forcing cs != 0
      DrvAlphaDtQ(i,j) = (maMax-rvQ)*source - (rvQ - maMin)/decayTimeQ;
      DrvAlphaDtL(i,j) = (maMax-rvL)*source - (rvL - maMin)/decayTimeL;
            
      // shock detection switch
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mClMultiplier, pathName + "/ClMultiplier");
  file.write(mCqMultiplier, pathName + "/CqMultiplier");
  file.write(mDrvAlphaDtQ, pathName + "/DrvAlphaDtQ");
  file.write(mDrvAlphaDtL, pathName + "/DrvAlphaDtL");
}
    
//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mClMultiplier, pathName + "/ClMultiplier");
  file.read(mCqMultiplier, pathName + "/CqMultiplier");
  file.read(mDrvAlphaDtQ, pathName + "/DrvAlphaDtQ");
  file.read(mDrvAlphaDtL, pathName + "/DrvAlphaDtL");
}
    
//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename MorrisMonaghanReducingViscosity<Dimension>::TimeStepType
MorrisMonaghanReducingViscosity<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/,
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return TimeStepType(1.0e100, "Rate of viscosity change -- NO VOTE.");
}
    
}
