//---------------------------------Spheral++----------------------------------//
// SmoothingScaleBase
//
// Abstract base class for packages that advance the smoothing scale.
//----------------------------------------------------------------------------//
#include "SmoothingScale/FixedSmoothingScale.hh"
#include "SmoothingScale/SmoothingScaleUtilities.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SmoothingScaleBase<Dimension>::
SmoothingScaleBase(const HEvolutionType HUpdate,
                   const bool fixShape,
                   const bool radialOnly):
  Physics<Dimension>(),
  mFixShape(fixShape),
  mRadialOnly(radialOnly),
  mHEvolution(HUpdate),
  mHideal(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mRadius0(FieldStorageType::CopyFields),
  mHidealFilterPtr(std::make_shared<ASPHSmoothingScaleUserFilter<Dimension>>()),
  mRadialFunctorPtr(std::make_shared<ASPHRadialFunctor<Dimension>>()),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Make sure our FieldLists are correctly sized.
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero(), ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero(), IncrementBoundedState<Dimension, Vector>::prefix() + HydroFieldNames::H, false);
  if (mRadialOnly) dataBase.resizeFluidFieldList(mRadius0, 0.0, "Start of step radius", false);
}

//------------------------------------------------------------------------------
// Register state
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  auto Hfields = dataBase.fluidHfield();
  const auto numFields = Hfields.numFields();
  for (auto k = 0u; k < numFields; ++k) {
    auto& Hfield = *Hfields[k];
    const auto& nodeList = Hfield.nodeList();
    const auto hmaxInv = 1.0/nodeList.hmax();
    const auto hminInv = 1.0/nodeList.hmin();
    switch (mHEvolution) {
      case HEvolutionType::IntegrateH:
        state.enroll(Hfield, make_policy<IncrementBoundedState<Dimension, SymTensor, Scalar>>(hmaxInv, hminInv));
        break;

      case HEvolutionType::IdealH:
        state.enroll(Hfield, make_policy<ReplaceBoundedState<Dimension, SymTensor, Scalar>>(hmaxInv, hminInv));
        break;

      case HEvolutionType::FixedH:
        state.enroll(Hfield);
        break;

       default:
         VERIFY2(false, "SmoothingScaleBase ERROR: Unknown Hevolution option ");
    }
  }
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  if (mHEvolution != HEvolutionType::FixedH) {
    derivs.enroll(mHideal);
    derivs.enroll(mDHDt);
  }
}

//------------------------------------------------------------------------------
// Initialize at the beginning of a timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // If we're using the radial H scaling, take a snapshot of the initial radius of
  // each point.
  if (mRadialOnly) {
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
    const auto numFields = pos.numFields();
    for (auto k = 0u; k < numFields; ++k) {
      const auto n = pos[k]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        mRadius0(k,i) = mRadialFunctorPtr->radialCoordinate(k, i, pos(k,i));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mHideal, pathName + "/Hideal");
  file.write(mDHDt, pathName + "/DHDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mHideal, pathName + "/Hideal");
  file.read(mDHDt, pathName + "/DHDt");
}

}
