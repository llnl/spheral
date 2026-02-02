//---------------------------------Spheral++----------------------------------//
// PointPotential -- Impose a potential from a point mass.
//
// Created by JMO, Sun Mar 30 22:08:55 PST 2003
//----------------------------------------------------------------------------//
#include "PointPotential.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PointPotential<Dimension>::
PointPotential(double G, double mass, double coreRadius, 
               const typename Dimension::Vector origin,
               const typename Dimension::Tensor metric):
  GenericBodyForce<Dimension>(),
  mG(G),
  mMass(mass),
  mCoreRadius2(coreRadius*coreRadius),
  mftimestep(0.1),
  mOrigin(origin),
  mMetric(metric),
  mDtMinAcc(0.0),
  mTotalPotentialEnergy(0.0),
  mPotential(FieldStorageType::CopyFields) {
  ENSURE(mG > 0.0);
  ENSURE(mMass >= 0.0);
  ENSURE(mCoreRadius2 >= 0.0);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PointPotential<Dimension>::
~PointPotential() {
}

//------------------------------------------------------------------------------
// Register some extra state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PointPotential<Dimension>::
registerState(DataBase<Dimension >& dataBase,
              State<Dimension >& state) {
  GenericBodyForce<Dimension >::registerState(dataBase, state);
  state.enroll(mPotential);
}

//------------------------------------------------------------------------------
// Calculate the acceleration due to the point mass on the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PointPotential<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the node positions from the state.
  const auto mnode = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto numNodeLists = position.numFields();

  // Get the acceleration and position change vectors we'll be modifying.
  auto DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero());
  auto DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero());

  // Loop over the internal nodes.
  mPotential = 0.0;
  mTotalPotentialEnergy = 0.0;
  mDtMinAcc = std::numeric_limits<Scalar>::max();
  const auto rcore = std::sqrt(mCoreRadius2);
  const Scalar thpt = mG*mMass;
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = mPotential[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto r = mMetric*(position(k, i) - mOrigin);
      const auto runit = r.unitVector();
      const auto rsoft2 = r.magnitude2() + mCoreRadius2;
      const auto rsoft  = sqrt(rsoft2);
      const auto rsoft3 = rsoft2*rsoft;
      const auto dphi = -thpt*mnode(k, i)/rsoft;
      const auto ai = -thpt*r.magnitude()/rsoft3*runit;
      CHECK(rsoft3 > 0.0);
      CHECK(rsoft  > 0.0);
      DxDt(k, i) += velocity(k, i);
      DvDt(k, i) += ai;
      mPotential(k, i) += dphi;
      mTotalPotentialEnergy += dphi;
      mDtMinAcc = std::min(mDtMinAcc, sqrt(rcore/ai.magnitude()));      // Similar to acceleration constraint from TreeGravity
    }
  }
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
PointPotential<Dimension>::
initializeProblemStartup(DataBase<Dimension>& db) {
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
PointPotential<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& db,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  // We need to make a dry run through setting derivatives and such
  // to set our initial vote on the time step.
  vector<Physics<Dimension>*> packages(1, this);
  this->initialize(0.0, 1.0, db, state, derivs);
  this->evaluateDerivatives(0.0, 1.0, db, state, derivs);
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename PointPotential<Dimension>::TimeStepType
PointPotential<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const {

  // A standard N-body approach -- just take the ratio of softening length/acceleration.
  const auto dt = mftimestep * mDtMinAcc;
  std::stringstream reasonStream;
  reasonStream << "PointPotential: f*sqrt(L/a) = " << dt << std::endl;
  return TimeStepType(dt, reasonStream.str());
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PointPotential<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mPotential, pathName + "/potential");
  file.write(mDtMinAcc, pathName + "/dtMinAcc");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PointPotential<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mPotential, pathName + "/potential");
  file.read(mDtMinAcc, pathName + "/dtMinAcc");
}

}
