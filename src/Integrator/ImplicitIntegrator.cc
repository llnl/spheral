//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrator -- Abstract base class for implicit in time integrators.
//
// Created by JMO, Fri Oct 18 16:34:11 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/ImplicitIntegrator.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/SpheralMessage.hh"
#include "FileIO/FileIO.hh"
#include "Distributed/Communicator.hh"

namespace Spheral {

using std::cout;
using std::endl;

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
ImplicitIntegrator(DataBase<Dimension>& dataBase,
                   const std::vector<Physics<Dimension>*>& physicsPackages,
                   const Scalar tol):
  Integrator<Dimension>(dataBase, physicsPackages),
  mTol(tol),
  mMaxAllowedDtMultiplier(100.0),
  mMaxGoodDtMultiplier(1.0),
  mNumExplicitSteps(0u),
  mNumImplicitSteps(0u) {
}

//------------------------------------------------------------------------------
// step
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ImplicitIntegrator<Dimension>::
step(const typename Dimension::Scalar maxTime) {

  DataBase<Dimension>& db = this->accessDataBase();
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  auto success = false;
  auto count = 0u;
  auto maxIterations = 10u;
  while (not success and count++ < maxIterations) {
    SpheralMessage("=============> Current dtMultiplier = " << mDtMultiplier << " " << mMaxGoodDtMultiplier << " " << mMaxAllowedDtMultiplier);
    
    // Try to advance using the current timestep multiplier
    success = this->step(maxTime, state, derivs);
    if (success) mMaxGoodDtMultiplier = mDtMultiplier;

    // Adjust the current timestep multiplier based on whether we succeeded or not
    mDtMultiplier *= (success ?
                      (mDtMultiplier < 0.9*mMaxGoodDtMultiplier ? 1.2 : 1.01) :
                      0.95);
    mDtMultiplier = min(mMaxAllowedDtMultiplier, mDtMultiplier);

    if (not success and this->verbose()) {
      SpheralMessage("ImplicitIntegrator::step did not converge with tolerance on iteration " << count << "/" << maxIterations << endl
                     << "                         reducing timestep multiplier to " << mDtMultiplier);
    }
  }
  return success;
}

//------------------------------------------------------------------------------
// Loop over the stored physics packages and pick the minimum timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ImplicitIntegrator<Dimension>::
selectDt(const typename Dimension::Scalar dtMin,
         const typename Dimension::Scalar dtMax,
         const State<Dimension>& state,
         const StateDerivatives<Dimension>& derivs) const {

  REQUIRE(dtMin >= 0 and dtMax > 0);
  REQUIRE(dtMin <= dtMax);

  // Get the current time and data base.
  auto        t = this->currentTime();
  const auto& db = this->dataBase();

  // Loop over each package, and pick their timesteps.
  TimeStepType dt(dtMax, ""), dtImp(dtMax, "");
  const auto& pkgs = this->physicsPackages();
  for (auto* pkg: pkgs) {
    auto dtVote = pkg->dt(db, state, derivs, t);
    if (dtVote.first > 0.0 and dtVote.first < dt.first) dt = dtVote;
    dtVote = pkg->dtImplicit(db, state, derivs, t);
    if (dtVote.first > 0.0 and dtVote.first < dtImp.first) dtImp = dtVote;
  }

  // In the parallel case we need to find the minimum timestep across all processors.
  const auto globalDt = allReduce(dt.first, SPHERAL_OP_MIN);
  const auto globalDtImp = allReduce(dtImp.first, SPHERAL_OP_MIN);

  // If we're verbose we have to know which rank(s) are limiting things
  const auto rank = Process::getRank();
  const auto numProcs = Process::getTotalNumberOfProcesses();
  auto dtRank    = (dt.first    == globalDt    ? rank : numProcs);
  auto dtRankImp = (dtImp.first == globalDtImp ? rank : numProcs);
  dtRank = allReduce(dtRank, SPHERAL_OP_MIN);
  dtRankImp = allReduce(dtRankImp, SPHERAL_OP_MIN);

  // Limit the curent mulitplier such that the timestep cannot exceed the implicit vote
  mDtMultiplier *= std::min(1.0, globalDtImp/globalDt);

  // Now set the dt to the global answer
  dt.first = globalDt;
#ifdef SPHERAL_ENABLE_MPI
  {
    int msgSize = dt.second.size();
    MPI_Bcast(&msgSize, 1, MPI_INT, dtRank, Communicator::communicator());
    dt.second.resize(msgSize);
    MPI_Bcast(&dt.second[0], msgSize, MPI_CHAR, dtRank, Communicator::communicator());
  }
#endif

  // Apply any dt scaling due to iteration
  dt.first *= mDtMultiplier;

  // We also require that the timestep is not allowed to grow faster than a
  // prescribed fraction.
  dt.first = std::min(dt.first, this->dtGrowth()*this->lastDt());

  // Enforce the timestep boundaries.
  dt.first = std::min(dtMax, std::max(dtMin, dt.first));

  CHECK(dt.first >= 0.0 and
        dt.first >= dtMin and dt.first <= dtMax);

  // Are we verbose?
  if (rank == dtRank and
      (this->verbose() or globalDt < this->dtMin())) {
    cout << "----------------------------------------" << endl
         << "Overall timestep chosen " << dt.first << endl
         << dt.second << endl;
  }
  if (rank == dtRankImp and
      (this->verbose() or globalDt < this->dtMin())) {
    cout << "Implicit limiting timestep " << dtImp.first << endl
         << dtImp.second << endl;
  }
  cout.flush();

  Barrier(Communicator::communicator());

  return dt.first;
}

//------------------------------------------------------------------------------
// Find the maximum residual difference in the states
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ImplicitIntegrator<Dimension>::
computeResiduals(const State<Dimension>& state1,
                 const State<Dimension>& state0,
                 const bool forceVerbose) const {

  // Get the local (to this MPI rank) answer
  const auto& db = this->dataBase();
  const auto& packages = this->physicsPackages();
  const auto  tol = this->convergenceTolerance();
  auto result = ResidualType(-1.0, "DEFAULT : You should not get this");
  for (const auto* pkg: packages) {
    const auto pres = pkg->maxResidual(db, state1, state0, tol);
    if (pres.first > result.first) result = pres;
  }

  // Reduce for the global result, and optionally print out some verbose info
  const auto globalMax = allReduce(result.first, SPHERAL_OP_MAX);
  if (result.first == globalMax and (forceVerbose or this->verbose())) {
    cout << "Global residual of "
         << result.first << endl
         << result.second << endl;
  }
  cout.flush();
  return globalMax;
}

//------------------------------------------------------------------------------
// Dump the current state of the Integrator to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ImplicitIntegrator<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  Integrator<Dimension>::dumpState(file, pathName);
  file.write(mMaxGoodDtMultiplier, pathName + "/maxGoodDtMultiplier");
  file.write(mNumExplicitSteps, pathName + "/numExplicitSteps");
  file.write(mNumImplicitSteps, pathName + "/numImplicitSteps");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ImplicitIntegrator<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  Integrator<Dimension>::restoreState(file, pathName);
  // When restarting from a non-implicit integrator these values may not be available.
  // This is not a fatal condition, so we make restarting these variables optional
  if (file.readIfAvailable(mMaxGoodDtMultiplier, pathName + "/maxGoodDtMultiplier") != 0)
    SpheralMessage("ImplicitIntegrator WARNING: unable to load old maxGoodDtMultiplier");
  if (file.readIfAvailable(mNumExplicitSteps, pathName + "/numExplicitSteps") != 0)
    SpheralMessage("ImplicitIntegrator WARNING: unable to load old numExplicitSteps");
  if (file.readIfAvailable(mNumImplicitSteps, pathName + "/numExplicitSteps") != 0)
    SpheralMessage("ImplicitIntegrator WARNING: unable to load old numImplicitSteps");
}

//------------------------------------------------------------------------------
// How should we query a physics package for the time step?
//------------------------------------------------------------------------------
template<typename Dimension>
typename ImplicitIntegrator<Dimension>::TimeStepType
ImplicitIntegrator<Dimension>::
dt(const Physics<Dimension>* pkg,
   const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return pkg->dt(dataBase, state, derivs, currentTime);
  // return pkg->dtImplicit(dataBase, state, derivs, currentTime);
}

}
