//---------------------------------Spheral++----------------------------------//
// Physics -- The topmost abstract base class for all physics packages in
// Spheral++.  Physics defines the methods generic to all physics classes,
// but leaves specialized interfaces for particular types of physics (radiative vs.
// hydro, etc.) to a descendent interface class, Physics.
//----------------------------------------------------------------------------//
#include "Physics.hh"
#include "Boundary/Boundary.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
Physics<Dimension>::
Physics():
  mBoundaryConditions(),
  mPreSubPackages(),
  mPostSubPackages() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
Physics<Dimension>::
~Physics() {
}

//------------------------------------------------------------------------------
// Add a Boundary condition to the end of the current boundary list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
appendBoundary(Boundary<Dimension>& boundary) {
  mBoundaryConditions.push_back(&boundary);
  for (auto* pkg: mPreSubPackages) pkg->appendBoundary(boundary);
  for (auto* pkg: mPostSubPackages) pkg->appendBoundary(boundary);
}

//------------------------------------------------------------------------------
// Add a Boundary condition to the beginning of the current boundary list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
prependBoundary(Boundary<Dimension>& boundary) {
  mBoundaryConditions.insert(mBoundaryConditions.begin(), &boundary);
  for (auto* pkg: mPreSubPackages) pkg->prependBoundary(boundary);
  for (auto* pkg: mPostSubPackages) pkg->prependBoundary(boundary);
}

//------------------------------------------------------------------------------
// Clear (erase) the boundary condition list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
clearBoundaries() {
  mBoundaryConditions = vector<Boundary<Dimension>*>();
  for (auto* pkg: mPreSubPackages) pkg->clearBoundaries();
  for (auto* pkg: mPostSubPackages) pkg->clearBoundaries();
}

//------------------------------------------------------------------------------
// Test if the given Boundary condition is listed in the physics package.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Physics<Dimension>::
haveBoundary(const Boundary<Dimension>& boundary) const {
  return std::count(mBoundaryConditions.begin(), mBoundaryConditions.end(), &boundary) > 0;
}

//------------------------------------------------------------------------------
// Provide a default no-op ghost boundary method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
applyGhostBoundaries(State<Dimension>& /*state*/,
                     StateDerivatives<Dimension>& /*derivs*/) {
}

//------------------------------------------------------------------------------
// Provide a default no-op enforce boundary method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
enforceBoundaries(State<Dimension>& /*state*/,
                  StateDerivatives<Dimension>& /*derivs*/) {
}

//------------------------------------------------------------------------------
// Add a physics package to be run after this one
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
appendSubPackage(Physics<Dimension>& package) {
  mPostSubPackages.push_back(&package);
}

//------------------------------------------------------------------------------
// Add a physics package to be run before this one
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
prependSubPackage(Physics<Dimension>& package) {
  mPreSubPackages.push_back(&package);
}

//------------------------------------------------------------------------------
// The set of packages to be run after this one
//------------------------------------------------------------------------------
template<typename Dimension>
const std::vector<Physics<Dimension>*>&
Physics<Dimension>::
postSubPackages() const {
  return mPostSubPackages;
}

//------------------------------------------------------------------------------
// The set of packages to be run before this one
//------------------------------------------------------------------------------
template<typename Dimension>
const std::vector<Physics<Dimension>*>&
Physics<Dimension>::
preSubPackages() const {
  return mPreSubPackages;
}

//------------------------------------------------------------------------------
// Provide a default no-op problem startup initialization method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
initializeProblemStartup(DataBase<Dimension>& /*dataBase*/) {
}

//------------------------------------------------------------------------------
// Provide a default no-op problem startup initialization dependencies method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& /*dataBase*/,
                                     State<Dimension>& /*state*/,
                                     StateDerivatives<Dimension>& /*derivs*/) {
}

//------------------------------------------------------------------------------
// Provide a default no-op pre-step initialization method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
preStepInitialize(const DataBase<Dimension>& /*dataBase*/, 
                  State<Dimension>& /*state*/,
                  StateDerivatives<Dimension>& /*derivs*/) {
}

//------------------------------------------------------------------------------
// Provide a default no-op initialization method.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Physics<Dimension>::
initialize(const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const DataBase<Dimension>& /*dataBase*/,
           State<Dimension>& /*state*/,
           StateDerivatives<Dimension>& /*derivs*/) {
  return false;
}

//------------------------------------------------------------------------------
// Provide a default no-op finalization method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
finalize(const typename Dimension::Scalar /*time*/,
         const typename Dimension::Scalar /*dt*/,
         DataBase<Dimension>& /*dataBase*/,
         State<Dimension>& /*state*/,
         StateDerivatives<Dimension>& /*derivs*/) {
}

//------------------------------------------------------------------------------
// Provide a default no-op finalizeDerivatives method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Physics<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& /*derivs*/) const {
}

//------------------------------------------------------------------------------
// Provide a default no-op postStateUpdate method.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Physics<Dimension>::
postStateUpdate(const Scalar /*time*/, 
                const Scalar /*dt*/,
                const DataBase<Dimension>& /*dataBase*/, 
                State<Dimension>& /*state*/,
                StateDerivatives<Dimension>& /*derivatives*/) {
  return false;
}

}
