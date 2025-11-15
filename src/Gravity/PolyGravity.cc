//---------------------------------Spheral++----------------------------------//
// PolyGravity -- Solve gravity on a polytope (2D or 3D).
//
// Based on Jason Pearl's approximate gravity model.  Currently 2D polygons
// are not really correct, as they use the 3D implementation and therefore
// don't properly reprsent the 2D infinite-rod logarithmic potential.
//
// Created by JMO, Fri Sep 16 15:04:54 PDT 2022
//----------------------------------------------------------------------------//
#include "PolyGravity.hh"
#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/packElement.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/PairComparisons.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::sqrt;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {   // anonymous

//------------------------------------------------------------------------------
// Dim::Vector -> Vector3d
//------------------------------------------------------------------------------
Dim<3>::Vector TO3VEC(const Dim<3>::Vector& p) { return p; }
Dim<3>::Vector TO3VEC(const Dim<2>::Vector& p) { return Dim<3>::Vector(p.x(), p.y(), 0.0); }

//------------------------------------------------------------------------------
// Dim::Vector -> Vector3d
//------------------------------------------------------------------------------
template<typename Dimension> typename Dimension::Vector FROM3VEC(const Dim<3>::Vector& p) { return p; }
template<> Dim<2>::Vector FROM3VEC<Dim<2>>(const Dim<3>::Vector& p) { return Dim<2>::Vector(p.x(), p.y()); }

}             // anonymous
  
//------------------------------------------------------------------------------
// Constructor (3D)
//------------------------------------------------------------------------------
template<>
PolyGravity<Dim<3>>::
PolyGravity(const PolyGravity<Dim<3>>::Polytope& poly,
            const double G,
            const double mass,
            const double ftimestep,
            const GravityTimeStepType timeStepChoice):
  GenericBodyForce<Dim<3>>(),
  mG(G),
  mMass(mass),
  mftimestep(ftimestep),
  mDynamicalTime(1.0/sqrt(G*mass/poly.volume())),
  mTimeStepChoice(timeStepChoice),
  mPoly(poly),
  mSolver(std::make_shared<ApproximatePolyhedralGravityModel>(poly, mass, G)),
  mPotential(FieldStorageType::CopyFields),
  mExtraEnergy(0.0),
  mDtMinAcc(0.0),
  mRestart(registerWithRestart(*this)) {
  VERIFY2(G > 0.0, "PolyGravity requires G > 0");
  VERIFY2(mass > 0.0,"PolyGravity requires mass > 0");
  VERIFY2(ftimestep > 0.0, "PolyGravity requires ftimestep > 0");
}

//------------------------------------------------------------------------------
// Constructor (2D)
//------------------------------------------------------------------------------
template<>
PolyGravity<Dim<2>>::
PolyGravity(const PolyGravity<Dim<2>>::Polytope& poly,
            const double G,
            const double mass,
            const double ftimestep,
            const GravityTimeStepType timeStepChoice):
  GenericBodyForce<Dim<2>>(),
  mG(G),
  mMass(mass),
  mftimestep(ftimestep),
  mDynamicalTime(1.0/sqrt(G*mass/poly.volume())),
  mTimeStepChoice(timeStepChoice),
  mPoly(),
  mSolver(),
  mPotential(FieldStorageType::CopyFields),
  mExtraEnergy(0.0),
  mDtMinAcc(0.0),
  mRestart(registerWithRestart(*this)) {
  VERIFY2(G > 0.0, "PolyGravity requires G > 0");
  VERIFY2(mass > 0.0,"PolyGravity requires mass > 0");
  VERIFY2(ftimestep > 0.0, "PolyGravity requires ftimestep > 0");

  // Extract the 2d geometry.
  const auto& points2d = poly.vertices();
  const auto& facets2d = poly.facets();

  // Build the 3D vertices by extruding the 2d points up and down in z.
  const auto length = sqrt(poly.volume());
  const unsigned npoints2d = points2d.size();
  vector<Dim<3>::Vector> points3d(2*npoints2d);
  for (auto i = 0u; i < npoints2d; ++i) {
    points3d[i].x(points2d[i].x());
    points3d[i].y(points2d[i].y());
    points3d[i].z(-length);
    const auto j = i + npoints2d;
    points3d[j].x(points2d[i].x());
    points3d[j].y(points2d[i].y());
    points3d[j].z(length);
  }

  // Build the 3D facets
  vector<vector<unsigned>> facets3d;
  for (const auto& facet: facets2d) {
    facets3d.push_back(vector<unsigned>({facet.ipoint1(),
                                         facet.ipoint2(),
                                         facet.ipoint2() + npoints2d,
                                         facet.ipoint1() + npoints2d}));
  }
  CHECK(facets3d.size() == facets2d.size());

  // Close the top and bottom.
  vector<unsigned> top_facet, bottom_facet;
  for (auto i = 0u; i < npoints2d; ++i) {
    top_facet.push_back(npoints2d + i);
    bottom_facet.push_back(npoints2d - i);
  }
  facets3d.push_back(top_facet);
  facets3d.push_back(bottom_facet);

  // Build a 3D extruded version of the polytope
  mPoly = Dim<3>::FacetedVolume(points3d, facets3d);
  mSolver =  std::make_shared<ApproximatePolyhedralGravityModel>(mPoly, mass, G);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PolyGravity<Dimension>::
~PolyGravity() {
}

//------------------------------------------------------------------------------
// Register some extra state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolyGravity<Dimension>::
registerState(DataBase<Dimension >& dataBase,
              State<Dimension >& state) {
  GenericBodyForce<Dimension >::registerState(dataBase, state);
  state.enroll(mPotential);
}

//------------------------------------------------------------------------------
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
PolyGravity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension >& state,
                    StateDerivatives<Dimension >& derivs) const {

  // Access the pertinent fields in the database.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto numNodeLists = pos.numFields();

  // Get the acceleration and position change vectors we'll be modifying.
  auto DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero());
  auto DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero());

  // Walk the points and compute the acceleration and potential
  mPotential = 0.0;
  mDtMinAcc = std::numeric_limits<Scalar>::max();
  mExtraEnergy = 0.0;
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = mPotential[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto ai = FROM3VEC<Dimension>(mSolver->acceleration(TO3VEC(pos(k, i))));
      DxDt(k, i) = vel(k, i);
      DvDt(k, i) += ai;
      mPotential(k, i) = mSolver->potential(TO3VEC(pos(k, i)));
      mExtraEnergy += mPotential(k, i);
      const auto hi = 1.0/(Dimension::nDim * H(k, i).Trace());
      mDtMinAcc = min(mDtMinAcc, sqrt(hi/ai.magnitude()));      // Similar to acceleration constraint from TreeGravity
    }
  }
  mExtraEnergy = allReduce(mExtraEnergy, SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
PolyGravity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
PolyGravity<Dimension>::
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
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename PolyGravity<Dimension>::TimeStepType
PolyGravity<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {

  // A standard N-body approach -- just take the ratio of softening length/acceleration.
  if (mTimeStepChoice == GravityTimeStepType::AccelerationRatio) {
    const auto dt = mftimestep * mDtMinAcc;
    std::stringstream reasonStream;
    reasonStream << "PolyGravity: f*sqrt(L/a) = " << dt << std::endl;
    return TimeStepType(dt, reasonStream.str());

  } else {

    // Dynamical time
    const auto dt = mftimestep * mDynamicalTime;
    std::stringstream reasonStream;
    reasonStream << "PolyGravity: dynamical time = " << mDynamicalTime << std::endl;
    return TimeStepType(dt, reasonStream.str());
  }
}

//------------------------------------------------------------------------------
// extraEnergy
//------------------------------------------------------------------------------
template<typename Dimension>
typename PolyGravity<Dimension>::Scalar 
PolyGravity<Dimension>::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldList<Dimension, typename PolyGravity<Dimension>::Scalar>&
PolyGravity<Dimension>::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// The surface model
//------------------------------------------------------------------------------
template<typename Dimension>
const Dim<3>::FacetedVolume&
PolyGravity<Dimension>::
poly() const {
  return mPoly;
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
template<typename Dimension>
double
PolyGravity<Dimension>::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// mass
//------------------------------------------------------------------------------
template<typename Dimension>
double
PolyGravity<Dimension>::
mass() const {
  return mMass;
}

//------------------------------------------------------------------------------
// ftimestep
//------------------------------------------------------------------------------
template<typename Dimension>
double
PolyGravity<Dimension>::
ftimestep() const {
  return mftimestep;
}

template<typename Dimension>
void
PolyGravity<Dimension>::
ftimestep(double x) {
  VERIFY(x > 0.0);
  mftimestep = x;
}

//------------------------------------------------------------------------------
// timeStepChoice
//------------------------------------------------------------------------------
template<typename Dimension>
GravityTimeStepType
PolyGravity<Dimension>::
timeStepChoice() const {
  return mTimeStepChoice;
}

template<typename Dimension>
void
PolyGravity<Dimension>::
timeStepChoice(GravityTimeStepType x) {
  mTimeStepChoice = x;
}

//------------------------------------------------------------------------------
// dynamicalTime
//------------------------------------------------------------------------------
template<typename Dimension>
double
PolyGravity<Dimension>::
dynamicalTime() const {
  return mDynamicalTime;
}

//------------------------------------------------------------------------------
// The polyhedral gravity solver
//------------------------------------------------------------------------------
template<typename Dimension>
const ApproximatePolyhedralGravityModel&
PolyGravity<Dimension>::
solver() const {
  return *mSolver;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolyGravity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mPotential, pathName + "/potential");
  file.write(mDtMinAcc, pathName + "/pairWiseDtMin");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolyGravity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mPotential, pathName + "/potential");
  file.read(mDtMinAcc, pathName + "/pairWiseDtMin");
}

}
