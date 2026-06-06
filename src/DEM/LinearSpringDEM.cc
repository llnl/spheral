
//---------------------------------Spheral++----------------------------------//
// DEM -- damped linear spring contact model based on pkdgrav immplementation
//----------------------------------------------------------------------------
// Schwartz, S.R. and Richards, D.C. "An implementation of the soft-sphere 
// discrete element method in a high-performance parallel gravity tree-code,"
// Granular Matter, (2012) 14:363–380, 10.1007/s10035-012-0346-z.
//
// Zhang et. al. "Rotational Failure of Rubble-pile Bodies: Influences of 
// Shear and Cohesive Strengths," The Astrophysical Journal, (2018) 857:15, 20
// 10.3847/1538-4357/aab5b2.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/MaxReplaceState.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"

#include "Boundary/Boundary.hh"

#include "DEM/ReplaceAndIncrementPairFieldList.hh"
#include "DEM/DEMFieldNames.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/LinearSpringDEM.hh"
#include "DEM/ContactStorageLocation.hh"
#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

#include "Utilities/Timer.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>

using std::string;
using std::make_pair;
using std::to_string;

namespace Spheral {

namespace {

// Provide to_string for Spheral Vector
template<typename Vector>
std::string
vec_to_string(const Vector& vec) {
  std::ostringstream oss;
  oss << vec << std::endl;
  return oss.str();
}

// function to calculate the contact durations
double
computeContactDuration(double m, 
                       double k, 
                       double B){
  CHECK(m >= 0);
  CHECK(k > 0);
  CHECK(B < 0);
  return M_PI*std::sqrt(0.5*m/k * (1.0 + 1.0/(B*B)));
}

// our lumped damping term were going to force the restitution coefficient to be less than
// 1 and more than zero with a small buffer to prevent singular behavior w/out resorting
// to if statements.
double
computeBeta(double restitutionCoefficient){
  CHECK(restitutionCoefficient >= 0);
  CHECK(restitutionCoefficient <= 1);
  const double tiny = 1.0e-8;
  return M_PI/std::log(std::min(std::max(restitutionCoefficient,tiny),1.0-tiny));
}

}



//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
LinearSpringDEM<Dimension>::
LinearSpringDEM(const DataBase<Dimension>& dataBase,
                const Scalar normalSpringConstant,
                const Scalar normalRestitutionCoefficient,
                const Scalar tangentialSpringConstant,
                const Scalar tangentialRestitutionCoefficient,
                const Scalar dynamicFrictionCoefficient,
                const Scalar staticFrictionCoefficient,
                const Scalar rollingFrictionCoefficient,
                const Scalar torsionalFrictionCoefficient,
                const Scalar cohesiveTensileStrength,
                const Scalar shapeFactor,
                const Scalar stepsPerCollision,
                const bool   enableFastTimeStepping,
                const Vector& xmin,
                const Vector& xmax):
  DEMBase<Dimension>(dataBase,stepsPerCollision,xmin,xmax),
  mEnableFastTimeStepping(enableFastTimeStepping),
  mNormalSpringConstant(normalSpringConstant),
  mNormalRestitutionCoefficient(normalRestitutionCoefficient),
  mTangentialSpringConstant(tangentialSpringConstant),
  mTangentialRestitutionCoefficient(tangentialRestitutionCoefficient),
  mDynamicFrictionCoefficient(dynamicFrictionCoefficient),
  mStaticFrictionCoefficient(staticFrictionCoefficient),
  mRollingFrictionCoefficient(rollingFrictionCoefficient),
  mTorsionalFrictionCoefficient(torsionalFrictionCoefficient),
  mCohesiveTensileStrength(cohesiveTensileStrength),
  mShapeFactor(shapeFactor),
  mNormalBeta(computeBeta(normalRestitutionCoefficient)),
  mTangentialBeta(computeBeta(tangentialRestitutionCoefficient)),
  mCollisionDuration(0.0),
  mMomentOfInertia(FieldStorageType::CopyFields),
  mMaximumOverlap(FieldStorageType::CopyFields),
  mNewMaximumOverlap(FieldStorageType::CopyFields) { 

    mMomentOfInertia = dataBase.newDEMFieldList(0.0, DEMFieldNames::momentOfInertia);
    mMaximumOverlap = dataBase.newDEMFieldList(0.0, DEMFieldNames::maximumOverlap);
    mNewMaximumOverlap = dataBase.newDEMFieldList(0.0,MaxReplaceState<Dimension, Scalar>::prefix() + DEMFieldNames::maximumOverlap);

    const auto mass = dataBase.DEMMass();
    const auto minMass = mass.min();
    mCollisionDuration = computeContactDuration(minMass,mNormalSpringConstant,mNormalBeta);

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
LinearSpringDEM<Dimension>::
~LinearSpringDEM() {}

//------------------------------------------------------------------------------
// time step -- constant for this model
//------------------------------------------------------------------------------
template<typename Dimension>
typename LinearSpringDEM<Dimension>::TimeStepType
LinearSpringDEM<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const typename Dimension::Scalar currentTime) const {

  TIME_BEGIN("LinearSpringDEMdt");

  auto dtMin = std::numeric_limits<Scalar>::max();
  TimeStepType result(dtMin, "DEM error, this message should not get to the end");

  if (this->enableFastTimeStepping()){
    result = this->variableTimeStep(dataBase,
                                    state,
                                    derivs,
                                    currentTime);
  }else{
    result = this->fixedTimeStep();
  }

  TIME_END("LinearSpringDEMdt");
  return result;
}

//------------------------------------------------------------------------------
// set our timestep (for now its a constant single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename LinearSpringDEM<Dimension>::TimeStepType
LinearSpringDEM<Dimension>::
fixedTimeStep() const {
    return make_pair(mCollisionDuration/this->stepsPerCollision(),("fixed-dt Linear Spring DEM vote for time step"));
}

template<typename Dimension>
void
LinearSpringDEM<Dimension>::
recomputeContactDuration()  {
  const auto& db = this->dataBase();
  const auto mass = db.DEMMass();
  const auto minMass = mass.min();
  mCollisionDuration = computeContactDuration(minMass,mNormalSpringConstant,mNormalBeta);
}

template<typename Dimension>
typename LinearSpringDEM<Dimension>::TimeStepType
LinearSpringDEM<Dimension>::
variableTimeStep(const DataBase<Dimension>& dataBase,
                 const State<Dimension>& state,
                 const StateDerivatives<Dimension>& /*derivs*/,
                 const typename Dimension::Scalar /*currentTime*/) const {
  
  // FieldLists we'll need for the timestep calc.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto r = state.fields(DEMFieldNames::particleRadius, 0.0);

  const auto& contacts = this->contactStorageIndices();
  const unsigned int numP2PContacts = this->numParticleParticleContacts();
  const unsigned int numTotContacts = this->numContacts();

  // vec of ptrs to our solid boundary conditions
  const auto& solidBoundaries = this->solidBoundaryConditions();

  // buffer distance used to set the max allowable timestep
  const auto bufferDistance = dataBase.maxNeighborSearchBuffer();

  // Compute the spring timestep constraint (except for the mass)
  const auto nsteps = this->stepsPerCollision();
  const auto dtSpring0 = computeContactDuration(1.0,mNormalSpringConstant,mNormalBeta)/nsteps;

  CHECK(nsteps > 0);

  // Check each pair to see how much they are closing.
  auto dtMin = std::numeric_limits<Scalar>::max();
  TimeStepType result(dtMin, "DEM error, this message should not get to the end");

  // start with the particle to particle contacts
  for (auto kk = 0u; kk < numP2PContacts; ++kk) {

    const auto nodeListi = contacts[kk].storeNodeList;
    const auto nodeListj = contacts[kk].pairNodeList;
    const auto  i = contacts[kk].storeNode;
    const auto  j = contacts[kk].pairNode;

    // node i
    const auto  mi = mass(nodeListi, i);
    const auto& ri = position(nodeListi, i);
    const auto& vi = velocity(nodeListi, i);
    const auto  Ri = r(nodeListi, i);

    // node j
    const auto  mj = mass(nodeListj, j);
    const auto& rj = position(nodeListj, j);
    const auto& vj = velocity(nodeListj, j);
    const auto  Rj = r(nodeListj, j);

    // Compare closing speed to separation
    const auto rji = rj - ri;
    const auto vji = vj - vi;
    const auto rhatji = rji.unitVector();
    const auto rmagji = rji.magnitude();
    const auto closing_speed = -min(0.0, vji.dot(rhatji));
    const auto overlap = (Ri + Rj) - rmagji;

    if (closing_speed > 0.0 or overlap > 0.0) {
      const auto dtSpringkk = dtSpring0*std::sqrt(std::min(mi, mj));
      const auto dtji = ( overlap > 0.0 ?
                          dtSpringkk :
                          max(dtSpringkk, -0.5*overlap*safeInvVar(closing_speed)));
      if (dtji < dtMin) {
        dtMin = dtji;
        result = make_pair(dtji,
                           (dtji > dtSpringkk ?
                            std::string("DEM pairwise closing rate:\n") :
                            std::string("DEM spring constraint:\n")) +
                           "  (listi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                           "  (listj, j, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(Process::getRank()) + ")\n" +
                           "        position_i = " + vec_to_string(ri) + "\n" +
                           "        velocity_i = " + vec_to_string(vi) +"\n" +
                           "          radius_i = " + to_string(Ri) + "\n" +
                           "        position_j = " + vec_to_string(rj) + "\n" +
                           "        velocity_j = " + vec_to_string(vj) + "\n" +
                           "          radius_j = " + to_string(Rj) + "\n" +
                           "         spring_dt = " + to_string(dtSpringkk) + "\n");
      }
    }
  }

  // particle to boundary contacts
  for (auto kk = numP2PContacts; kk < numTotContacts; ++kk) {

    const auto nodeListi = contacts[kk].storeNodeList;
    const auto i = contacts[kk].storeNode;
    const auto bci = contacts[kk].solidBoundary;

    // node i
    const auto  mi = mass(nodeListi, i);
    const auto& ri = position(nodeListi, i);
    const auto& vi = velocity(nodeListi, i);
    const auto  Ri = r(nodeListi, i);

    //solid boundary and distance vector to particle i
    const auto& solidBoundary = solidBoundaries[bci];
    const auto rib = solidBoundary->distance(ri);
    const auto vb = solidBoundary->localVelocity(ri);

    // Compare closing speed to separation
    const auto vib = vi-vb;
    const auto rhatib = rib.unitVector();
    const auto rmagib = rib.magnitude();

    const auto closing_speed = -min(0.0, vib.dot(rhatib));
    const auto overlap =  Ri - rmagib;
    
    if (closing_speed > 0.0 or overlap > 0.0) {

      const auto dtSpringkk = dtSpring0*std::sqrt(mi);
      // Spring constant timestep for this pair
      const auto dtji = (overlap > 0.0 ?
                         dtSpringkk:
                         max(dtSpringkk, -0.5*overlap*safeInvVar(closing_speed)));
      if (dtji < dtMin) {
        dtMin = dtji;
        result = make_pair(dtji,
                           (dtji > dtSpringkk ?
                            std::string("DEM pairwise closing rate:\n") :
                            std::string("DEM spring constraint:\n")) +
                           "  (listi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                           "        position_i = " + vec_to_string(ri) + "\n" +
                           "        velocity_i = " + vec_to_string(vi) +"\n" +
                           "          radius_i = " + to_string(Ri) + "\n" +
                           "    distance to bc = " + vec_to_string(rib) + "\n" +
                           "    velocity of bc = " + vec_to_string(vb) + "\n" +
                           "         spring_dt = " + to_string(dtSpringkk) + "\n");
      }
    }
  }


  // Ensure no point moves further than the buffer distance in one timestep
  //--------------------------------------------------------------------------------------
  // NOTE: it would be nice if this wasn't based on the absolute velocity for cases
  //       where we have a blob of dem particles moving at elevated speeds
  //--------------------------------------------------------------------------------------
  const auto numNodeLists = position.size();
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = position[k]->size();
    for (auto i = 0u; i < n; ++i) {
      const auto dti = 0.5* bufferDistance *r(k,i)*safeInvVar(velocity(k,i).magnitude());
      if (dti < dtMin) {
        dtMin = dti;
        result = make_pair(dti,
                           std::string("DEM single particle velocity limit:\n") +
                           "  (listi, i, rank) = (" + to_string(k) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                           "        position_i = " + vec_to_string(position(k,i)) + "\n" +
                           "        velocity_i = " + vec_to_string(velocity(k,i)) +"\n" +
                           "          radius_i = " + to_string(r(k,i)) + "\n");
      }
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// method that fires once on startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs){
  TIME_BEGIN("preStepInitialize");
  DEMBase<Dimension>::preStepInitialize(dataBase,state,derivs);
  this->recomputeContactDuration();
  TIME_END("preStepInitialize");
}

//------------------------------------------------------------------------------
// method that fires once on startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase){
  TIME_BEGIN("LinearSpringDEMinitializeProblemStartup");
  DEMBase<Dimension>::initializeProblemStartup(dataBase);
  this->setMomentOfInertia();
  TIME_END("LinearSpringDEMinitializeProblemStartup");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("LinearSpringDEMregisterState");

  DEMBase<Dimension>::registerState(dataBase,state);

  dataBase.resizeDEMFieldList(mMomentOfInertia, 0.0, DEMFieldNames::momentOfInertia, false);
  dataBase.resizeDEMFieldList(mMaximumOverlap, 0.0, DEMFieldNames::maximumOverlap, false);

  auto maxOverlapPolicy = make_policy<MaxReplaceState<Dimension,Scalar>>();

  state.enroll(mMomentOfInertia);
  state.enroll(mMaximumOverlap, maxOverlapPolicy);

  TIME_END("LinearSpringDEMregisterState");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
              StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("LinearSpringDEMregisterDerivs");

  DEMBase<Dimension>::registerDerivatives(dataBase,derivs);

  dataBase.resizeDEMFieldList(mNewMaximumOverlap, 0.0, DEMFieldNames::maximumOverlap, false);

  derivs.enroll(mNewMaximumOverlap);

  TIME_END("LinearSpringDEMregisterDerivs");
}

//------------------------------------------------------------------------------
// evaluate the derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const{
  TIME_BEGIN("LinearSpringDEMevaluateDerivatives");                    
  this->resizeDerivativePairFieldLists(derivatives);

  // A few useful constants we'll use in the following loop.
  const double tiny = std::numeric_limits<double>::epsilon();
  const auto shapeFactor = this->shapeFactor();
  const auto shapeFactor2 = shapeFactor*shapeFactor; 

  const auto muD = this->dynamicFrictionCoefficient();
  const auto muS = this->staticFrictionCoefficient();
  const auto muT = this->torsionalFrictionCoefficient() * shapeFactor * muS;
  const auto muR = this->rollingFrictionCoefficient() * shapeFactor;

  const auto Cc = this->cohesiveTensileStrength();

  // spring constants
  const auto kn = this->normalSpringConstant();        // normal
  const auto ks = this->tangentialSpringConstant();    // sliding
  const auto kt = 2.0 * ks * shapeFactor2;
  const auto kr = kn * shapeFactor2;

  const auto invKs = 1.0/max(ks,tiny);
  const auto invKt = 1.0/max(kt,tiny);
  const auto invKr = 1.0/max(kr,tiny);
  
  const auto normalDampingTerms = 4.0*kn/(1.0+mNormalBeta*mNormalBeta);
  const auto tangentialDampingTerms = 4.0*ks/(1.0+mTangentialBeta*mTangentialBeta);
 
  // The connectivity.
  const auto& nodeLists = dataBase.DEMNodeListPtrs();
  const auto  numNodeLists = nodeLists.size();
  
  // contact contacts for different types
  const unsigned int numP2PContacts = this->numParticleParticleContacts();
  //const unsigned int numP2BContacts = this->numParticleBoundaryContacts();
  const unsigned int numTotContacts = this->numContacts();

  // Get the state FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto momentOfInertia = state.fields(DEMFieldNames::momentOfInertia, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto omega = state.fields(DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  const auto radius = state.fields(DEMFieldNames::particleRadius, 0.0);
  const auto uniqueIndices = state.fields(DEMFieldNames::uniqueIndices, (size_t)0);
  const auto compositeIndex = state.fields(DEMFieldNames::compositeParticleIndex, (int)0);
  
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(radius.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(uniqueIndices.size() == numNodeLists);
  CHECK(compositeIndex.size() == numNodeLists);

  // Get the state pairFieldLists.
  const auto equilibriumOverlap = state.fields(DEMFieldNames::equilibriumOverlap, std::vector<Scalar>());
  const auto shearDisplacement = state.fields(DEMFieldNames::shearDisplacement, std::vector<Vector>());
  const auto rollingDisplacement = state.fields(DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  const auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  const auto neighborIds = state.fields(DEMFieldNames::neighborIndices, std::vector<size_t>());

  CHECK(equilibriumOverlap.size() == numNodeLists);
  CHECK(shearDisplacement.size() == numNodeLists);
  CHECK(rollingDisplacement.size() == numNodeLists);
  CHECK(torsionalDisplacement.size() == numNodeLists);
  CHECK(neighborIds.size() == numNodeLists);

  // Get the deriv FieldLists
  auto DxDt = derivatives.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto DomegaDt = derivatives.fields(IncrementState<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  auto newMaximumOverlap = derivatives.fields(MaxReplaceState<Dimension,Scalar>::prefix() + DEMFieldNames::maximumOverlap, 0.0);
  
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DomegaDt.size() == numNodeLists);
  CHECK(newMaximumOverlap.size() == numNodeLists);
  
  // Get the deriv pairFieldLists
  auto DDtShearDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto newShearDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto DDtRollingDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto newRollingDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto DDtTorsionalDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::incrementPrefix() +  DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  auto newTorsionalDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::replacePrefix() +  DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());

  CHECK(DDtShearDisplacement.size() == numNodeLists);
  CHECK(newShearDisplacement.size() == numNodeLists);
  CHECK(DDtRollingDisplacement.size() == numNodeLists);
  CHECK(newRollingDisplacement.size() == numNodeLists);
  CHECK(DDtTorsionalDisplacement.size() == numNodeLists);
  CHECK(newTorsionalDisplacement.size() == numNodeLists);

  // storage locations in our pairwise FieldLists
  const auto& contacts = this->contactStorageIndices();

  // vec of ptrs to our solid boundary conditions
  const auto& solidBoundaries = this->solidBoundaryConditions();

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj, contacti,bci;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DomegaDt_thread = DomegaDt.threadCopy(threadStack);
    auto newMaxOverlap_thread = newMaximumOverlap.threadCopy(threadStack, ThreadReduction::MAX);

    //------------------------------------
    // particle-particle contacts
    //------------------------------------
#pragma omp for
    for (auto kk = 0u; kk < numP2PContacts; ++kk) {
      
      CHECK(contacts[kk].storeNodeList >= 0)
      CHECK(contacts[kk].storeNode >= 0)
      CHECK(contacts[kk].pairNodeList >= 0)
      CHECK(contacts[kk].pairNode >= 0)
      CHECK(contacts[kk].solidBoundary == -1)

      nodeListi = contacts[kk].storeNodeList;
      nodeListj = contacts[kk].pairNodeList;
      i = contacts[kk].storeNode;
      j = contacts[kk].pairNode;
      contacti = contacts[kk].storeContact;

      // Get the state vars for node i needed for prox check
      const auto& ri = position(nodeListi, i);
      const auto  Ri = radius(nodeListi, i);
      CHECK(Ri > 0.0);

      // Get the state vars for node j needed for prox check
      const auto& rj = position(nodeListj, j);
      const auto  Rj = radius(nodeListj, j);
      CHECK(Rj > 0.0);

      // are we overlapping ? 
      const auto rij = ri-rj;
      const auto rijMag = rij.magnitude();
      const auto delta0 = (Ri+Rj) - rijMag;  // raw delta
      CHECK2(rijMag > 0.0, "ERROR: particles (" << nodeListi << " " << i << ") and (" << nodeListj << " " << j << ") share the same position @ " << ri);
      
      // if so do the things
      if (delta0 > 0.0){
        
        // get remaining state for node i
        const auto  cIdi = compositeIndex(nodeListi,i);
        const auto  uIdi = uniqueIndices(nodeListi,i);
        const auto& vi = velocity(nodeListi, i);
        const auto& omegai = omega(nodeListi, i);
        const auto  mi = mass(nodeListi, i);
        const auto  Ii = momentOfInertia(nodeListi,i);
        CHECK(mi > 0.0);
        CHECK(Ii > 0.0);

        // get remaining state for node j
        const auto  cIdj = compositeIndex(nodeListj,j);
        const auto  uIdj = uniqueIndices(nodeListj,j);
        const auto& vj = velocity(nodeListj, j);
        const auto& omegaj = omega(nodeListj, j);
        const auto  mj = mass(nodeListj, j);
        const auto  Ij = momentOfInertia(nodeListj,j);
        CHECK(mj > 0.0);
        CHECK(Ij > 0.0);

        // Get the derivs from node i
        auto& DvDti = DvDt_thread(nodeListi, i);
        auto& DomegaDti = DomegaDt_thread(nodeListi, i);
        auto& maxOverlapi = newMaxOverlap_thread(nodeListi,i);

        // Get the derivs from node j
        auto& DvDtj = DvDt_thread(nodeListj, j);
        auto& DomegaDtj = DomegaDt_thread(nodeListj, j);
        auto& maxOverlapj = newMaxOverlap_thread(nodeListj,j);

        // storage sign, this makes pairwise values i-j independent
        const int storageSign = (uIdi <= uIdj ? 1 : -1);

        // stored pair-wise values
        const auto  overlapij   = equilibriumOverlap(nodeListi,i)[contacti];
        const auto& deltaSlidij = shearDisplacement(nodeListi,i)[contacti]*storageSign;
        const auto& deltaRollij = rollingDisplacement(nodeListi,i)[contacti];
        const auto  deltaTorsij = torsionalDisplacement(nodeListi,i)[contacti];

        // boolean checks
        const auto isBondedParticle = (cIdi == cIdj);
        const auto allowSliding = (!isBondedParticle);

        // effective delta
        const auto delta = delta0 - overlapij;

        // line of action for the contact
        const auto rhatij = rij.unitVector();

        // reduced radii
        const auto li = (Ri*Ri-Rj*Rj + rijMag*rijMag)*safeInv(2.0*rijMag);
        const auto lj = rijMag-li;

        // effective quantities (mass, reduced radius) respectively, mij->mi for mi=mj
        const auto mij = mi*mj/(mi+mj);
        const auto lij = li*lj/(li+lj);

        // damping constants -- ct and cr derived quantities ala Zhang 2017
        const auto Cn = std::sqrt(mij*normalDampingTerms);
        const auto Cs = std::sqrt(mij*tangentialDampingTerms);
        const auto Ct = 2.0 * Cs * shapeFactor2;
        const auto Cr = Cn * shapeFactor2;

        // our velocities
        const Vector vroti = -DEMDimension<Dimension>::cross(omegai,rhatij);
        const Vector vrotj =  DEMDimension<Dimension>::cross(omegaj,rhatij);

        const Vector vij = vi-vj + li*vroti - lj*vrotj;

        const Scalar vn = vij.dot(rhatij);                   // normal velocity
        const Vector vs = vij - vn*rhatij;                   // sliding velocity
        const Vector vr = -li*vroti - lj*vrotj;              // rolling velocity
        const Scalar vt = -lij*DEMDimension<Dimension>::dot(omegai-omegaj,rhatij); // torsion velocity

        // normal forces 
        //------------------------------------------------------------
        const Vector fn = (kn*delta - Cn*vn)*rhatij;          // normal spring
        const Vector fc = 4.0*Cc*shapeFactor2*lij*lij*rhatij; // normal cohesion
        const Scalar fnMag = fn.magnitude();                  // magnitude of normal spring force

        // sliding
        //------------------------------------------------------------
        Vector newDeltaSlidij, fs;
        this->slidingSpringDamper(ks,Cs,muS,muD,deltaSlidij,vs,fnMag,invKs,rhatij,allowSliding,
                                  newDeltaSlidij,fs); // outputs

        // torsion
        //------------------------------------------------------------
        Scalar newDeltaTorsij, ft;
        this->slidingSpringDamper(kt,Ct,muT,muT,deltaTorsij,vt,fnMag,invKt,allowSliding,
                                  newDeltaTorsij,ft); // output

        // rolling
        //------------------------------------------------------------
        Vector newDeltaRollij, fr;
        this->slidingSpringDamper(kr,Cr,muR,muR,deltaRollij,vr,fnMag,invKr,rhatij,allowSliding,
                                  newDeltaRollij,fr); // outputs
        // accelerations
        //------------------------------------------------------------
        // Rectilinear Acceleration 
        const auto fij = fn - fc + fs;
        DvDti += fij/mi;
        DvDtj -= fij/mj;

        // angular acceleration
        const auto Msliding = -DEMDimension<Dimension>::cross(rhatij,fij);
        const auto Mrolling = -DEMDimension<Dimension>::cross(rhatij,fr);
        const auto Mtorsion = ft * this->torsionMoment(rhatij,omegai,omegaj); // rename torsionDirection
        DomegaDti += (Msliding*li - (Mtorsion + Mrolling) * lij)/Ii;
        DomegaDtj += (Msliding*lj + (Mtorsion + Mrolling) * lij)/Ij;

        // update max overlaps
        maxOverlapi = max(maxOverlapi,delta);
        maxOverlapj = max(maxOverlapj,delta);

        // for spring updates
        newShearDisplacement(nodeListi,i)[contacti] = storageSign*newDeltaSlidij;
        DDtShearDisplacement(nodeListi,i)[contacti] = storageSign*vs;
        newTorsionalDisplacement(nodeListi,i)[contacti] = newDeltaTorsij;
        DDtTorsionalDisplacement(nodeListi,i)[contacti] = vt;
        newRollingDisplacement(nodeListi,i)[contacti] = newDeltaRollij;
        DDtRollingDisplacement(nodeListi,i)[contacti] = vr;
    
      } // if contacing
    }   // loop over pairs


    //?????????????????????????????????????????????????????????????????//
    //------------------------------------
    // Loop the particle-boundary contacts
    //------------------------------------
    //?????????????????????????????????????????????????????????????????//
    #pragma omp for
    for (auto kk = numP2PContacts; kk < numTotContacts; ++kk) {
    
      CHECK(contacts[kk].storeNodeList >= 0)
      CHECK(contacts[kk].storeNode >= 0)
      CHECK(contacts[kk].pairNodeList == -1)
      CHECK(contacts[kk].pairNode == -1)
      CHECK(contacts[kk].solidBoundary >= 0)

      nodeListi = contacts[kk].storeNodeList;
      i = contacts[kk].storeNode;
      contacti = contacts[kk].storeContact;
      bci = contacts[kk].solidBoundary;

      //Get the state vars for node i needed for prox check
      const auto& ri = position(nodeListi, i);
      const auto& Ri = radius(nodeListi, i);
      CHECK(Ri > 0.0);

      //solid boundary and distance vector to particle i
      const auto& solidBoundary = solidBoundaries[bci];
      const auto rib = solidBoundary->distance(ri);

      //effective delta
      const auto delta = Ri-rib.magnitude();

      // are overlapping?
      if (delta > 0.0){

        // get remaining state for node i
        const auto& vi = velocity(nodeListi, i);
        const auto& omegai = omega(nodeListi, i);
        const auto& mi = mass(nodeListi, i);
        const auto& Ii = momentOfInertia(nodeListi, i);
        CHECK(Ii > 0.0);
        CHECK(mi > 0.0);

        // Get the derivs from node i
        auto& DvDti = DvDt_thread(nodeListi, i);
        auto& DomegaDti = DomegaDt_thread(nodeListi, i);
        auto& maxOverlapi = newMaxOverlap_thread(nodeListi,i);

        // get pairwise variables
        const auto& deltaSlidib = shearDisplacement(nodeListi,i)[contacti];
        const auto& deltaRollib = rollingDisplacement(nodeListi,i)[contacti];
        const auto  deltaTorsib = torsionalDisplacement(nodeListi,i)[contacti];

        // velocity of boundary @ ri
        const auto vb = solidBoundary->localVelocity(ri);

        // line of action for the contact
        const auto rhatib = rib.unitVector();

        // lever arm 
        const auto li = rib.magnitude();

        // effective quantities
        const auto mib = mi;
        const auto lib = li;

        // damping coefficients
        const auto Cn = std::sqrt(mib*normalDampingTerms);
        const auto Cs = std::sqrt(mib*tangentialDampingTerms);
        const auto Ct = 2.0 * Cs * shapeFactor2;
        const auto Cr = Cn * shapeFactor2;

        // velocities
        const Vector vroti = -DEMDimension<Dimension>::cross(omegai,rhatib);

        const Vector vib = vi-vb + li*vroti;

        const Scalar vn = vib.dot(rhatib);                                  // normal velocity
        const Vector vs = vib - vn*rhatib;                                  // sliding velocity
        const Vector vr = -li*vroti;                                        // rolling velocity
        const Scalar vt = -lib*DEMDimension<Dimension>::dot(omegai,rhatib); // torsion velocity

        // normal forces 
        //------------------------------------------------------------
        const Vector fn = (kn*delta - Cn*vn)*rhatib;          // normal spring
        const Vector fc = 4.0*Cc*shapeFactor2*lib*lib*rhatib; // normal cohesion
        const Scalar fnMag = fn.magnitude();                  // magnitude of normal spring force

        // sliding
        //------------------------------------------------------------
        Vector newDeltaSlidib, ft;
        this->slidingSpringDamper(ks,Cs,muS,muD,deltaSlidib,vs,fnMag,invKs,rhatib,true,
                                  newDeltaSlidib,ft);  // outputs

        // torsion
        //------------------------------------------------------------
        Scalar newDeltaTorsib, MtorsionMag;
        this->slidingSpringDamper(kt,Ct,muT,muT,deltaTorsib,vt,fnMag,invKt,true,
                                  newDeltaTorsib, MtorsionMag); // output

        // rolling
        //------------------------------------------------------------
        Vector newDeltaRollib, fr;
        this->slidingSpringDamper(kr,Cr,muR,muR,deltaRollib,vr,fnMag,invKr,rhatib,true,
                                  newDeltaRollib, fr); // outputs
        // accelerations
        //------------------------------------------------------------
        // Rectilinear Acceleration 
        const Vector fib = fn - fc + ft;
        DvDti += fib/mi;

        // angular acceleration
        const auto Msliding = -DEMDimension<Dimension>::cross(rhatib,fib);
        const auto Mrolling = -DEMDimension<Dimension>::cross(rhatib,fr);
        const auto Mtorsion = MtorsionMag * this->torsionMoment(rhatib,omegai,0*omegai); // rename torsionDirection
        DomegaDti += (Msliding*li - (Mtorsion + Mrolling) * lib)/Ii;

        // update max overlaps
        maxOverlapi = max(maxOverlapi,delta);

        // for spring updates
        newShearDisplacement(nodeListi,i)[contacti] = newDeltaSlidib;
        DDtShearDisplacement(nodeListi,i)[contacti] = vs;
        newTorsionalDisplacement(nodeListi,i)[contacti] = newDeltaTorsib;
        DDtTorsionalDisplacement(nodeListi,i)[contacti] = vt;
        newRollingDisplacement(nodeListi,i)[contacti] = newDeltaRollib;
        DDtRollingDisplacement(nodeListi,i)[contacti] = vr;

     } // if statement
   }   // loop pairs
  threadReduceFieldLists<Dimension>(threadStack);
  }    // omp parfor

  // finish with loop over nodelists
  //-----------------------------------------------------------
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = nodeLists[nodeListi];
    const auto ni = nodeList->numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
        const auto veli = velocity(nodeListi,i);
        DxDt(nodeListi,i) = veli;
    }   // loop nodes
  }     // loop nodelists

  TIME_END("LinearSpringDEMevaluateDerivatives");
}       // method


//------------------------------------------------------------------------------
// method that fires once on startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
setMomentOfInertia() {

  const auto& db = this->dataBase();
  const auto  mass   = db.DEMMass();
  const auto  radius = db.DEMParticleRadius();
  const auto  numNodeLists = db.numNodeLists();
  const auto& nodeLists = db.DEMNodeListPtrs();

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = nodeLists[nodeListi];
    const auto ni = nodeList->numInternalNodes();

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
        const auto mi = mass(nodeListi,i);
        const auto Ri = radius(nodeListi,i);
        mMomentOfInertia(nodeListi,i) = this->momentOfInertia(mi,Ri); 

    }   // loop nodes
  }     // loop nodelists
}       // method



//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for linear dem state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  DEMBase<Dimension>::applyGhostBoundaries(state,derivs);
  auto I = state.fields(DEMFieldNames::momentOfInertia,0.0);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(I);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for linear dem state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  DEMBase<Dimension>::enforceBoundaries(state,derivs);

  auto I = state.fields(DEMFieldNames::momentOfInertia,0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(I);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DEMBase<Dimension>::dumpState(file,pathName);
  file.write(mMomentOfInertia, pathName + "/momentOfInertia");
  file.write(mMaximumOverlap, pathName + "/maximumOverlap");
  file.write(mNewMaximumOverlap, pathName + "/newMaximumOverlap");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DEMBase<Dimension>::restoreState(file,pathName);
  file.read(mMomentOfInertia, pathName + "/momentOfInertia");
  file.read(mMaximumOverlap, pathName + "/maximumOverlap");
  file.read(mNewMaximumOverlap, pathName + "/newMaximumOverlap");
}
} // namespace

