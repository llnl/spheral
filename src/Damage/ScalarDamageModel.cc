//---------------------------------Spheral++----------------------------------//
// ScalarDamageModel -- Base class for the scalar damage physics models.
// This class does not know how to seed the flaw distribution -- that is 
// required of descendent classes.
//
// Created by JMO, Sun Oct 10 17:22:05 PDT 2004
//----------------------------------------------------------------------------//
#include "ScalarDamageModel.hh"
#include "SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "StrainPolicy.hh"
#include "ScalarDamagePolicy.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"
#include "Field/FieldList.hh"

#include <string>

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ScalarDamageModel<Dimension>::
ScalarDamageModel(SolidNodeList<Dimension>& nodeList,
                  FluidNodeList<Dimension>& damagedNodeList,
                  const double kernelExtent,
                  const double crackGrowthMultiplier,
                  const FlawStorageType& flaws):
  DamageModel<Dimension>(nodeList, kernelExtent, crackGrowthMultiplier, flaws),
  mDamagedNodeListPtr(&damagedNodeList),
  mSplitFailedNodes(true),
  mStrain(SolidFieldNames::strain, nodeList),
  mDamage(SolidFieldNames::scalarDamage, nodeList),
  mDdamageDt(ScalarDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mMass0(HydroFieldNames::mass + "0", nodeList.mass()),
  mUndamagedToDamagedIndex("UD -> D index", nodeList, -1),
  mBoundPlanes() {

    // Make sure the damaged NodeList is initially empty.
    damagedNodeList.numInternalNodes(0);
    damagedNodeList.numGhostNodes(0);

    ENSURE(mDamagedNodeListPtr->numInternalNodes() == 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ScalarDamageModel<Dimension>::
~ScalarDamageModel() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the scalar damage field.
// We set the total damage derivative assuming nothing about what cracks are
// activated -- this is scaled to the correct value for the active cracks 
// in the damage update policy.
// Note also that this is actually the derivative of the 1/3 power of the 
// damage.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamageModel<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Fields we're going to set.
  typedef typename State<Dimension>::FieldKeyType KeyType;
  const KeyType DdamageDtKey(&(this->nodeList()),
                             ScalarDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage);
  Field<Dimension, Scalar>& DDDt = derivatives.scalarField(DdamageDtKey);
  
  // DamageModel actually knows how to do this.
  this->computeScalarDDDt(dataBase,
                          state,
                          time,
                          dt,
                          DDDt);

}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename ScalarDamageModel<Dimension>::TimeStepType
ScalarDamageModel<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {

//   // Look at how quickly we're trying to change the damage.
//   double dt = DBL_MAX;
//   for (int i = 0; i != this->nodeList().numInternalNodes(); ++i) {
//     dt = min(dt, 0.2*max(mDamage(i), 1.0 - mDamage(i))/
//              std::sqrt(mDdamageDt(i)*mDdamageDt(i) + 1.0e-20));
//   }

  return TimeStepType(1.0e100, "Rate of damage change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::FieldKeyType KeyType;
  typedef typename State<Dimension>::ScalarPolicyPointerType ScalarPolicyPointer;
  typedef typename State<Dimension>::SymTensorPolicyPointerType SymTensorPolicyPointer;

  // Register the strain.
  ScalarPolicyPointer strainPolicy(new StrainPolicy<Dimension>());
  state.registerField(mStrain, strainPolicy);

  // Register the damage field.
  ScalarPolicyPointer damagePolicy(new ScalarDamagePolicy<Dimension>(*this));
  state.registerField(mDamage, damagePolicy);

  // Get the current set of position weights.
  FieldList<Dimension, Scalar> positionWeight = state.scalarFields(HydroFieldNames::positionWeight);
  const SolidNodeList<Dimension>& nodeList = this->nodeList();
  CHECK(positionWeight.fieldForNodeList(nodeList) < positionWeight.end());
  CHECK(positionWeight.fieldForNodeList(*mDamagedNodeListPtr) < positionWeight.end());
  Field<Dimension, Scalar>& undamagedPositionWeight = **(positionWeight.fieldForNodeList(nodeList));
  Field<Dimension, Scalar>& damagedPositionWeight = **(positionWeight.fieldForNodeList(*mDamagedNodeListPtr));

  // Scale up the position weight of the damaged nodes.  This is a crude method
  // for causing the h's of surrounding nodes to decrease.
  damagedPositionWeight = 5.0;

  // If there are damaged nodes overlapping undamaged, we have to modify the 
  // position weights for the H calculation.
  if (mUndamagedToDamagedIndex.size() > 0) {
    const SolidNodeList<Dimension>& nodeList = this->nodeList();
    for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
      const int j = mUndamagedToDamagedIndex(i);
      if (j >= 0) {
        CHECK(j < mDamagedNodeListPtr->numInternalNodes());
        undamagedPositionWeight(i) = 0.5;
        damagedPositionWeight(j) = 0.5;
      }
    }

  }

  // Register the base classes stuff.
  DamageModel<Dimension>::registerState(dataBase, state);

}

//------------------------------------------------------------------------------
// Register the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamageModel<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  // Register the damage time derivative.
  derivs.registerField(mDdamageDt);
}

//------------------------------------------------------------------------------
// Finalize method, called at the end of a time step.
// We use this opportunity to create mass in the "damaged" material, removing
// it from the undamaged in proportion to the scalar damage measure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamageModel<Dimension>::
finalize(const typename Dimension::Scalar time, 
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& db, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  const double tiny = 1.0e-20;

  // Get the state fields we care about.
  typedef typename State<Dimension>::FieldKeyType KeyType;
  SolidNodeList<Dimension>& nodeList = this->nodeList();
  const KeyType umKey(&nodeList, HydroFieldNames::mass);
  const KeyType urKey(&nodeList, HydroFieldNames::position);
  const KeyType uvKey(&nodeList, HydroFieldNames::velocity);
  const KeyType uHKey(&nodeList, HydroFieldNames::H);
  const KeyType uRhoKey(&nodeList, HydroFieldNames::massDensity);
  const KeyType uuKey(&nodeList, HydroFieldNames::specificThermalEnergy);
  const KeyType uwKey(&nodeList, HydroFieldNames::weight);
  const KeyType DKey(&nodeList, SolidFieldNames::scalarDamage);
  const KeyType PKey(&nodeList, HydroFieldNames::pressure);
  const KeyType SKey(&nodeList, SolidFieldNames::deviatoricStress);
  const KeyType clKey(&nodeList, SolidFieldNames::longitudinalSoundSpeed);

  const KeyType dmKey(mDamagedNodeListPtr, HydroFieldNames::mass);
  const KeyType drKey(mDamagedNodeListPtr, HydroFieldNames::position);
  const KeyType dvKey(mDamagedNodeListPtr, HydroFieldNames::velocity);
  const KeyType dHKey(mDamagedNodeListPtr, HydroFieldNames::H);
  const KeyType dRhoKey(mDamagedNodeListPtr, HydroFieldNames::massDensity);
  const KeyType duKey(mDamagedNodeListPtr, HydroFieldNames::specificThermalEnergy);
  const KeyType dwKey(mDamagedNodeListPtr, HydroFieldNames::weight);

  Field<Dimension, Scalar>& um = state.scalarField(umKey);
  Field<Dimension, Vector>& ur = state.vectorField(urKey);
  Field<Dimension, Vector>& uv = state.vectorField(uvKey);
  Field<Dimension, SymTensor>& uH = state.symTensorField(uHKey);
  Field<Dimension, Scalar>& urho = state.scalarField(uRhoKey);
  Field<Dimension, Scalar>& uu = state.scalarField(uuKey);
  Field<Dimension, Scalar>& uw = state.scalarField(uwKey);
  Field<Dimension, Scalar>& D = state.scalarField(DKey);
  Field<Dimension, Scalar>& P = state.scalarField(PKey);
  Field<Dimension, SymTensor>& S = state.symTensorField(SKey);
  Field<Dimension, Scalar>&  cl = state.scalarField(clKey);

  Field<Dimension, Scalar>& dm = state.scalarField(dmKey);
  Field<Dimension, Vector>& dr = state.vectorField(drKey);
  Field<Dimension, Vector>& dv = state.vectorField(dvKey);
  Field<Dimension, SymTensor>& dH = state.symTensorField(dHKey);
  Field<Dimension, Scalar>& drho = state.scalarField(dRhoKey);
  Field<Dimension, Scalar>& du = state.scalarField(duKey);
  Field<Dimension, Scalar>& dw = state.scalarField(dwKey);

  // Iterate over all the (internal) nodes in the undamaged material.
  vector<int> nodesToKill;
  vector<int> nodesToSplit;
  for (int i = 0; i != nodeList.numInternalNodes(); ++i) {

    // Is this node damaged?
    if (D(i) > 0.0) {
      const Scalar f = 1.0 - D(i);
      CHECK(f >= 0.0 && f <= 1.0);

      // Get the damaged node index corresponding to this undamaged node.
      // If there is no damaged node yet, then create one.
      int j = mUndamagedToDamagedIndex(i);
      if (j < 0) {
        j = mDamagedNodeListPtr->numInternalNodes();
        mDamagedNodeListPtr->numInternalNodes(j + 1);
        CHECK(j >= 0 && j < mDamagedNodeListPtr->numInternalNodes());
        mUndamagedToDamagedIndex(i) = j;
        dm(j) = 0.0;
        dr(j) = ur(i);
        dv(j) = uv(i);
        dH(j) = uH(i);
        drho(j) = urho(i);
        du(j) = uu(i);
      }
      CHECK(j >= 0 && j < mDamagedNodeListPtr->numInternalNodes());
      CHECK(mUndamagedToDamagedIndex(i) == j);

      // Split the mass up between the damaged and undamaged nodes.
      const Scalar um0 = um(i);
      const Scalar dm0 = dm(j);
      um(i) = mMass0(i)*f;
      dm(j) = mMass0(i)*(1.0 - f);
      CHECK(fuzzyEqual(um0 + dm0, mMass0(i)));
      CHECK(fuzzyEqual(um(i) + dm(j), mMass0(i)));
      CHECK(fuzzyLessThanOrEqual(um(i), um0));
      CHECK(fuzzyGreaterThanOrEqual(dm(j), dm0));

      // Assign the positions to the center of mass value.
      const double thpt = 1.0/mMass0(i);
      const Vector rcm = thpt*(um0*ur(i) + dm0*dr(j));
      ur(i) = rcm;
      dr(j) = rcm;

      // Combine the velocities, preserving the linear momentum.
      const Vector vcm = thpt*(um0*uv(i) + dm0*dv(j));
      uv(i) = vcm;
      dv(j) = vcm;

      // Combine the H's.
      const SymTensor uH0 = uH(i);
      const SymTensor dH0 = dH(j);
      const SymTensor Hcm = thpt*(um0*uH(i) + dm0*dH(j));
      uH(i) = Hcm;
      dH(j) = Hcm;

      // Divvy up the mass density.
      const Scalar rhocm = thpt*(um0*urho(i) + dm0*drho(j));
      urho(i) = rhocm;
      drho(j) = rhocm;

      // Split up the thermal energy, conserving the total.
//       const Scalar ucm = thpt*(um0*uu(i) + dm0*du(j));
//       uu(i) = ucm;
//       du(j) = ucm;
      const Scalar deltau = max(0.0, um0 - um(i))*uu(i);
      CHECK(deltau*uu(i) >= 0.0);
      uu(i) = (uu(i)*um0 - deltau)*um(i)/(um(i)*um(i) + tiny);
      du(j) = (du(j)*dm0 + deltau)*dm(j)/(dm(j)*dm(j) + tiny);

      // Correct the node weights.
      uw(i) = um(i)*urho(i)/(urho(i)*urho(i) + tiny);
      dw(j) = dm(j)*drho(j)/(drho(j)*drho(j) + tiny);

      // Should we delete the undamaged node?
      if (fuzzyEqual(D(i), 1.0)) {
        nodesToKill.push_back(i);
        nodesToSplit.push_back(j);
      }
    }
  }

  // If any nodes became fully damaged, we delete the original node
  // entirely.
  // Optionally at this point the damaged node can be split into two
  // half-mass nodes, in a simple attempt at dynamic refinement near
  // cracks.
  CHECK(nodesToKill.size() == nodesToSplit.size());
  const int num = nodesToKill.size();
  if (num > 0) {
    const double nPerh = mDamagedNodeListPtr->nodesPerSmoothingScale();
    CHECK(distinctlyGreaterThan(nPerh, 0.0));

    // Are we splitting the fully failed nodes?
    if (mSplitFailedNodes) {
      for (int k = 0; k != num; ++k) {
        CHECK(k < nodesToKill.size() && k < nodesToSplit.size());
        const int i = nodesToKill[k];
        CHECK(i >= 0 && i < nodeList.numInternalNodes());

        // Allocate a new node.
        mDamagedNodeListPtr->numInternalNodes(mDamagedNodeListPtr->numInternalNodes() + 1);
        const int j1 = nodesToSplit[k];
        const int j2 = mDamagedNodeListPtr->numInternalNodes() - 1;
        CHECK(j1 >= 0 && j1 < mDamagedNodeListPtr->numInternalNodes());
        CHECK(j2 >= 0 && j2 < mDamagedNodeListPtr->numInternalNodes());

        // Identify the direction we're going to fail in (the crack orientation).
        Vector strainDirection;
        {
          const SymTensor sigmai = S(i) - P(i)*SymTensor::one();
          const typename SymTensor::EigenStructType eigen = sigmai.eigenVectors();
          double maxValue = -DBL_MAX;
          int k = 0;
          while (k < Dimension::nDim) {
            if (eigen.eigenValues(k) > maxValue) {
              maxValue = eigen.eigenValues(k);
              strainDirection = eigen.eigenVectors.getColumn(k);
            }
            ++k;
          }
          CHECK(maxValue > -DBL_MAX);
        }
        CHECK(fuzzyEqual(strainDirection.magnitude2(), 1.0, 1.0e-10));

        // Split the mass up.
        const double m0 = dm(j1);
        dm(j1) *= 0.5;
        dm(j2) = dm(j1);
        CHECK(fuzzyEqual(dm(j1) + dm(j2), m0));

        // Assign the positions.
        const Vector rcm = dr(j1);
        const SymTensor Hi = dH(j1).Inverse();
        const double hi = (Hi*strainDirection).magnitude();
        const Vector dx = 0.5/nPerh*hi*strainDirection;
        dr(j1) = rcm + dx;
        dr(j2) = rcm - dx;

        // Check if either of the daughter nodes is out of bounds.  If so, set them
        // back to the center of mass.
        if (positionOutOfBounds(dr(j1)) || positionOutOfBounds(dr(j2))) {
          dr(j1) = rcm;
          dr(j2) = rcm;
        }

        // Assign the velocities such that the nodes are initially moving apart
        // at the crack growth speed.  We choose to align this relative motion with
        // the direction of maximum strain.
        const Vector vcm = dv(j1);
        //       const Vector vdelta = mCrackGrowthMultiplier*cl(i)*strainDirection;
        dv(j1) = vcm; //  + vdelta;
        dv(j2) = vcm; //  - vdelta;

        // Assign the H's.
        dH(j2) = dH(j1);

        // Assign the mass density.
        drho(j2) = drho(j1);

        // Split up the thermal energy, conserving the total.
        du(j2) = du(j1);

        // Correct the node weights.
        dw(j1) *= 0.5;
        dw(j2) = dw(j1);
      }
    }

    // Force the update of neighbor info on the damaged NodeList.
    mDamagedNodeListPtr->neighborPtr()->updateNodes();

    // Now remove any of the "undamaged" nodes that have become fully damaged.
    nodeList.deleteNodes(nodesToKill);
    nodeList.neighborPtr()->updateNodes();
  }

}

//------------------------------------------------------------------------------
// Access the internal parameters of the model.
//------------------------------------------------------------------------------
template<typename Dimension>
const FluidNodeList<Dimension>&
ScalarDamageModel<Dimension>::
damagedNodeList() const {
  CHECK(mDamagedNodeListPtr != 0);
  return *mDamagedNodeListPtr;
}

template<typename Dimension>
bool
ScalarDamageModel<Dimension>::
splitFailedNodes() const {
  return mSplitFailedNodes;
}

template<typename Dimension>
void
ScalarDamageModel<Dimension>::
splitFailedNodes(const bool x) {
  mSplitFailedNodes = x;
}

//------------------------------------------------------------------------------
// Access the state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
ScalarDamageModel<Dimension>::
strain() const {
  return mStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
ScalarDamageModel<Dimension>::
damage() const {
  return mDamage;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
ScalarDamageModel<Dimension>::
DdamageDt() const {
  return mDdamageDt;
}

//------------------------------------------------------------------------------
// Return the current undamaged to damaged NodeList indexing.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
ScalarDamageModel<Dimension>::
undamagedToDamagedNodeIndicies() const {
  vector<int> result;
  const SolidNodeList<Dimension>& nodeList = this->nodeList();
  result.reserve(nodeList.numInternalNodes());
  for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
    CHECK(mUndamagedToDamagedIndex(i) >= -1 &&
          mUndamagedToDamagedIndex(i) < mDamagedNodeListPtr->numInternalNodes());
    result.push_back(mUndamagedToDamagedIndex(i));
  }
  ENSURE(result.size() == nodeList.numInternalNodes());
  return result;
}

//------------------------------------------------------------------------------
// The current set of boundary planes.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<GeomPlane<Dimension> >&
ScalarDamageModel<Dimension>::
boundPlanes() {
  return mBoundPlanes;
}

//------------------------------------------------------------------------------
// Determine if the given position is out of the bounds as defined by a set of
// user defined planes.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ScalarDamageModel<Dimension>::
positionOutOfBounds(const Vector& r) const {
  bool result = false;
  typename vector<Plane>::const_iterator planeItr = mBoundPlanes.begin();
  while (!result && (planeItr != mBoundPlanes.end())) {
    const Vector dr = r - planeItr->point();
    const double thpt = dr.dot(planeItr->normal());
    result = (thpt <= 0.0);
    ++planeItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mStrain, pathName + "/strain");
  file.write(mDamage, pathName + "/damage");
  file.write(mDdamageDt, pathName + "/DdamageDt");
  file.write(mMass0, pathName + "/mass0");
  file.write(mUndamagedToDamagedIndex, pathName + "/undamagedToDamagedIndex");
  file.write(mBoundPlanes, pathName + "/boundPlanes");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DamageModel<Dimension>::restoreState(file, pathName);
  file.read(mStrain, pathName + "/strain");
  file.read(mDamage, pathName + "/damage");
  file.read(mDdamageDt, pathName + "/DdamageDt");
  file.read(mMass0, pathName + "/mass0");
  file.read(mUndamagedToDamagedIndex, pathName + "/undamagedToDamagedIndex");
  file.read(mBoundPlanes, pathName + "/boundPlanes");
}

}

//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class ScalarDamageModel<Dim<1> >;
  template class ScalarDamageModel<Dim<2> >;
  template class ScalarDamageModel<Dim<3> >;
}
