//---------------------------------Spheral++----------------------------------//
// CompatibleFaceSpecificThermalEnergyPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#include "CompatibleFaceSpecificThermalEnergyPolicy.hh"
#include "computeSVPHCorrectionsOnFaces.hh"
#include "Hydro/HydroFieldNames.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Kernel/TableKernel.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"

#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Define the weighting function deciding how to divvy up the work between
// nodes.
//------------------------------------------------------------------------------
// The old method.
// inline
// double standardWeighting(const double& ui,
//                          const double& uj,
//                          const double& mi,
//                          const double& mj,
//                          const double& duij) {
//   const double uji = uj - ui;
//   const double result = 0.5*(1.0 + uji/(abs(uji) + 1.0/(1.0 + abs(uji)))*sgn0(duij));
//   ENSURE(result >= 0.0 and result <= 1.0);
//   return result;
// }

inline
double standardWeighting(const double& ui,
                         const double& uj,
                         const double& /*mi*/,
                         const double& /*mj*/,
                         const double& duij) {
  const double uji = uj - ui;
  const double result = 0.5*(1.0 + uji*sgn0(duij)/(abs(uji) + 0.5));
  ENSURE(result >= 0.0 and result <= 1.0);
  return result;
}

// A strictly monotonic scheme that *will* pull the specific thermal energies
// together.
inline
double monotonicWeighting(const double& ui,
                          const double& uj,
                          const double& mi,
                          const double& mj,
                          const double& DEij) {
  REQUIRE(mi > 0.0);
  REQUIRE(mj > 0.0);
  const double uji = uj - ui;
  const double miInv = 1.0/mi;
  const double mjInv = 1.0/mj;

  double result;
  if (uji*DEij >= 0.0) {
    const double A = min(1.0, max(0.0, uji*mi*safeInv(DEij, 1.0e-50)));
    const double B = min(1.0, max(0.0, mjInv/(miInv + mjInv)));
    CHECK(A >= 0.0 and A <= 1.0);
    CHECK(B >= 0.0 and B <= 1.0);
    result = A + B*(1.0 - A);
  } else {
    const double A = min(1.0, max(0.0, -uji*mj*safeInv(DEij, 1.0e-50)));
    const double B = min(1.0, max(0.0, miInv/(miInv + mjInv)));
    CHECK(A >= 0.0 and A <= 1.0);
    CHECK(B >= 0.0 and B <= 1.0);
    result = 1.0 - A - B*(1.0 - A);
  }
  ENSURE(result >= 0.0 and result <= 1.0);
  return result;
}

inline
double weighting(const double& ui,
                 const double& uj,
                 const double& mi,
                 const double& mj,
                 const double& duij,
                 const double& dt) {

  // First use our standard weighting algorithm.
  const double fi = standardWeighting(ui, uj, mi, mj, duij);
  CHECK(fi >= 0.0 and fi <= 1.0);

  // Now the monotonic weighting.
  const double mfi = monotonicWeighting(ui, uj, mi, mj, mi*duij*dt);

  // Combine them for the final answer.
  const double chi = abs(ui - uj)/(abs(ui) + abs(uj) + 1.0e-50);
  CHECK(chi >= 0.0 and chi <= 1.0);
  return chi*mfi + (1.0 - chi)*fi;
}

// inline
// double weighting(const double& ui,
//                  const double& uj,
//                  const double& mi,
//                  const double& mj,
//                  const double& duij,
//                  const double& dt) {

//   // First use our standard weighting algorithm.
//   const double fi = standardWeighting(ui, uj, mi, mj, duij);
//   const double fj = 1.0 - fi;
//   CHECK(fi >= 0.0 and fi <= 1.0);
//   CHECK(fj >= 0.0 and fj <= 1.0);
//   CHECK(fuzzyEqual(fi + fj, 1.0));

//   // Check if this would result in one (but not both) energies
//   // flipping sign.  If so, we revert to the strictly monotonic, but not
//   // quite as accurate scheme.
//   const double ui1 = ui + fi*duij*dt;
//   const double uj1 = uj + fj*duij*dt;
//   if (sgn(ui*uj) == sgn(ui1*uj1)) {
//     return fi;
//   } else {
//     return monotonicWeighting(ui, uj, mi, mj, mi*duij*dt);
//   }
// }

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
CompatibleFaceSpecificThermalEnergyPolicy(const TableKernel<Dimension>& W,
                                          const DataBase<Dimension>& dataBase,
                                          ConstBoundaryIterator boundaryBegin,
                                          ConstBoundaryIterator boundaryEnd):
  IncrementState<Dimension, typename Dimension::Scalar>(),
  mW(W),
  mDataBase(dataBase),
  mBoundaryBegin(boundaryBegin),
  mBoundaryEnd(boundaryEnd) {
  mFired = false;
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double dt) {
  using Zone = typename Mesh<Dimension>::Zone;
  using Face = typename Mesh<Dimension>::Face;

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  if (not mFired) {
    mFired = true;

    // Get the state fields.
    FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
    const FieldList<Dimension, Scalar> specificEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
    const FieldList<Dimension, vector<Vector> > faceForce = derivs.fields(HydroFieldNames::faceForce, vector<Vector>());

    // The acceleration is special -- we need the ghost state for this derivative value so we'll force
    // it to be updated now.
    FieldList<Dimension, Vector> acceleration = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero());
    for (ConstBoundaryIterator itr = mBoundaryBegin;
         itr != mBoundaryEnd;
         ++itr) {
      (*itr)->applyFieldListGhostBoundary(acceleration);
    }
    for (ConstBoundaryIterator itr = mBoundaryBegin;
         itr != mBoundaryEnd;
         ++itr) {
      (*itr)->finalizeGhostBoundary();
    }

    const double hdt = 0.5*multiplier;
    const size_t numNodeLists = eps.numFields();

    const Mesh<Dimension>& mesh = state.mesh();

    // Walk the nodes of this NodeList
    unsigned i, k, n, nodeListi;
    Scalar DepsDti, duij;
    Vector vi12, vj12, vji12;
    for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      n = eps[nodeListi]->numInternalElements();
      for (i = 0; i != n; ++i) {
        const Zone& zonei = mesh.zone(nodeListi, i);
        const vector<int>& faceIDs = zonei.faceIDs();
        const unsigned nfaces = faceIDs.size();

        // State for node i.
        DepsDti = 0.0;
        const Scalar& mi = mass(nodeListi, i);
        const Vector& vi = velocity(nodeListi, i);
        const Scalar& ui = specificEnergy0(nodeListi, i);
        const Vector& ai = acceleration(nodeListi, i);
        const vector<Vector>& fforce = faceForce(nodeListi, i);
        CHECK2(fforce.size() == nfaces, fforce.size() << " " << nfaces);
        vi12 = vi + ai*hdt;

        // Walk the faces of this SVPH node.
        for (k = 0; k != nfaces; ++k) {
          const Face& face = mesh.face(faceIDs[k]);
          unsigned nodeListj = nodeListi, j = i;
          const unsigned oppZoneID = Mesh<Dimension>::positiveID(face.oppositeZoneID(zonei.ID()));
          if (oppZoneID != Mesh<Dimension>::UNSETID) {
            mesh.lookupNodeListID(oppZoneID, nodeListj, j);
          }

          // State for opposite node.
          const Scalar& mj = mass(nodeListj, j);
          const Vector& vj = velocity(nodeListj, j);
          const Scalar& uj = specificEnergy0(nodeListj, j);
          const Vector& aj = acceleration(nodeListj, j);
          vj12 = vj + aj*hdt;
          vji12 = vj12 - vi12;

          duij = vji12.dot(fforce[k])/mi;
          const Scalar wi = weighting(ui, uj, mi, mj, duij, dt);
          // if (abs(duij) > 1.0e-10) {
          //   cerr << "    " << i << " " << Mesh<Dimension>::positiveID(faceIDs[k]) << " : "
          //        << duij << " " << wi << " " << wi*duij << " " << (1.0 - wi)*duij << " : "
          //        << endl;
          // }
          CHECK(wi >= 0.0 and wi <= 1.0);
          CHECK(fuzzyEqual(wi + weighting(uj, ui, mj, mi, duij*mi/mj, dt), 1.0, 1.0e-10));
          DepsDti += wi*duij;
        }

        // Update the work on this node.
        eps(nodeListi, i) += multiplier*DepsDti;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a CompatibleSpecificThermalEnergy operator.
  const CompatibleFaceSpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const CompatibleFaceSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename Dimension> bool CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::mFired = false;

}

