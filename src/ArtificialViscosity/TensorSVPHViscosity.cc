//---------------------------------Spheral++----------------------------------//
// A version of our tensor viscosity specialized for the SVPHFacetedHydro
// algorithm.
//
// Created by J. Michael Owen, Sat Aug 31 13:31:51 PDT 2013
//----------------------------------------------------------------------------//
#include "TensorSVPHViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/EigenStruct.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"

#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorSVPHViscosity<Dimension>::
TensorSVPHViscosity(const Scalar Clinear,
                    const Scalar Cquadratic,
                    const TableKernel<Dimension>& WT,
                    const Scalar fslice):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic, WT),
  mfslice(fslice),
  mDvDx(),
  mShearCorrection(),
  mQface() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TensorSVPHViscosity<Dimension>::
initialize(const Scalar t,
           const Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  using Face = typename Mesh<Dimension>::Face;
  const auto& W = this->kernel();

  // The set of NodeLists.
  const vector<const NodeList<Dimension>*> nodeLists(dataBase.fluidNodeListBegin(),
                                                     dataBase.fluidNodeListEnd());
  const unsigned numNodeLists = nodeLists.size();

  // Grab the state we need.
  const Mesh<Dimension>& mesh = state.mesh();
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero());
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());

  // Prepare the state we're going to compute.
  const unsigned numFaces = mesh.numFaces();
  mDvDx = vector<Tensor>(numFaces, Tensor::zero());
  mShearCorrection = vector<Scalar>(numFaces, 1.0);
  mQface = vector<Tensor>(numFaces, Tensor::zero());

  // Walk the faces.
  const Scalar Cl = this->Cl();
  const Scalar Cq = this->Cq();
  const bool balsaraCorrection = this->balsaraShearCorrection();
  const Scalar csMin = this->negligibleSoundSpeed();
  const Scalar eps2 = this->epsilon2();
  const SymTensor H0 = 1.0e100*SymTensor::one();
  for (unsigned iface = 0; iface != numFaces; ++iface) {
    const Face& face = mesh.face(iface);
    const Vector posFace = face.position();
    Scalar volSum = 1.0e-30;

    // Look up the SVPH nodes either side of the Face.
    unsigned z1id = Mesh<Dimension>::positiveID(face.zone1ID());
    unsigned z2id = Mesh<Dimension>::positiveID(face.zone2ID());
    if (z1id > z2id) std::swap(z1id, z2id);
    CHECK(z1id != Mesh<Dimension>::UNSETID);
    unsigned nodeListi, nodeListj, i, j;
    mesh.lookupNodeListID(z1id, nodeListi, i);
    if (z2id == Mesh<Dimension>::UNSETID) {
      nodeListj = nodeListi;
      j = i;
    } else {
      mesh.lookupNodeListID(z2id, nodeListj, j);
    }
    const Vector vface = 0.5*(velocity(nodeListi, i) + velocity(nodeListj, j));
    const Scalar rhoFace = 0.5*(massDensity(nodeListi, i) + massDensity(nodeListj, j));
    const Scalar csFace = 0.5*(soundSpeed(nodeListi, i) + soundSpeed(nodeListj, j));
    const Scalar hface = 0.5*(H(nodeListi, i).Inverse() + H(nodeListj, j).Inverse()).Trace() / Dimension::nDim;
    const Scalar hface2 = hface*hface;

    // Set the neighbors for this face.
    vector<vector<int>> masterLists, coarseNeighbors;
    Neighbor<Dimension>::setMasterNeighborGroup(posFace, H0,
                                                nodeLists.begin(), nodeLists.end(),
                                                W.kernelExtent(),
                                                masterLists,
                                                coarseNeighbors);

    // Iterate over the NodeLists.
    for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
      Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
      vector<int> refineNeighbors;
      neighbor.setRefineNeighborList(posFace, H0, coarseNeighbors[nodeListj], refineNeighbors);
      for (auto neighborItr = refineNeighbors.begin(); neighborItr < refineNeighbors.end(); ++neighborItr) {
        const unsigned j = *neighborItr;
      
        // Get the state for node j
        const Vector& rj = position(nodeListj, j);
        const Vector& vj = velocity(nodeListj, j);
        const SymTensor& Hj = H(nodeListj, j);
        const Scalar& Vj = volume(nodeListj, j);
        const Scalar Hdetj = Hj.Determinant();
        CHECK(Vj > 0.0);
        CHECK(Hdetj > 0.0);

        // Pair-wise kernel type stuff.
        const Vector rij = posFace - rj;
        const Vector etaj = Hj*rij;
        const Vector Hetaj = Hj*etaj.unitVector();
        Scalar Wj, gWj;
        W.kernelAndGradValue(etaj.magnitude(), Hdetj, Wj, gWj);
        const Vector gradWj = gWj*Hetaj;

        // Increment the face fluid properties.
        volSum += Vj*Wj;
        mDvDx[iface] += (vface - vj)*Vj*gradWj;
      }
    }

    // Finish the background velocity gradient.
    CHECK(volSum > 0.0);
    mDvDx[iface] /= -volSum;

    // Add in any of the slice term from these two nodes.
    if (mfslice > 0.0 and (nodeListi != nodeListj or i != j)) {
      const Vector rij = position(nodeListi, i) - position(nodeListj, j);
      const Vector vij = velocity(nodeListi, i) - velocity(nodeListj, j);
      const Vector rijUnit = rij.unitVector();
      const Scalar rij2 = rij.magnitude2();
      const Tensor R = rotationMatrix(rijUnit);
      const Tensor Rinverse = R.Transpose();
      const Vector thpt1 = sqrt(rij2)*(R*vij);
      const Vector deltaSigma = thpt1/(rij2 + eps2*hface2);
      mDvDx[iface].rotationalTransform(R);
      mDvDx[iface].setColumn(0, mfslice*deltaSigma + (1.0 - mfslice)*mDvDx[iface].getColumn(0));
      mDvDx[iface].rotationalTransform(Rinverse);
    }

    // Now compute the shear correction.
    if (balsaraCorrection) {
      const Scalar div = abs(mDvDx[iface].Trace());
      const Scalar curl = this->curlVelocityMagnitude(mDvDx[iface]);
      const Scalar cs = max(csMin, csFace);
      CHECK(div >= 0.0);
      CHECK(curl >= 0.0);
      CHECK(cs > 0.0);
      mShearCorrection[iface] = div/(div + curl + eps2*cs/hface);
      CHECK(mShearCorrection[iface] >= 0.0 and mShearCorrection[iface] <= 1.0);
    }

    // Compute the symmetric, compressing portion of the velocity gradient.
    SymTensor muface = mDvDx[iface].Symmetric();
    const typename SymTensor::EigenStructType es = muface.eigenVectors();
    muface.Zero();
    for (unsigned j = 0; j != Dimension::nDim; ++j) muface(j,j) = min(0.0, es.eigenValues(j));
    muface.rotationalTransform(es.eigenVectors);
    muface *= hface;
    
    // Now we can compute the Q on this face.
    CHECK(fuzzyLessThanOrEqual(muface.Trace(), 0.0, 1.0e-8));
    mQface[iface] = mShearCorrection[iface]*rhoFace*(-Cl*csFace*muface + Cq*muface*muface);
  }
  return false;
}

}

