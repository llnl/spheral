//------------------------------------------------------------------------------
// Compute the Voronoi cell mass density summation.
//------------------------------------------------------------------------------
#include "computeSumVoronoiCellMassDensityFromFaces.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/safeInv.hh"
#include "Mesh/Mesh.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dimension>& mesh,
                                          const TableKernel<Dimension>& W,
                                          const DataBase<Dimension>& dataBase,
                                          FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Face Face;

  // State fields.
  const FieldList<Dimension, typename Dimension::Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, typename Dimension::Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, typename Dimension::SymTensor> H = dataBase.fluidHfield();

  const unsigned numFaces = mesh.numFaces();
  const SymTensor H0 = 1.0e100*SymTensor::one();
  unsigned nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = **itr;
    Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
    const Scalar rhoMin = nodeList.rhoMin();
    const Scalar rhoMax = nodeList.rhoMax();
    const unsigned n = nodeList.numInternalNodes();

    // Prepare the face fields for summing.
    vector<Scalar> massFace(numFaces, 0.0), volFace(numFaces, 0.0);
    vector<Vector> posFace(numFaces, Vector::zero());

    // Walk the faces and sum the mass and volume properties to faces for this material.
    for (unsigned iface = 0; iface != numFaces; ++iface) {
      const Face& face = mesh.face(iface);
      posFace[iface] = face.position();

      // Set the neighbors for this face.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors;
      Neighbor<Dimension>::setMasterNeighborGroup(posFace[iface], H0,
                                                  dataBase.fluidNodeListBegin(),
                                                  dataBase.fluidNodeListEnd(),
                                                  W.kernelExtent(),
                                                  masterLists,
                                                  coarseNeighbors);
      neighbor.setRefineNeighborList(posFace[iface], H0, coarseNeighbors[nodeListi], refineNeighbors[nodeListi]);

      // Walk the neighbors in this NodeList.
      for (auto neighborItr = refineNeighbors[nodeListi].begin(); neighborItr != refineNeighbors[nodeListi].end(); ++neighborItr) {
        unsigned j = *neighborItr;
      
        // Get the state for node j
        const Vector& rj = position(nodeListi, j);
        const Scalar& mj = mass(nodeListi, j);
        const SymTensor& Hj = H(nodeListi, j);
        const Scalar Vj = mesh.zone(nodeListi, j).volume();
        const Scalar Hdetj = Hj.Determinant();
        CHECK(Vj > 0.0);
        CHECK(Hdetj > 0.0);

        // Pair-wise kernel type stuff.
        const Vector rij = posFace[iface] - rj;
        const Vector etaj = Hj*rij;
        const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);

        // Contribution to this face.
        massFace[iface] += mj*Wj;
        volFace[iface] += Vj*Wj;
      }
    }

    // Interpolate from faces to SVPH nodes.
    for (unsigned i = 0; i != n; ++i) {
      const Zone& zone = mesh.zone(nodeListi, i);
      const vector<int>& faceIDs = zone.faceIDs();
      const Vector rc = position(nodeListi, i); // zone.position();
      Scalar Meff = 0.0, Veff = 0.0;
      for (unsigned k = 0; k != faceIDs.size(); ++k) {
        const unsigned iface = Mesh<Dimension>::positiveID(faceIDs[k]);
        const Face& face = mesh.face(iface);
        const Vector rf = face.position();
        const Scalar Af = face.area();
        const Scalar wf = Af*(rf - rc).magnitude();
        Meff += massFace[iface] * wf;
        Veff += volFace[iface] * wf;
      }
      CHECK(Meff > 0.0);
      CHECK(Veff > 0.0);
      massDensity(nodeListi, i) = max(rhoMin, min(rhoMax, Meff/Veff));
    }
  }
}

}

