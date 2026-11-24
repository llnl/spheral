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
#include "Utilities/range.hh"
#include "Utilities/safeInv.hh"
#include "Mesh/Mesh.hh"

#include <algorithm>

using std::vector;

namespace Spheral {

template<typename Dimension>
void
computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dimension>& mesh,
                                          const TableKernel<Dimension>& W,
                                          const DataBase<Dimension>& dataBase,
                                          FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using Zone = typename Mesh<Dimension>::Zone;
  using Face = typename Mesh<Dimension>::Face;
  using NodeID = std::pair<size_t, size_t>;

  // State fields.
  const FieldList<Dimension, typename Dimension::Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, typename Dimension::Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, typename Dimension::SymTensor> H = dataBase.fluidHfield();

  const auto  numFaces = mesh.numFaces();
  const auto  H0 = 1.0e100*SymTensor::one();
  const auto& neighbor = dataBase.neighbor();
  const auto  numNodeLists = dataBase.numFluidNodeLists();
  REQUIRE(dataBase.numNodeLists() == numNodeLists);

  // Prepare the face fields for summing.
  vector<vector<Scalar>> massFace(numNodeLists, vector<Scalar>(numFaces, 0.0)), volFace(numNodeLists, vector<Scalar>(numFaces, 0.0));
  vector<Vector> posFace(numFaces, Vector::zero());

  // Walk the faces and sum the mass and volume properties to faces
  for (auto iface = 0u; iface < numFaces; ++iface) {
    const Face& face = mesh.face(iface);
    posFace[iface] = face.position();

    // Set the neighbors for this face.
    vector<NodeID> masterLists, coarseNeighbors, refineNeighbors;
    auto [ok, numMaster, numCoarse] = visitNeighbor(neighbor, [&](const auto& x) { return x.setMasterList(posFace[iface], H0,
                                                                                                          masterLists,
                                                                                                          coarseNeighbors); });
    CHECK(ok);
    CONTRACT_VAR(numMaster);
    refineNeighbors.resize(numCoarse);
    const auto numRefine = visitNeighbor(neighbor, [&](const auto& x) { return x.setRefineNeighborList(posFace[iface], H0, numCoarse, coarseNeighbors, refineNeighbors); });

    // Walk the neighbors
    for (auto k = 0u; k < numRefine; ++k) {
      const auto [nodeListj, j] = refineNeighbors[k];
      
      // Get the state for node j
      const Vector& rj = position(nodeListj, j);
      const Scalar& mj = mass(nodeListj, j);
      const SymTensor& Hj = H(nodeListj, j);
      const Scalar Vj = mesh.zone(nodeListj, j).volume();
      const Scalar Hdetj = Hj.Determinant();
      CHECK(Vj > 0.0);
      CHECK(Hdetj > 0.0);

      // Pair-wise kernel type stuff.
      const Vector rij = posFace[iface] - rj;
      const Vector etaj = Hj*rij;
      const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);

      // Contribution to this face.
      massFace[nodeListj][iface] += mj*Wj;
      volFace[nodeListj][iface] += Vj*Wj;
    }
  }

  // Interpolate from faces to SVPH nodes.
  for (auto [nodeListi, nptr]: enumerate(dataBase.fluidNodeListPtrs())) {
    const auto n = nptr->numInternalNodes();
    const auto rhoMin = nptr->rhoMin();
    const auto rhoMax = nptr->rhoMax();
    for (auto i = 0u; i < n; ++i) {
      const Zone& zone = mesh.zone(nodeListi, i);
      const vector<int>& faceIDs = zone.faceIDs();
      const Vector rc = position(nodeListi, i); // zone.position();
      Scalar Meff = 0.0, Veff = 0.0;
      for (auto k = 0u; k < faceIDs.size(); ++k) {
        const unsigned iface = Mesh<Dimension>::positiveID(faceIDs[k]);
        const Face& face = mesh.face(iface);
        const Vector rf = face.position();
        const Scalar Af = face.area();
        const Scalar wf = Af*(rf - rc).magnitude();
        Meff += massFace[nodeListi][iface] * wf;
        Veff += volFace[nodeListi][iface] * wf;
      }
      CHECK(Meff > 0.0);
      CHECK(Veff > 0.0);
      massDensity(nodeListi, i) = std::clamp(Meff/Veff, rhoMin, rhoMax);
    }
  }
}

}
