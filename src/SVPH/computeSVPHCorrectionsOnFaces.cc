//------------------------------------------------------------------------------
// Modified ideas from CRKSPH using the SVPH normalized kernel estimates.
// This version computes the corrections on the faces of a mesh, rather than on
// the points themselves.
//------------------------------------------------------------------------------
#include "computeSVPHCorrectionsOnFaces.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"
#include "Mesh/Mesh.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Compute corrections.
//------------------------------------------------------------------------------
template<typename Dimension, typename BoundaryIterator>
void
computeSVPHCorrectionsOnFaces(const Mesh<Dimension>& mesh,
                              const DataBase<Dimension>& dataBase,
                              const TableKernel<Dimension>& W,
                              const FieldList<Dimension, typename Dimension::Scalar>& volume,
                              const FieldList<Dimension, typename Dimension::Vector>& position,
                              const FieldList<Dimension, typename Dimension::SymTensor>& H,
                              const BoundaryIterator& /*boundaryBegin*/,
                              const BoundaryIterator& /*boundaryEnd*/,
                              vector<typename Dimension::Scalar>& A,
                              vector<typename Dimension::Vector>& B) {

  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using Face = typename Mesh<Dimension>::Face;
  using NodeID = std::pair<size_t, size_t>;

  // Pre-conditions.
  const size_t numNodeLists = dataBase.numFluidNodeLists();
  REQUIRE(dataBase.numNodeLists() == numNodeLists);
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  const auto& neighbor = dataBase.neighbor();

  // Zero out the result.
  const unsigned nfaces = mesh.numFaces();
  A = vector<Scalar>(nfaces, 0.0);
  B = vector<Vector>(nfaces, Vector::zero());

  // We can derive everything in terms of the first and second moments 
  // of the local positions.
  vector<Vector> m1(nfaces, Vector::zero());
  vector<SymTensor> m2(nfaces, SymTensor::zero());

  // Walk the faces.
  unsigned iface;
  Vector posFace;
  const SymTensor Hface = 1.0e100*SymTensor::one();
  for (iface = 0; iface != nfaces; ++iface) {
    const Face& face = mesh.face(iface);
    posFace = face.position();

    // Set the neighbors for this face.
    vector<NodeID> masterLists, coarseNeighbors, refineNeighbors;
    auto [ok, numMaster, numCoarse] = visitNeighbor(neighbor, [&](const auto& x) { return x.setMasterList(posFace, Hface,
                                                                                                          masterLists,
                                                                                                          coarseNeighbors); });
    CHECK(ok);
    CONTRACT_VAR(numMaster);
    refineNeighbors.resize(numCoarse);
    const auto numRefine = visitNeighbor(neighbor, [&](const auto& x) { return x.setRefineNeighborList(posFace, Hface, numCoarse, coarseNeighbors, refineNeighbors); });

    // Iterate over the neighbors
    for (auto k = 0u; k < numRefine; ++k) {
      const auto [nodeListj, j] = refineNeighbors[k];

      // Get the state for node j
      const Vector& rj = position(nodeListj, j);
      const Scalar& Vj = volume(nodeListj, j);
      const SymTensor& Hj = H(nodeListj, j);
      const Scalar Hdetj = Hj.Determinant();
      CHECK(Vj > 0.0);
      CHECK(Hdetj > 0.0);

      // Pair-wise kernel type stuff.
      const Vector rij = posFace - rj;
      const Vector etaj = Hj*rij;
      const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);

      // Normalization (zeroth moment).
      A[iface] += Vj*Wj;

      // First moment.
      m1[iface] += Vj*Wj * rij;

      // Second moment.
      const SymTensor thpt = rij.selfdyad();
      m2[iface] += Vj*Wj * thpt;
    }
  }

  // // Apply any boundary action to the sums.
  // for (BoundaryIterator itr = boundaryBegin;
  //      itr != boundaryEnd;
  //      ++itr) {
  //   (*itr)->enforceBoundary(A, mesh);
  //   (*itr)->enforceBoundary(m1, mesh);
  //   (*itr)->enforceBoundary(m2, mesh);
  // }

  // Based on the moments we can calculate the SVPH corrections terms and their gradients.
  for (iface = 0; iface != nfaces; ++iface) {
    CHECK(A[iface] >= 0.0);
    if (A[iface] > 0.0) {
      A[iface] = 1.0/A[iface];

      // CHECK2(abs(m2[iface].Determinant()) > 1.0e-30, iface << " " << m2[iface].Determinant());
      if (abs(m2[iface].Determinant()) > 1.0e-30) {
        const SymTensor m2inv = m2[iface].Inverse();
        B[iface] = -(m2inv*m1[iface]);
      }
    }
  }
}

}
