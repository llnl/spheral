//------------------------------------------------------------------------------
// Modified ideas from CRKSPH using the SVPH normalized kernel estimates.
// This version computes the corrections on the faces of a mesh, rather than on
// the points themselves.
//------------------------------------------------------------------------------
#include "computeSVPHCorrectionsOnFaces.hh"
#include "Field/Field.hh"
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
                              const TableKernel<Dimension>& W,
                              const FieldList<Dimension, typename Dimension::Scalar>& volume,
                              const FieldList<Dimension, typename Dimension::Vector>& position,
                              const FieldList<Dimension, typename Dimension::SymTensor>& H,
                              const BoundaryIterator& /*boundaryBegin*/,
                              const BoundaryIterator& /*boundaryEnd*/,
                              vector<typename Dimension::Scalar>& A,
                              vector<typename Dimension::Vector>& B) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Mesh<Dimension>::Face Face;

  // Pre-conditions.
  const size_t numNodeLists = volume.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Extract the NodeLists we're working on.
  const vector<NodeList<Dimension>*> nodeLists = volume.nodeListPtrs();

  // Zero out the result.
  const unsigned nfaces = mesh.numFaces();
  A = vector<Scalar>(nfaces, 0.0);
  B = vector<Vector>(nfaces, Vector::zero());

  // We can derive everything in terms of the first and second moments 
  // of the local positions.
  vector<Vector> m1(nfaces, Vector::zero());
  vector<SymTensor> m2(nfaces, SymTensor::zero());

  // Walk the faces.
  unsigned iface, nodeListj, j;
  Vector posFace;
  const SymTensor Hface = 1.0e100*SymTensor::one();
  for (iface = 0; iface != nfaces; ++iface) {
    const Face& face = mesh.face(iface);
    posFace = face.position();

    // Set the neighbors for this face.
    vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors;
    Neighbor<Dimension>::setMasterNeighborGroup(posFace, Hface,
                                                nodeLists.begin(), nodeLists.end(),
                                                W.kernelExtent(),
                                                masterLists,
                                                coarseNeighbors);

    // Iterate over the NodeLists.
    for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
      Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
      neighbor.setRefineNeighborList(posFace, Hface, coarseNeighbors[nodeListj], refineNeighbors[nodeListj]);
      for (auto neighborItr = refineNeighbors[nodeListj].begin(); neighborItr != refineNeighbors[nodeListj].end(); ++neighborItr) {
        j = *neighborItr;

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

