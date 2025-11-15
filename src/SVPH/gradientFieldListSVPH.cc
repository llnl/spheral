//------------------------------------------------------------------------------
// Use SVPH to take the gradient of a FieldList.
//------------------------------------------------------------------------------
#include "sampleFieldListSVPH.hh"
#include "computeSVPHCorrections.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/safeInv.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/innerProduct.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Local utility methods.
//------------------------------------------------------------------------------
namespace {

// 1D
static inline
Dim<1>::Vector
delta_grad(const Dim<1>::Vector& xij,
           const Dim<1>::Vector& Bi,
           const Dim<1>::Tensor& gradBi,
           const Dim<1>::Scalar& Wj,
           const Dim<1>::Vector& gradWj) {
  return Dim<1>::Vector((1.0 + Bi(0)*xij(0))*gradWj(0) + (Bi(0) + xij(0)*gradBi(0,0))*Wj);
}

// 2D
static inline
Dim<2>::Vector
delta_grad(const Dim<2>::Vector& xij,
           const Dim<2>::Vector& Bi,
           const Dim<2>::Tensor& gradBi,
           const Dim<2>::Scalar& Wj,
           const Dim<2>::Vector& gradWj) {
  return Dim<2>::Vector((1.0 + Bi(0)*xij(0) + Bi(1)*xij(1))*gradWj(0) + (Bi(0) + xij(0)*gradBi(0,0) + xij(1)*gradBi(1,0))*Wj,
                        (1.0 + Bi(0)*xij(0) + Bi(1)*xij(1))*gradWj(1) + (Bi(1) + xij(0)*gradBi(0,1) + xij(1)*gradBi(1,1))*Wj);
}

// 3D
static inline
Dim<3>::Vector
delta_grad(const Dim<3>::Vector& /*xij*/,
           const Dim<3>::Vector& /*Bi*/,
           const Dim<3>::Tensor& /*gradBi*/,
           const Dim<3>::Scalar& /*Wj*/,
           const Dim<3>::Vector& /*gradWj*/) {
  return Dim<3>::Vector();
}

}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientFieldListSVPH(const FieldList<Dimension, DataType>& fieldList,
                      const FieldList<Dimension, typename Dimension::Vector>& position,
                      const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                      const ConnectivityMap<Dimension>& connectivityMap,
                      const TableKernel<Dimension>& W,
                      const Mesh<Dimension>& mesh,
                      const bool firstOrderConsistent) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(Hfield.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Prepare the result and some work fields.
  FieldList<Dimension, GradientType> result(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> volume(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> A(FieldStorageType::CopyFields);
  FieldList<Dimension, Vector> B(FieldStorageType::CopyFields);
  FieldList<Dimension, Tensor> gradB(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();
    result.appendNewField("SVPH gradient of " + fieldList[nodeListi]->name(), nodeList, GradientType());
    volume.appendNewField("volume", nodeList, 0.0);
    A.appendNewField("SVPH normalization for " + fieldList[nodeListi]->name(), nodeList, 0.0);
    B.appendNewField("SVPH linear correction for " + fieldList[nodeListi]->name(), nodeList, Vector::zero());
    gradB.appendNewField("SVPH linear correction gradient for " + fieldList[nodeListi]->name(), nodeList, Tensor::zero());
  }

  // If we're enforcing first-order consistency, compute the correction fields.
  if (firstOrderConsistent) {

    // Copy the cell volumes to the FieldList.
    for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
      const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();
      for (auto i = 0u; i != nodeList.numNodes(); ++i) {
        volume(nodeListi, i) = mesh.zone(nodeListi, i).volume();
      }
    }

    // We have a handy utility method to compute the corrections.
    computeSVPHCorrections(connectivityMap,
                           W,
                           volume,
                           position,
                           Hfield,
                           A,
                           B,
                           gradB);
  }

  // Walk the NodeLists to build our answer.
  const Scalar W0 = W.kernelValue(0.0, 1.0);
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = Hfield(nodeListi, i);
      const Scalar Vi = mesh.zone(nodeListi, i).volume();
      const Vector& Bi = B(nodeListi, i);
      const Tensor& gradBi = gradB(nodeListi, i);
      const DataType& Fi = fieldList(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      Scalar norm = Vi*W0*Hdeti;

      // Walk the neighbors for this node.
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      for (auto nodeListj = 0u; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          // Get the state for node j.
          const DataType& Fj = fieldList(nodeListj, j);
          const Vector& rj = position(nodeListj, j);
          const SymTensor& Hj = Hfield(nodeListj, j);
          const Scalar Vj = mesh.zone(nodeListj, j).volume();
          const Scalar Hdetj = Hj.Determinant();

          // Pair-wise kernel type stuff.
          const Vector rij = ri - rj;
          const Vector etaj = Hj*rij;
          const Vector Hetaj = Hj*etaj.unitVector();
          Scalar Wj, gWj;
          W.kernelAndGradValue(etaj.magnitude(), Hdetj, Wj, gWj);
          const Vector gradWj = gWj*Hetaj;

          // Increment the result.
          norm += Vj*(1.0 + Bi.dot(rij))*Wj;
          result(nodeListi, i) += Vj*(Fj - Fi)*delta_grad(rij, Bi, gradBi, Wj, gradWj);
          // result(nodeListi, i) += Vj*(Fj - Fi)*((1.0 + Bi.dot(rij))*gradWj + graddot(rij, Bi, gradBi)*Wj);
          // result(nodeListi, i) += Vj*(Fj - Fi)*((1.0 + Bi.dot(rij))*gradWj + (Bi + innerProduct<Dimension>(rij, gradBi))*Wj);
        }
      }

      // Normalize our ith value.
      CHECK(norm > 0.0);
      result(nodeListi, i) /= norm;
    }
  }

  // That's it.
  return result;
}

}
