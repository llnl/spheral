//------------------------------------------------------------------------------
// Use SVPH to sample a FieldList.
//------------------------------------------------------------------------------
#include "sampleFieldListSVPH.hh"
#include "computeSVPHCorrections.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/safeInv.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/MathTraits.hh"

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

namespace {
//------------------------------------------------------------------------------
// Convert 1D types to scalars.
//------------------------------------------------------------------------------
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::Scalar& val) { return val; }
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::Vector& val) { return val.x(); }
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::Tensor& val) { return val.xx(); }
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::SymTensor& val) { return val.xx(); }
  
}

//------------------------------------------------------------------------------
// 1-D gradient in cell.
//------------------------------------------------------------------------------
template<typename DataType>
typename MathTraits<Dim<1>, DataType>::GradientType
computeCellGradient(const FieldList<Dim<1>, DataType>& fieldList,
                    const FieldList<Dim<1>, Dim<1>::Vector>& position,
                    const int nodeListi,
                    const int i,
                    const Mesh<Dim<1> >& mesh) {
  typedef Dim<1> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Mesh<Dimension>::Zone Zone;
  typedef Mesh<Dimension>::Face Face;
  typedef typename MathTraits<Dim<1>, DataType>::GradientType GradientType;

  // Find the neighboring mesh cells.
  const Zone& zonei = mesh.zone(nodeListi, i);
  const std::vector<int>& faces = zonei.faceIDs();
  CHECK(faces.size() == 2);
  const int face1ID = Mesh<Dimension>::positiveID(faces[0]),
            face2ID = Mesh<Dimension>::positiveID(faces[1]);
  const Face face1 = mesh.face(face1ID),
             face2 = mesh.face(face2ID);
  const unsigned zone1ID = Mesh<Dimension>::positiveID(face1.oppositeZoneID(zonei.ID())),
                 zone2ID = Mesh<Dimension>::positiveID(face2.oppositeZoneID(zonei.ID()));

  // Look up the indices into the FieldList for the two neighboring cells.
  unsigned nodeList1ID = nodeListi, node1ID = i, 
           nodeList2ID = nodeListi, node2ID = i;
  if (zone1ID != Mesh<Dimension>::UNSETID) mesh.lookupNodeListID(zone1ID, nodeList1ID, node1ID);
  if (zone2ID != Mesh<Dimension>::UNSETID) mesh.lookupNodeListID(zone2ID, nodeList2ID, node2ID);

  // Now compute the weighted intrazonal gradient.
  const Scalar fi = scalarAdapter(fieldList(nodeListi, i)),
               f1 = scalarAdapter(fieldList(nodeList1ID, node1ID)),
               f2 = scalarAdapter(fieldList(nodeList2ID, node2ID));
  const Scalar dx01 = (position(nodeList1ID, node1ID) - position(nodeListi, i)).x(),
               dx02 = (position(nodeList2ID, node2ID) - position(nodeListi, i)).x();
  CHECK(std::abs(dx01) + std::abs(dx02) > 0.0);
  const GradientType result((sgn(dx01)*(f1 - fi) + sgn(dx02)*(f2 - fi))*
                            safeInv(std::abs(dx01) + std::abs(dx02)));
  return result;
}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
sampleFieldListSVPH(const FieldList<Dimension, DataType>& fieldList,
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

  // Prepare the result and some work fields.
  FieldList<Dimension, DataType> result(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> volume(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> A(FieldStorageType::CopyFields);
  FieldList<Dimension, Vector> B(FieldStorageType::CopyFields);
  FieldList<Dimension, Tensor> gradB(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();
    result.appendNewField("SVPH sample of " + fieldList[nodeListi]->name(), nodeList, DataType(0));
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
      const Scalar Hdeti = Hi.Determinant();
      Scalar norm = Vi*W0*Hdeti;
      result(nodeListi, i) = Vi*W0*Hdeti * fieldList(nodeListi, i);

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
          const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);

          // Increment the result.
          const Scalar VWRj = Vj*(1.0 + Bi.dot(rij))*Wj;
          norm += VWRj;
          result(nodeListi, i) += VWRj * Fj;
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
