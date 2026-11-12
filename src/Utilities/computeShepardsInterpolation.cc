//------------------------------------------------------------------------------
// Compute the Shepards interpolation of a FieldList.
//------------------------------------------------------------------------------
#include "computeShepardsInterpolation.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/DataTypeTraits.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
computeShepardsInterpolation(const FieldList<Dimension, DataType>& fieldList,
                             const ConnectivityMap<Dimension>& connectivityMap,
                             const TableKernel<Dimension>& W,
                             const FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldList<Dimension, typename Dimension::SymTensor>& H,
                             const FieldList<Dimension, typename Dimension::Scalar>& weight) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);

  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;

  // Some scratch variables.
  size_t nodeListi, nodeListj;
  Scalar Wi, Wj;
  Vector rij, etai, etaj;

  // Prepare the return value.
  FieldList<Dimension, DataType> result(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> wsum(FieldStorageType::CopyFields);
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    result.appendNewField(fieldList[nodeListi]->name(), fieldList[nodeListi]->nodeList(), DataTypeTraits<DataType>::zero());
    wsum.appendNewField("weight sum", fieldList[nodeListi]->nodeList(), 0.0);
  }

  // Walk all the pairs
  const auto& pairs = connectivityMap.nodePairList();
  for (const auto p: pairs) {
    const auto ki = p.i_list;
    const auto kj = p.j_list;
    const auto i = p.i_node;
    const auto j = p.j_node;

    // State for i
    const auto& ri = position(ki, i);
    const auto  weighti = weight(ki, i);
    const auto& Hi = H(ki, i);
    const auto  Hdeti = Hi.Determinant();

    // State for j
    const auto& rj = position(kj, j);
    const auto  weightj = weight(kj, j);
    const auto& Hj = H(kj, j);
    const auto  Hdetj = Hj.Determinant();

    // Kernel weighting and gradient.
    rij = ri - rj;
    etai = Hi*rij;
    etaj = Hj*rij;
    Wi = W.kernelValue(etai.magnitude(), Hdeti);
    Wj = W.kernelValue(etaj.magnitude(), Hdetj);

    // Sum the pair-wise contributions.
    wsum(ki, i) += weightj*Wj;
    wsum(kj, j) += weighti*Wi;
    result(ki, i) += fieldList(ki, i) * weightj*Wj;
    result(kj, j) += fieldList(kj, j) * weighti*Wi;
  }
  
  // Finalize the interpolation 
  for (auto ki = 0u; ki < numNodeLists; ++ki) {
    const auto n = result[ki]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto  weighti = weight(ki, i);
      const auto& Hi = H(ki, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  wW0 = weighti*W.kernelValue(0.0, Hdeti);
      result(ki, i) = ((result(ki, i) + fieldList(ki, i) * wW0)/
                       (wsum(ki, i) + wW0));
    }
  }

  // That's it.
  return result;
}

}

