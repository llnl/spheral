//---------------------------------Spheral++----------------------------------//
// FieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "FieldListFunctions.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

using std::vector;

//------------------------------------------------------------------------------
// Return a smoothed estimate of the given FieldList.
//------------------------------------------------------------------------------
namespace Spheral {

template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFields(const FieldList<Dimension, DataType>& fieldList,
             const FieldList<Dimension, typename Dimension::Vector>& position,
             const FieldList<Dimension, typename Dimension::Scalar>& weight,
             const FieldList<Dimension, typename Dimension::Scalar>& /*mass*/,
             const FieldList<Dimension, typename Dimension::Scalar>& /*rho*/,
             const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
             const TableKernel<Dimension>& kernel) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;

  const auto numNodeLists = fieldList.numFields();
  REQUIRE(position.numFields() == numNodeLists);
  REQUIRE(weight.numFields() == numNodeLists);
  REQUIRE(Hfield.numFields() == numNodeLists);

  // Return FieldList.
  FieldList<Dimension, DataType> result;
  result.copyFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto* fieldPtr = fieldList[nodeListi];
    REQUIRE(position[nodeListi]->nodeListPtr() == fieldPtr->nodeListPtr());
    REQUIRE(weight[nodeListi]->nodeListPtr() == fieldPtr->nodeListPtr());
    REQUIRE(Hfield[nodeListi]->nodeListPtr() == fieldPtr->nodeListPtr());
    result.appendField(Field<Dimension, DataType>("smooth " + fieldPtr->name(), fieldPtr->nodeList()));
  }

  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto& connectivityMap = connectivity.connectivityMap();
  const auto searchType = connectivity.searchType();

  // NodePairList does not contain self-pairs, so initialize the self
  // contribution for every internal node.
  for (auto fieldi = 0u; fieldi < numNodeLists; ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& Hi = Hfield(fieldi, i);
      result(fieldi, i) += weight(fieldi, i)*fieldList(fieldi, i)*kernel(0.0, Hi.Determinant());
    }
  }

  // Visit each interacting pair once and accumulate its contribution in both
  // directions.  Ghost endpoints contribute to internal endpoints, but are
  // not themselves smoothed when ghost connectivity was not requested.
  for (const auto& pair: connectivityMap.nodePairList()) {
    REQUIRE(pair.i_list < numNodeLists and pair.j_list < numNodeLists);
    const auto fieldi = connectivity.fieldIndex(pair.i_list);
    const auto fieldj = connectivity.fieldIndex(pair.j_list);
    REQUIRE(fieldi < numNodeLists and fieldj < numNodeLists);
    const auto i = pair.i_node;
    const auto j = pair.j_node;

    const auto& ri = position(fieldi, i);
    const auto& rj = position(fieldj, j);
    const auto& Hi = Hfield(fieldi, i);
    const auto& Hj = Hfield(fieldj, j);
    const auto rij = ri - rj;
    const auto Wi = kernel((Hi*rij).magnitude(), Hi.Determinant());
    const auto Wj = kernel((Hj*rij).magnitude(), Hj.Determinant());

    Scalar Wij = 0.0, Wji = 0.0;
    switch(searchType) {
    case NeighborSearchType::GatherScatter:
      Wij = Wji = 0.5*(Wi + Wj);
      break;
    case NeighborSearchType::Gather:
      Wij = Wi;
      Wji = Wj;
      break;
    case NeighborSearchType::Scatter:
      Wij = Wj;
      Wji = Wi;
      break;
    default:
      VERIFY2(false, "Unhandled neighbor search type.");
    }

    if (i < fieldList[fieldi]->nodeListPtr()->numInternalNodes()) {
      result(fieldi, i) += weight(fieldj, j)*fieldList(fieldj, j)*Wij;
    }
    if (j < fieldList[fieldj]->nodeListPtr()->numInternalNodes()) {
      result(fieldj, j) += weight(fieldi, i)*fieldList(fieldi, i)*Wji;
    }
  }

  return result;
}

}
