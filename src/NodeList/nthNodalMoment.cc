//---------------------------------Spheral++----------------------------------//
// nthNodalMoment
// nodalMoments
//
// Compute the nth moment of the local nodal distribution in \eta space:
//
//    \sum_j  (\eta_i)^n W_ij
//    -----------------------
//         \sum_j W_ij
//
// Nodal moments is specialized to simultaneously compute the non-normalized 
// zeroth and normalized first moments.
//
// Created by JMO, Mon May  9 16:20:18 PDT 2011
//----------------------------------------------------------------------------//
#include "nthNodalMoment.hh"
#include "Geometry/Dimension.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Specialized method to compute the pair-wise contribution to the moment.
//------------------------------------------------------------------------------
template<typename Dimension, unsigned moment> struct nthMomentKernel;
template<typename Dimension> struct nthMomentKernel<Dimension, 0U> {
  typename Dimension::Scalar operator()(const typename Dimension::Vector& /*eta*/) { return 1.0; }
};
template<typename Dimension> struct nthMomentKernel<Dimension, 1U> {
  typename Dimension::Vector operator()(const typename Dimension::Vector& eta) { return eta; }
};
template<typename Dimension> struct nthMomentKernel<Dimension, 2U> {
  typename Dimension::SymTensor operator()(const typename Dimension::Vector& eta) { return eta.selfdyad(); }
};

//------------------------------------------------------------------------------
// Generalized moment.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator, unsigned moment>
FieldList<Dimension, typename MomentTraits<Dimension, moment>::Moment>
nthNodalMoment(const NodeListIterator nodeListBegin,
               const NodeListIterator nodeListEnd,
               const TableKernel<Dimension>& W,
               const bool renormalize) {

  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using Moment = typename MomentTraits<Dimension, moment>::Moment;

  // Build a connectivity map for walking nodes.  This relies on the NodeLists 
  // Neighbor objects being up to date.
  const ConnectivityMap<Dimension> cm(nodeListBegin, nodeListEnd, false, false, false);
  const auto& pairs = cm.nodePairList();

  // Build up the FieldLists of positions, H's, and the first moment that we're going
  // to build.
  FieldList<Dimension, Vector> pos(FieldStorageType::ReferenceFields);
  FieldList<Dimension, SymTensor> H(FieldStorageType::ReferenceFields);
  FieldList<Dimension, Scalar> wsum(FieldStorageType::CopyFields);
  FieldList<Dimension, Moment> result(FieldStorageType::CopyFields);
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr) {
    const NodeList<Dimension>& nodes = **itr;
    pos.appendField(nodes.positions());
    H.appendField(nodes.Hfield());
    wsum.appendNewField("wsum", nodes, W(0.0, 1.0));
    result.appendNewField("moment", nodes, DataTypeTraits<Moment>::zero());
  }

  // Find the moment of the node distribution in eta coordinates.
  for (const auto& p: pairs) {
    const auto ki = p.i_list;
    const auto kj = p.j_list;
    const auto i = p.i_node;
    const auto j = p.j_node;

    const auto& ri = pos(ki, i);
    const auto& Hi = H(ki, i);
    const auto& rj = pos(kj, j);
    const auto& Hj = H(kj, j);

    const auto etai = Hi*(rj - ri);
    const auto Wi = W(etai.magnitude(), 1.0);
    wsum(ki, i) += Wi;
    result(ki, i) += Wi * nthMomentKernel<Dimension, moment>()(etai);

    const auto etaj = Hj*(ri - rj);
    const auto Wj = W(etaj.magnitude(), 1.0);
    wsum(kj, j) += Wj;
    result(kj, j) += Wj * nthMomentKernel<Dimension, moment>()(etaj);
  }

  if (renormalize) result /= wsum;

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Zeroth and first moment.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator>
void
zerothAndFirstNodalMoments(const NodeListIterator nodeListBegin,
                           const NodeListIterator nodeListEnd,
                           const TableKernel<Dimension>& W,
                           const bool useGradientAsKernel,
                           FieldList<Dimension, typename Dimension::Scalar>& zerothMoment,
                           FieldList<Dimension, typename Dimension::Vector>& firstMoment) {

  // Preconditions.
  VERIFY(zerothMoment.numFields() == 0);
  VERIFY(firstMoment.numFields() == 0);

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // The total number of NodeLists we're working on.
  const size_t numNodeLists = distance(nodeListBegin, nodeListEnd);

  // Build a connectivity map for walking nodes.  This relies on the NodeLists 
  // Neighbor objects being up to date.
  const ConnectivityMap<Dimension> cm(nodeListBegin, nodeListEnd, false, false, false);
  const auto& pairs = cm.nodePairList();

  // Value of the kernel at the center.
  const double W0 = 0.0; // useGradientAsKernel ?  abs(W.gradValue(0.0, 1.0)) : W.kernelValue(0.0, 1.0);

  // Build up the FieldLists of positions, H's, and the moments that we're going
  // to build.
  FieldList<Dimension, Vector> pos(FieldStorageType::ReferenceFields);
  FieldList<Dimension, SymTensor> H(FieldStorageType::ReferenceFields);
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr) {
    const NodeList<Dimension>& nodes = **itr;
    pos.appendField(nodes.positions());
    H.appendField(nodes.Hfield());
    zerothMoment.appendNewField("zeroth moment", nodes, W0);
    firstMoment.appendNewField("first moment", nodes, Vector::zero());
  }

  // Find the moment of the node distribution in eta coordinates.
  for (const auto& p: pairs) {
    const auto ki = p.i_list;
    const auto kj = p.j_list;
    const auto i = p.i_node;
    const auto j = p.j_node;

    const auto& ri = pos(ki, i);
    const auto& Hi = H(ki, i);
    const auto& rj = pos(kj, j);
    const auto& Hj = H(kj, j);

    const auto etai = Hi*(rj - ri);
    const auto Wi = W(etai.magnitude(), 1.0);
    zerothMoment(ki, i) += Wi;
    firstMoment(ki, i) += Wi * etai;

    const auto etaj = Hj*(ri - rj);
    const auto Wj = W(etaj.magnitude(), 1.0);
    zerothMoment(kj, j) += Wj;
    firstMoment(kj, j) += Wj * etaj;
  }

  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = zerothMoment[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      firstMoment(k, i) *= safeInv(zerothMoment(k, i));
      zerothMoment(k, i) = Dimension::rootnu(zerothMoment(k, i));
    }
  }
}

}
