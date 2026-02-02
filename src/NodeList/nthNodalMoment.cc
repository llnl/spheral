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

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MomentTraits<Dimension, moment>::Moment Moment;

  // The total number of NodeLists we're working on.
  const size_t numNodeLists = distance(nodeListBegin, nodeListEnd);

  // Build a connectivity map for walking nodes.  This relies on the NodeLists 
  // Neighbor objects being up to date.
  const ConnectivityMap<Dimension> cm(nodeListBegin, nodeListEnd, false, false, false);

  // Build up the FieldLists of positions, H's, and the first moment that we're going
  // to build.
  FieldList<Dimension, Vector> pos(FieldStorageType::ReferenceFields);
  FieldList<Dimension, SymTensor> H(FieldStorageType::ReferenceFields);
  FieldList<Dimension, Moment> result(FieldStorageType::CopyFields);
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr) {
    const NodeList<Dimension>& nodes = **itr;
    pos.appendField(nodes.positions());
    H.appendField(nodes.Hfield());
    result.appendNewField("moment", nodes, DataTypeTraits<Moment>::zero());
  }

  // Find the moment of the node distribution in eta coordinates.
  const double W0 = W(0.0, 1.0);
  unsigned nodeListi = 0;
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodes = **itr;
    for (unsigned i = 0; i != nodes.numInternalNodes(); ++i) {
      const vector<vector<int> >& allNeighbors = cm.connectivityForNode(nodeListi, i);
      CHECK(allNeighbors.size() == numNodeLists);
      const Vector ri = pos(nodeListi, i);
      const SymTensor Hi = H(nodeListi, i);
      double wsum = W0;
      for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& neighbors = allNeighbors[nodeListj];
        for (unsigned k = 0; k != neighbors.size(); ++k) {
          const unsigned j = neighbors[k];
          const Vector etai = Hi*(pos(nodeListj, j) - ri);
          const double Wi = W(etai.magnitude(), 1.0);
          wsum += Wi;
          result(nodeListi, i) += Wi * nthMomentKernel<Dimension, moment>()(etai);
        }
      }
      if (renormalize) result(nodeListi, i) *= safeInv(wsum);
    }
  }

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
  unsigned nodeListi = 0;
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodes = **itr;
    for (unsigned i = 0; i != nodes.numInternalNodes(); ++i) {
      const vector<vector<int> >& allNeighbors = cm.connectivityForNode(nodeListi, i);
      CHECK(allNeighbors.size() == numNodeLists);
      const Vector ri = pos(nodeListi, i);
      const SymTensor Hi = H(nodeListi, i);
      for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& neighbors = allNeighbors[nodeListj];
        for (unsigned k = 0; k != neighbors.size(); ++k) {
          const unsigned j = neighbors[k];
          const Vector rj = pos(nodeListj, j);
          const Vector etai = Hi*(rj - ri);
          const double Wi = useGradientAsKernel ? abs(W.gradValue(etai.magnitude(), 1.0)) : W.kernelValue(etai.magnitude(), 1.0);
          zerothMoment(nodeListi, i) += Wi;
          firstMoment(nodeListi, i) += Wi * etai;
        }
      }
      firstMoment(nodeListi, i) *= safeInv(zerothMoment(nodeListi, i));
      zerothMoment(nodeListi, i) = Dimension::rootnu(zerothMoment(nodeListi, i));
    }
  }
}

}
