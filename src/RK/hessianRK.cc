//---------------------------------Spheral++------------------------------------
// Compute the RK hessian.
//------------------------------------------------------------------------------
#include "RK/hessianRK.hh"
#include "RK/RKUtilities.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Geometry/outerProduct.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::HessianType>
hessianRK(const FieldList<Dimension, DataType>& fieldList,
          const FieldList<Dimension, typename Dimension::Vector>& position,
          const FieldList<Dimension, typename Dimension::Scalar>& weight,
          const FieldList<Dimension, typename Dimension::SymTensor>& H,
          const ConnectivityMap<Dimension>& connectivityMap,
          const ReproducingKernel<Dimension>& WR,
          const FieldList<Dimension, RKCoefficients<Dimension>>& corrections,
          const NodeCoupling& nodeCoupling) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(corrections.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::HessianType HessianType;

  // Prepare the result.
  FieldList<Dimension, HessianType> result;
  result.copyFields();
  for (auto fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, HessianType>("hessian ", (*fieldItr)->nodeList()));
  }

  // Walk all the interacting pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {

    // Thread private stuff
    auto result_thread = result.threadCopy();
    int i, j, nodeListi, nodeListj;
    SymTensor hessWi, hessWj;

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& xi = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      const auto& Fi = fieldList(nodeListi, i);
      auto&       gradFi = result_thread(nodeListi, i);

      // The coupling between these nodes.
      const auto fij = nodeCoupling(pairs[kk]);
      if (fij > 0.0) {

        // Find the effective weights of i->j and j->i.
        const Scalar wi = fij*weight(nodeListi, i);
        const Scalar wj = fij*weight(nodeListj, j);

        // Get the state for node j.
        const auto& xj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto& correctionsj = corrections(nodeListj, j);
        const auto& Fj = fieldList(nodeListj, j);
        auto&       gradFj = result(nodeListj, j);

        // Pair contributions
        const auto xij = xi - xj;
        hessWj = WR.evaluateHessian( xij, Hj, correctionsi);
        hessWi = WR.evaluateHessian(-xij, Hi, correctionsj);
        gradFi += wj*outerProduct<Dimension>(hessWj, Fj);
        gradFj += wi*outerProduct<Dimension>(hessWi, Fi);
      }
    }

#pragma omp critical
    {
      result_thread.threadReduce();
    }
  }

  // Add the self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = position[nodeListi]->nodeList().numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto& Hi = H(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      result(nodeListi, i) += weight(nodeListi, i)*outerProduct<Dimension>(WR.evaluateHessian(Vector::zero(), Hi, correctionsi), fieldList(nodeListi, i));
    }
  }

  // That's it!
  return result;
}

}
