//---------------------------------Spheral++------------------------------------
// Compute the RK gradient.
//------------------------------------------------------------------------------
#include "RK/gradientRK.hh"
#include "RK/RKUtilities.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"

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

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientRK(const FieldList<Dimension, DataType>& fieldList,
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
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Prepare the result.
  FieldList<Dimension, GradientType> result;
  result.copyFields();
  for (auto fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad ", (*fieldItr)->nodeList()));
  }

  // Walk all the interacting pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {

    // Thread private stuff
    auto result_thread = result.threadCopy();
    int i, j, nodeListi, nodeListj;
    Vector gradWi, gradWj;

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
        // const Scalar wi = fij*2.0*weight(nodeListi, i)*weight(nodeListj, j)/(weight(nodeListi, i) + weight(nodeListj, j));
        // const Scalar wi = fij*0.5*(weight(nodeListi, i) + weight(nodeListj, j));
        // const Scalar wj = wi;
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
        gradWj = WR.evaluateGradient( xij, Hj, correctionsi);
        gradWi = WR.evaluateGradient(-xij, Hi, correctionsj);
        gradFi += wj*Fj*gradWj;
        gradFj += wi*Fi*gradWi;
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
      result(nodeListi, i) += weight(nodeListi, i)*fieldList(nodeListi, i)*WR.evaluateGradient(Vector::zero, Hi, correctionsi);
    }
  }

  // That's it!
  return result;
}

template<typename Dimension, typename DataType>
FieldList<Dimension, std::vector<typename MathTraits<Dimension, DataType>::GradientType>>
gradientRK(const FieldList<Dimension, std::vector<DataType>>& fieldList,
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
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Get size of vector
  const auto vectorSize = (fieldList.numInternalElements() > 0 ?
                           fieldList(fieldList.internalNodeBegin()).size() :
                           0);
  
  // Prepare the result.
  FieldList<Dimension, std::vector<GradientType>> result;
  result.copyFields();
  for (auto fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, std::vector<GradientType>>("grad ",
                                                                   (*fieldItr)->nodeList(),
                                                                   std::vector<GradientType>(vectorSize, GradientType::zero)));
  }
  
  // Walk all the interacting pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {

    // Thread private stuff
    auto result_thread = result.threadCopy();
    int i, j, nodeListi, nodeListj;
    Vector gradWi, gradWj;

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
      CHECK(Fi.size() == vectorSize &&
            gradFi.size() == vectorSize);

      // The coupling between these nodes.
      const auto fij = nodeCoupling(pairs[kk]);
      if (fij > 0.0) {

        // Find the effective weights of i->j and j->i.
        // const Scalar wi = fij*2.0*weight(nodeListi, i)*weight(nodeListj, j)/(weight(nodeListi, i) + weight(nodeListj, j));
        // const Scalar wi = fij*0.5*(weight(nodeListi, i) + weight(nodeListj, j));
        // const Scalar wj = wi;
        const Scalar wi = fij*weight(nodeListi, i);
        const Scalar wj = fij*weight(nodeListj, j);

        // Get the state for node j.
        const auto& xj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto& correctionsj = corrections(nodeListj, j);
        const auto& Fj = fieldList(nodeListj, j);
        auto&       gradFj = result(nodeListj, j);
        CHECK(Fj.size() == vectorSize &&
              gradFj.size() == vectorSize);

        // Pair contributions
        const auto xij = xi - xj;
        gradWj = WR.evaluateGradient( xij, Hj, correctionsi);
        gradWi = WR.evaluateGradient(-xij, Hi, correctionsj);

        for (auto m = 0u; m < vectorSize; ++m) {
          CHECK(gradFi.size() == vectorSize &&
                gradFj.size() == vectorSize &&
                Fj.size() == vectorSize &&
                Fi.size() == vectorSize);
          gradFi[m] += wj*Fj[m]*gradWj;
          gradFj[m] += wi*Fi[m]*gradWi;
        }
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
      const auto wi = weight(nodeListi, i);
      const auto gradWi = WR.evaluateGradient(Vector::zero, Hi, correctionsi);
      const auto& Fi = fieldList(nodeListi, i);
      auto& gradFi = result(nodeListi, i);
      for (auto m = 0u; m < vectorSize; ++m) {
        gradFi[m] += wi*Fi[m]*gradWi;
      }
    }
  }

  // That's it!
  return result;
}

}
