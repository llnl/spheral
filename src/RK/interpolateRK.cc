//---------------------------------Spheral++------------------------------------
// Compute the RK interpolate.
//------------------------------------------------------------------------------
#include "RK/interpolateRK.hh"
#include "RK/RKUtilities.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/Timer.hh"
#include <variant>

namespace Spheral {

using std::variant;
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

namespace {
//------------------------------------------------------------------------------
// Verify FieldList sizes.
//------------------------------------------------------------------------------
struct FieldListSize {
  template<typename FieldListType>
  inline
  int operator()(const FieldListType& x) const {
    return x.size();
  }
};

template<typename Dimension> 
inline
bool checkSizes(const vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
                                     FieldList<Dimension, typename Dimension::Vector>,
                                     FieldList<Dimension, typename Dimension::Tensor>,
                                     FieldList<Dimension, typename Dimension::SymTensor>,
                                     FieldList<Dimension, typename Dimension::ThirdRankTensor>>>& fieldLists,
                const size_t numNodeLists) {
  for (const auto& f: fieldLists) {
    if (std::visit(FieldListSize(), f) != (int)numNodeLists) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// copyFields
//------------------------------------------------------------------------------
struct CopyFields {
  template<typename FieldListType>
  inline
  void operator()(FieldListType& x) const { x.copyFields(); }
};

//------------------------------------------------------------------------------
// Zero
//------------------------------------------------------------------------------
struct ZeroFields {
  template<typename FieldListType>
  inline
  void operator()(FieldListType& x) const { x.Zero(); }
};

//------------------------------------------------------------------------------
// PrependNameFields
//------------------------------------------------------------------------------
struct PrependNameFields {
  std::string prefix;
  
  PrependNameFields(const std::string prefix_):
    prefix(prefix_) {}

  template<typename FieldListType>
  inline
  void operator()(FieldListType& x) const {
    for (auto fieldPtr: x) fieldPtr->name(prefix + fieldPtr->name());
  }
};

//------------------------------------------------------------------------------
// IncrementElement
//------------------------------------------------------------------------------
struct IncrementElement {
  int nodeListi, i, nodeListj, j;
  double mult;

  IncrementElement(int nodeListi_, int i_, int nodeListj_, int j_, double mult_):
    nodeListi(nodeListi_),
    i(i_),
    nodeListj(nodeListj_),
    j(j_),
    mult(mult_) {}

  template<typename A, typename B>
  inline
  void operator()(A&, const B&) const { VERIFY(false); }

  template<typename FieldListType>
  inline
  void operator()(FieldListType& x, const FieldListType& y) const {
    CHECK2(nodeListi < (int)x.size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    CHECK2(i < (int)x[nodeListi]->size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    CHECK2(nodeListj < (int)y.size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    CHECK2(j < (int)y[nodeListj]->size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    x(nodeListi, i) += mult*y(nodeListj, j);
  }
};

//------------------------------------------------------------------------------
// IncrementFieldList
//------------------------------------------------------------------------------
struct IncrementFieldList {

  template<typename A, typename B>
  inline
  void operator()(A&, const B&) const { VERIFY(false); }

  template<typename FieldListType>
  inline
  void operator()(FieldListType& x, const FieldListType& y) const {
    CHECK(x.size() == y.size());
    x += y;
  }
};

//------------------------------------------------------------------------------
// Work around the templating of the correction order
//------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
// interplateRK
//------------------------------------------------------------------------------
template<typename Dimension>
vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
               FieldList<Dimension, typename Dimension::Vector>,
               FieldList<Dimension, typename Dimension::Tensor>,
               FieldList<Dimension, typename Dimension::SymTensor>,
               FieldList<Dimension, typename Dimension::ThirdRankTensor>>>
interpolateRK(const vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
                                   FieldList<Dimension, typename Dimension::Vector>,
                                   FieldList<Dimension, typename Dimension::Tensor>,
                                   FieldList<Dimension, typename Dimension::SymTensor>,
                                   FieldList<Dimension, typename Dimension::ThirdRankTensor>>>& fieldLists,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  const ReproducingKernel<Dimension>& WR,
                  const FieldList<Dimension, RKCoefficients<Dimension>>& corrections,
                  const NodeCoupling& nodeCoupling) {

  TIME_FUNCTION;

  typedef typename Dimension::Vector Vector;
  typedef vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
                         FieldList<Dimension, typename Dimension::Vector>,
                         FieldList<Dimension, typename Dimension::Tensor>,
                         FieldList<Dimension, typename Dimension::SymTensor>,
                         FieldList<Dimension, typename Dimension::ThirdRankTensor>>> FieldListArray;

  // Pre-conditions.
  const auto numNodeLists = position.size();
  const auto numFieldLists = fieldLists.size();
  REQUIRE(checkSizes<Dimension>(fieldLists, numNodeLists));
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(corrections.size() == numNodeLists);

  // Prepare the result.
  FieldListArray result;
  for (const auto& fieldList: fieldLists) {
    result.push_back(fieldList);
    std::visit(CopyFields(), result.back());
    std::visit(ZeroFields(), result.back());
    std::visit(PrependNameFields("interpolate "), result.back());
  }

  // Walk all the interacting pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {

    // Thread private result
    FieldListArray localResult;
    #pragma omp critical
    for (const auto& fieldList: fieldLists) {
      localResult.push_back(fieldList);
      std::visit(CopyFields(), localResult.back());
      std::visit(ZeroFields(), localResult.back());
      std::visit(PrependNameFields("local interpolate "), localResult.back());
    }

    int i, j, nodeListi, nodeListj;

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);

      // The coupling between these nodes.
      const auto fij = nodeCoupling(pairs[kk]);
      if (fij > 0.0) {

        // Find the effective weights of i->j and j->i.
        const auto wi = fij*weight(nodeListi, i);
        const auto wj = fij*weight(nodeListj, j);

        // Get the state for node j.
        const auto& rj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto& correctionsj = corrections(nodeListj, j);

        // Kernel weight.
        const auto rij = ri - rj;
        const auto Wj = WR.evaluateKernel( rij, Hj, correctionsi);
        const auto Wi = WR.evaluateKernel(-rij, Hi, correctionsj);

        // Increment the pair-wise values.
        for (auto k = 0u; k < numFieldLists; ++k) {
          std::visit(IncrementElement(nodeListi, i, nodeListj, j, wj*Wj), localResult[k], fieldLists[k]);
          std::visit(IncrementElement(nodeListj, j, nodeListi, i, wi*Wi), localResult[k], fieldLists[k]);
        }
      }
    }

    // Merge the local to global result
#pragma omp critical
    {
      for (auto k = 0u; k < numFieldLists; ++k) {
        std::visit(IncrementFieldList(), result[k], localResult[k]);
      }
    }
  }

  // Add the self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = position[nodeListi]->nodeList().numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {

      // Get the state for node i.
      const auto& Hi = H(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      const auto  Wj = WR.evaluateKernel(Vector::zero, Hi, correctionsi);

      // Add the self-contribution to each FieldList.
      for (auto k = 0u; k < numFieldLists; ++k) {
        std::visit(IncrementElement(nodeListi, i, nodeListi, i, weight(nodeListi, i)*Wj), result[k], fieldLists[k]);
      }
    }
  }


  // That's it!
  return result;
}

}
