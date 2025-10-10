//------------------------------------------------------------------------------
// Provide dummy OpenMP methods for when we compile without OpenMP support.
//------------------------------------------------------------------------------
#ifndef __Spheral_OpenMP_wrapper__
#define __Spheral_OpenMP_wrapper__

#ifdef _OPENMP
#include "omp.h"
#else
inline int  omp_get_num_threads() { return 1; }
inline int  omp_get_thread_num()  { return 0; }
inline void omp_set_num_threads() {}
#endif

#include <variant>
#include <vector>

//------------------------------------------------------------------------------
// += as vector extension
// Necessary to support reduction operations with FieldList<vector<Vector>>
//------------------------------------------------------------------------------
template<typename T>
inline
std::vector<T>&
operator+=(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

namespace Spheral {

// Forward declarations
template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class FieldList;

//------------------------------------------------------------------------------
// An enum type to help with thread reductions
//------------------------------------------------------------------------------
enum class ThreadReduction {
  MIN = 0,
  MAX = 1,
  SUM = 2
};

// Put the helpers for threadReduceFields in an enclosing struct
template<typename Dimension>
struct SpheralThreads {

  typedef std::variant<FieldList<Dimension, int>*,
                       FieldList<Dimension, typename Dimension::Scalar>*,
                       FieldList<Dimension, typename Dimension::Vector>*,
                       FieldList<Dimension, typename Dimension::Tensor>*,
                       FieldList<Dimension, typename Dimension::SymTensor>*,
                       FieldList<Dimension, typename Dimension::ThirdRankTensor>*,
                       FieldList<Dimension, typename Dimension::FourthRankTensor>*,
                       FieldList<Dimension, typename Dimension::FifthRankTensor>*,
                       FieldList<Dimension, std::vector<typename Dimension::Vector>>*> FieldListReductionVariant;
  typedef std::vector<FieldListReductionVariant> FieldListStack;

};

//------------------------------------------------------------------------------
// Do an appropriate reduction across threads for a set of FieldLists.
// This is identical to calling FieldList::threadReduce, just fusing the
// loops together.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
threadReduceFieldLists(typename SpheralThreads<Dimension>::FieldListStack& stack) {
  if (stack.size() > 0) {
#pragma omp critical
    {
      if (omp_get_num_threads() > 1) {
        const auto& nodeListPtrs = std::visit([](const auto* threadValue) -> const std::vector<NodeList<Dimension>*>& { return threadValue->nodeListPtrs(); }, stack[0]);
        const auto  numNodeLists = nodeListPtrs.size();
        for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
          const auto ni = nodeListPtrs[nodeListi]->numInternalNodes();
          for (auto i = 0u; i < ni; ++i) {
            for (auto& flv: stack) {
              std::visit([=](const auto* threadValue) {
                CHECK(nodeListi < threadValue->size());
                CHECK(i < (*threadValue)[nodeListi]->size());
                switch (threadValue->reductionType) {
                case ThreadReduction::SUM:
                  (*(threadValue->threadMasterPtr))(nodeListi, i) += (*threadValue)(nodeListi,i);
                  break;

                case ThreadReduction::MIN:
                  (*(threadValue->threadMasterPtr))(nodeListi, i) = std::min((*threadValue)(nodeListi, i), (*(threadValue->threadMasterPtr))(nodeListi, i));
                  break;

                case ThreadReduction::MAX:
                  (*(threadValue->threadMasterPtr))(nodeListi, i) = std::max((*threadValue)(nodeListi, i), (*(threadValue->threadMasterPtr))(nodeListi, i));
                }
              }, flv);
            }
          }
        }
      }
    } // OMP crtical
  }
}

}

#endif