#include "NodePairList.hh"
#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {
  
//------------------------------------------------------------------------------
// index
//------------------------------------------------------------------------------
size_t
NodePairList::index(const NodePairIdxType& x) const {
  if (mPair2Index.size() != mNodePairList.size()) computeLookup();  // Lazy evaluation
  auto itr = mPair2Index.find(x);
  CHECK(itr != mPair2Index.end());
  return itr->second;
}

//------------------------------------------------------------------------------
// Recompute the lookup table for NodePair->index
//------------------------------------------------------------------------------
void
NodePairList::computeLookup() const {
  mPair2Index.clear();
  const auto n = this->size();
  for (size_t k = 0u; k < n; ++k) {
    mPair2Index[mNodePairList[k]] = k;
  }
}

  void NodePairList::makeUnique() {
    // Make sure the node pairs are ordered correctly
    for (auto kk = 0u; kk < mNodePairList.size(); ++kk)
    {
      auto& pair = mNodePairList[kk];
      if (pair.i_list > pair.j_list)
      {
        const auto temp_list = pair.j_list;
        const auto temp_node = pair.j_node;
        pair.j_list = pair.i_list;
        pair.j_node = pair.i_node;
        pair.i_list = temp_list;
        pair.i_node = temp_node;
      }
      if (pair.i_list == pair.j_list && pair.i_node > pair.j_node)
      {
        const auto temp = pair.j_node;
        pair.j_node = pair.i_node;
        pair.i_node = temp;
      }
    }
    
    // Remove duplicates
    std::sort(mNodePairList.begin(), mNodePairList.end());
    mNodePairList.erase(std::unique(mNodePairList.begin(), mNodePairList.end()), mNodePairList.end());
  }

}
