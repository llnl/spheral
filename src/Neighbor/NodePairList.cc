#include "NodePairList.hh"
#include "Utilities/DBC.hh"

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

//------------------------------------------------------------------------------
// Data operations
//------------------------------------------------------------------------------

// Warning: Not performant if called frequently
void
NodePairList::push_back(NodePairIdxType nodePair) {
  mNodePairList.push_back(nodePair);
  initMA();
}

void
NodePairList::clear() {
  mNodePairList.clear();
  mPair2Index.clear();
  mData.free();
}

void
NodePairList::reserve(const size_t n) {
  mNodePairList.reserve(n);
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
NodePairList::NodePairList(const NodePairList& rhs) {
  mNodePairList = rhs.mNodePairList;
  initMA();
}

//------------------------------------------------------------------------------
// Assignment constructor
//------------------------------------------------------------------------------
NodePairList& NodePairList::operator=(const NodePairList& rhs) {
  if (this != &rhs) {
    mNodePairList = rhs.mNodePairList;
    initMA();
  }
  return *this;
}

}
