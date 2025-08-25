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
  initializeMA();
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
// Initialize ManagedArray
//------------------------------------------------------------------------------
template<typename F>
void
NodePairList::initializeMA(F callback) {
  if (!(mNodePairList.data() == mData.data(chai::CPU, false)
        && mNodePairList.size() == mData.size())) {
    mData.free();
    mData = chai::makeManagedArray(mNodePairList.data(), mNodePairList.size(), chai::CPU, false);
    mData.setUserCallback(callback);
  }
}

void
NodePairList::initializeMA() {
  this->initializeMA([](const chai::PointerRecord*,
                        chai::Action action,
                        chai::ExecutionSpace) { });
}

//------------------------------------------------------------------------------
// Initialize ManagedArray
//------------------------------------------------------------------------------
template<typename F>
NodePairListView NodePairList::view(F callback) {
  initializeMA(callback);
  return NodePairListView(mData);
}

NodePairListView NodePairList::view() {
  initializeMA();
  return NodePairListView(mData);
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
NodePairList::NodePairList(const NodePairList& rhs)
  :
  NodePairListView(rhs) {
  mNodePairList = rhs.mNodePairList;
  initializeMA();
}

//------------------------------------------------------------------------------
// Assignment constructor
//------------------------------------------------------------------------------
NodePairList& NodePairList::operator=(const NodePairList& rhs) {
  if (this !=&rhs) {
    NodePairListView::operator=(rhs);
    mNodePairList = rhs.mNodePairList;
    initializeMA();
  }
  return *this;
}
}
