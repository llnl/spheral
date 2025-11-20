#include "NodeList/NodeList.hh"
#include "Field/FieldBase.hh"

#include <algorithm>
#include <sstream>

namespace Spheral {

//------------------------------------------------------------------------------
// The number of NodeLists currently registered.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NodeListRegistrar<Dimension>::
numNodeLists() const {
  return mNodeLists.size();
}

template<typename Dimension>
inline
int
NodeListRegistrar<Dimension>::
numFluidNodeLists() const {
  return mFluidNodeLists.size();
}

//------------------------------------------------------------------------------
// Iterators over the registered NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::iterator
NodeListRegistrar<Dimension>::
begin() {
  return mNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::iterator
NodeListRegistrar<Dimension>::
end() {
  return mNodeLists.end();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_iterator
NodeListRegistrar<Dimension>::
begin() const {
  return mNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_iterator
NodeListRegistrar<Dimension>::
end() const {
  return mNodeLists.end();
}

// FluidNodeLists.
template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::fluid_iterator
NodeListRegistrar<Dimension>::
fluidBegin() {
  return mFluidNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::fluid_iterator
NodeListRegistrar<Dimension>::
fluidEnd() {
  return mFluidNodeLists.end();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_fluid_iterator
NodeListRegistrar<Dimension>::
fluidBegin() const {
  return mFluidNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_fluid_iterator
NodeListRegistrar<Dimension>::
fluidEnd() const {
  return mFluidNodeLists.end();
}

//------------------------------------------------------------------------------
// Generic method to find the proper place in a sequence to insert a Field
// or NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename IteratorType, typename ThingyType>
inline
IteratorType
NodeListRegistrar<Dimension>::
findInsertionPoint(const ThingyType& thingy,
                   const IteratorType begin,
                   const IteratorType end) {

  // If the input iterator sequence is empty, then the answer is easy!
  const auto containerSize = std::distance(begin, end);
  if (containerSize == 0) return end;

  // Get the sequence of NodeLists represented by the input.
  std::vector<NodeList<Dimension>*> nodeListPtrs;
  nodeListPtrs.reserve(containerSize);
  for (auto itr = begin; itr != end; ++itr) {
    auto nodeListPtr = getNodeListPtr(*itr);
    CHECK(itr == begin or (nodeListPtr->name() > getNodeListPtr(*(itr - 1))->name()));
    nodeListPtrs.push_back(nodeListPtr);
  }
  CHECK((int)nodeListPtrs.size() == containerSize);

  // Now we can find where the specified thingy should go.
  auto nodeListPtr = getNodeListPtr(thingy);
  auto orderItr = std::upper_bound(nodeListPtrs.begin(),
                                   nodeListPtrs.end(),
                                   nodeListPtr,
                                   NodeListComparator());
  const auto displacement = std::distance(nodeListPtrs.begin(), orderItr);
  CHECK(displacement >= 0 && displacement <= containerSize);
  auto result = begin + displacement;
  ENSURE(result >= begin && result <= end);
  ENSURE((result > begin and getNodeListPtr(thingy)->name() > getNodeListPtr(*(result - 1))->name()) or
                             getNodeListPtr(thingy)->name() < getNodeListPtr(*begin)->name());
  return result;
}

//------------------------------------------------------------------------------
// Similar method to directly sort a set of stuff into NodeList ordering.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename IteratorType>
void
NodeListRegistrar<Dimension>::
sortInNodeListOrder(IteratorType begin, IteratorType end) {
  std::sort(begin, end, NodeListRegistrar<Dimension>::NodeListComparator());
  // std::sort(begin, end, [](const IteratorType a, const IteratorType b, NodeListRegistrar<Dimension>::NodeListComparator())
}

}
