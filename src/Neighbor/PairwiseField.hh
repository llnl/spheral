//---------------------------------Spheral++----------------------------------//
// PairwiseField
//
// Stores a value per node pair in a NodePairList. Because connectivity is
// allowed to change step to step in our meshfree methods, PairwiseField is
// ephemeral and will be invalidated when topology is updated.
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_PairwiseField_
#define _Spheral_NeighborSpace_PairwiseField_

#include "Neighbor/PairwiseFieldView.hh"

#include <memory>

namespace Spheral {

// Forward declarations
template<typename Dimension> class ConnectivityMap;
struct NodePairIdxType;
class NodePairList;

template<typename Dimension, typename Value, size_t numElements=1>
class PairwiseField: public PairwiseFieldView<Value, numElements> {
public:
  //--------------------------- Public Interface ---------------------------//
  using ViewType = PairwiseFieldView<Value, numElements>;

  using ContainerType = std::vector<Value>;
  using value_type = typename ContainerType::value_type;
  using SelfType = PairwiseField<Dimension, Value, numElements>;
  using Accessor = PairwiseFieldDetail::ElementAccessor<SelfType, numElements>;
  using reference = typename Accessor::reference;
  using const_reference = typename Accessor::const_reference;
  using iterator = StrideIterator<Value, numElements>;
  using const_iterator = StrideIterator<const Value, numElements>;

  // Bring in various methods hidden in PairwiseFieldView
  using ViewType::operator();
  using ViewType::operator[];

  // Constructors, destructors
  PairwiseField()                                               = default;
  PairwiseField(const ConnectivityMap<Dimension>& connectivity);
  PairwiseField(std::shared_ptr<NodePairList> pairsPtr);
  PairwiseField(const PairwiseField& rhs);
  virtual ~PairwiseField();
  PairwiseField& operator=(const PairwiseField& rhs);

  // Access the data using Node pair information
  const_reference operator()(const NodePairIdxType& x) const;
  reference       operator()(const NodePairIdxType& x);

  // Iterators
  const_iterator begin() const                                  { REQUIRE(!mPairsPtr.expired()); return const_iterator(&(*mArray.begin())); }
  const_iterator end() const                                    { REQUIRE(!mPairsPtr.expired()); return const_iterator(&(*mArray.end())); }

  iterator begin()                                              { REQUIRE(!mPairsPtr.expired()); return iterator(&(*mArray.begin())); }
  iterator end()                                                { REQUIRE(!mPairsPtr.expired()); return iterator(&(*mArray.end())); }

  // Comparators
  bool operator==(const PairwiseField& rhs) const               { REQUIRE(!mPairPtr.expired()); return mArray == rhs.mArray; }
  bool operator!=(const PairwiseField& rhs) const               { REQUIRE(!mPairPtr.expired()); return mArray != rhs.mArray; }

  // Other methods
  const NodePairList& pairs() const;
  
  // Get the view (for trivially copyable types)
  ViewType view()                                               { return static_cast<ViewType>(*this); }

private:
  //--------------------------- Private Interface ---------------------------//
  using ViewType::mSpan;
  std::weak_ptr<NodePairList> mPairsPtr;
  ContainerType mArray;

  void assignDataSpan();
};

}

#include "Neighbor/PairwiseFieldInline.hh"

#endif

