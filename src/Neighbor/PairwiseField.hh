//---------------------------------Spheral++----------------------------------//
// PairwiseField
//
// Stores a value per node pair in a NodePairList. Because connectivity is
// allowed to change step to step in our meshfree methods, PairwiseField is
// ephemeral and will be invalidated when topology is updated.
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_PairwiseField_hh_
#define _Spheral_NeighborSpace_PairwiseField_hh_

#include "Neighbor/PairwiseFieldElementAccessor.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/StrideIterator.hh"
#include "Utilities/DBC.hh"

#include <vector>
#include <memory>

namespace Spheral {

// Forward declarations
template<typename Dimension> class ConnectivityMap;
struct NodePairIdxType;
class NodePairList;

template<typename Dimension, typename Value, size_t numElements=1>
class PairwiseField {
public:
  //--------------------------- Public Interface ---------------------------//
  using ContainerType = std::vector<Value>;
  using value_type = typename ContainerType::value_type;
  using SelfType = PairwiseField<Dimension, Value, numElements>;
  using Accessor = PairwiseFieldDetail::ElementAccessor<SelfType, numElements>;
  using reference = typename Accessor::reference;
  using const_reference = typename Accessor::const_reference;
  using iterator = StrideIterator<Value, numElements>;
  using const_iterator = StrideIterator<const Value, numElements>;

  // Constructors, destructors
  PairwiseField(const ConnectivityMap<Dimension>& connectivity);
  PairwiseField(const PairwiseField& rhs)                       = default;
  ~PairwiseField()                                              = default;
  PairwiseField& operator=(const PairwiseField& rhs)            = default;

  // Access the data
  const_reference operator[](const size_t k) const              { REQUIRE(!mPairsPtr.expired()); return Accessor::at(mValues, k); }
  const_reference operator()(const size_t k) const              { return (*this)[k]; }
  const_reference operator()(const NodePairIdxType& x) const;

  reference       operator[](const size_t k)                    { REQUIRE(!mPairsPtr.expired()); return Accessor::at(mValues, k); }
  reference       operator()(const size_t k)                    { return (*this)[k]; }
  reference       operator()(const NodePairIdxType& x);

  // Comparators
  bool operator==(const PairwiseField& rhs) const               { REQUIRE(!mPairsPtr.expired()); return mValues == rhs.mValues; }
  bool operator!=(const PairwiseField& rhs) const               { REQUIRE(!mPairsPtr.expired()); return mValues != rhs.mValues; }

  // Iterators
  const_iterator begin() const                                  { REQUIRE(!mPairsPtr.expired()); return const_iterator(&(*mValues.begin())); }
  const_iterator end() const                                    { REQUIRE(!mPairsPtr.expired()); return const_iterator(&(*mValues.end())); }

  iterator begin()                                              { REQUIRE(!mPairsPtr.expired()); return iterator(&(*mValues.begin())); }
  iterator end()                                                { REQUIRE(!mPairsPtr.expired()); return iterator(&(*mValues.end())); }

  // Other methods
  const NodePairList& pairs() const;
  size_t size() const                                           { REQUIRE(!mPairsPtr.expired()); return mValues.size()/numElements; }

  // Zero the Field
  void Zero()                                                   { for (auto& x: mValues) x = DataTypeTraits<Value>::zero(); }

  // Forbidden methods
  PairwiseField()                                               = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  std::weak_ptr<NodePairList> mPairsPtr;
  ContainerType mValues;
};

}

#include "Neighbor/PairwiseFieldInline.hh"

#endif

