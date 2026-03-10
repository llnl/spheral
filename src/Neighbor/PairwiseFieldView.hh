//---------------------------------Spheral++----------------------------------//
// PairwiseFieldView -- view class for PairwiseField
//
// Stores a value per node pair in a NodePairList. Because connectivity is
// allowed to change step to step in our meshfree methods, PairwiseField is
// ephemeral and will be invalidated when topology is updated.
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_PairwiseFieldView_
#define _Spheral_NeighborSpace_PairwiseFieldView_

#include "Neighbor/PairwiseFieldElementAccessor.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/StrideIterator.hh"
#include "Utilities/DBC.hh"

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

#ifdef SPHERAL_UNIFIED_MEMORY
#include "Utilities/span.hh"
#endif

namespace Spheral {

// Forward declarations
struct NodePairIdxType;

template<typename Value, size_t numElements=1>
class PairwiseFieldView {
public:
  //--------------------------- Public Interface ---------------------------//
#ifdef SPHERAL_UNIFIED_MEMORY
  using ContainerType = SPHERAL_SPAN_TYPE<Value>;
  using value_type = typename ContainerType::value_type;
#else
  using ContainerType = typename chai::ManagedArray<Value>;
  using value_type = Value;
#endif
  
  using SelfType = PairwiseFieldView<Value, numElements>;
  using Accessor = PairwiseFieldDetail::ElementAccessor<SelfType, numElements>;
  using reference = typename Accessor::reference;
  using const_reference = typename Accessor::const_reference;
  using iterator = StrideIterator<Value, numElements>;
  using const_iterator = StrideIterator<const Value, numElements>;

  // Constructors, destructors
  SPHERAL_HOST_DEVICE PairwiseFieldView()                                           = default;
  SPHERAL_HOST_DEVICE PairwiseFieldView(const PairwiseFieldView& rhs)               = default;
  SPHERAL_HOST_DEVICE PairwiseFieldView(PairwiseFieldView&& rhs)                    = default;
  SPHERAL_HOST_DEVICE virtual ~PairwiseFieldView()                                  = default;
  SPHERAL_HOST_DEVICE PairwiseFieldView& operator=(const PairwiseFieldView& rhs)    = default;

  // Access the data
  SPHERAL_HOST_DEVICE reference operator[](const size_t k) const                    { return Accessor::at(mSpan, k); }
  SPHERAL_HOST_DEVICE reference operator()(const size_t k) const                    { return (*this)[k]; }

  SPHERAL_HOST_DEVICE reference operator[](const size_t k)                          { return Accessor::at(mSpan, k); }
  SPHERAL_HOST_DEVICE reference operator()(const size_t k)                          { return (*this)[k]; }

  // Comparators
  SPHERAL_HOST_DEVICE bool operator==(const PairwiseFieldView& rhs) const           { return mSpan == rhs.mSpan; }
  SPHERAL_HOST_DEVICE bool operator!=(const PairwiseFieldView& rhs) const           { return mSpan != rhs.mSpan; }

  // Other methods
  SPHERAL_HOST_DEVICE size_t size() const                                           { return mSpan.size()/numElements; }
#ifdef SPHERAL_UNIFIED_MEMORY
  SPHERAL_HOST_DEVICE bool empty() const                                            { return mSpan.empty(); }
#else
  SPHERAL_HOST        bool empty() const                                            { return mSpan.size() == 0u; }
#endif

  // Zero the Field
  SPHERAL_HOST_DEVICE void Zero()                                                   { for (auto& x: mSpan) x = DataTypeTraits<Value>::zero(); }

protected:
  //--------------------------- Protected Interface ---------------------------//
  ContainerType mSpan;
};

}

#endif

