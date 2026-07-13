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

#include "Neighbor/PairwiseFieldViewBase.hh"
#include "Neighbor/PairwiseFieldElementAccessor.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/StrideIterator.hh"
#include "Utilities/DBC.hh"

#include "Threading/GPUUtils.hh"

namespace Spheral {

// Forward declarations
struct NodePairIdxType;

template<typename Value, size_t numElements=1>
class PairwiseFieldView:
    public PairwiseFieldViewBase {

public:
  //--------------------------- Public Interface ---------------------------//
#ifdef SPHERAL_UNIFIED_MEMORY
  using ContainerType = SPHERAL_SPAN_TYPE<Value>;
  using value_type = typename ContainerType::value_type;
#else
  using ContainerType = typename chai::ManagedArray<Value>;
  using value_type = Value;
#endif
  
  using SelfType = PairwiseFieldView<value_type, numElements>;
  using Accessor = PairwiseFieldDetail::ElementAccessor<SelfType, numElements>;
  using reference = typename Accessor::reference;
  using const_reference = typename Accessor::const_reference;
  using iterator = StrideIterator<value_type, numElements>;
  using const_iterator = StrideIterator<const value_type, numElements>;

  // Constructors, destructors
  SPHERAL_HOST_DEVICE PairwiseFieldView()                                        = default;
  SPHERAL_HOST_DEVICE PairwiseFieldView(const PairwiseFieldView& rhs)            = default;
  SPHERAL_HOST_DEVICE PairwiseFieldView(PairwiseFieldView&& rhs)                 = default;
  SPHERAL_HOST_DEVICE virtual ~PairwiseFieldView()                               = default;
  SPHERAL_HOST_DEVICE PairwiseFieldView& operator=(const PairwiseFieldView& rhs) = default;

  // Access the data
  SPHERAL_HOST_DEVICE reference operator[](const size_t k) const                 { return Accessor::at(mSpan, k); }
  SPHERAL_HOST_DEVICE reference operator()(const size_t k) const                 { return (*this)[k]; }

  SPHERAL_HOST_DEVICE reference operator[](const size_t k)                       { return Accessor::at(mSpan, k); }
  SPHERAL_HOST_DEVICE reference operator()(const size_t k)                       { return (*this)[k]; }
  SPHERAL_HOST_DEVICE Value* data() const                                        { return mSpan.data(); }

  // Comparators
  SPHERAL_HOST_DEVICE bool operator==(const PairwiseFieldView& rhs) const        { return data() == rhs.data(); }
  SPHERAL_HOST_DEVICE bool operator!=(const PairwiseFieldView& rhs) const        { return data() != rhs.data(); }

  // Other methods
  SPHERAL_HOST_DEVICE size_t size() const                                        { return mSpan.size()/numElements; }
#ifdef SPHERAL_UNIFIED_MEMORY
  SPHERAL_HOST_DEVICE bool empty() const                                         { return mSpan.empty(); }
#else
  SPHERAL_HOST        bool empty() const                                         { return mSpan.size() == 0u; }
#endif

  // Zero the Field
  SPHERAL_HOST_DEVICE void Zero()                                                { for (auto& x: mSpan) x = DataTypeTraits<value_type>::zero(); }

  SPHERAL_HOST virtual void move(chai::ExecutionSpace space) override            { GPUUtils::move(mSpan, space); }

  SPHERAL_HOST virtual void touch(chai::ExecutionSpace space) override           { GPUUtils::touch(mSpan, space); }

protected:
  //--------------------------- Protected Interface ---------------------------//
  ContainerType mSpan;
};

}

#endif

