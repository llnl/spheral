//---------------------------------Spheral++----------------------------------//
// PairwiseFieldElement Accessor
//
// Allows us to access either single element values without the index 
// operator (by reference), or multiple value elements requiring an index
// operator (by pointer)
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_PairwiseFieldElementAccessor__
#define _Spheral_NeighborSpace_PairwiseFieldElementAccessor__

#include "Utilities/DBC.hh"
#include "chai/ManagedArray.hpp"

namespace Spheral {
namespace PairwiseFieldDetail {

// Trait class to work around missing information in ManagedArray
template<typename Container>
struct ContainerTraits {
  using value_type = typename Container::value_type;
};

template<typename Value>
struct ContainerTraits<chai::ManagedArray<Value>> {
  using value_type = Value;
};

// General case, returns by pointer
template<typename T, size_t stride>
struct ElementAccessor {
  using ContainerType   = typename T::ContainerType;
  using value_type      = typename ContainerTraits<ContainerType>::value_type*;
  using reference       = typename ContainerTraits<ContainerType>::value_type*;
  using const_reference = const typename ContainerTraits<ContainerType>::value_type*;
  SPHERAL_HOST_DEVICE static reference        at(      ContainerType& values, const size_t k) { REQUIRE(stride*k < values.size()); return &values[stride*k]; }
  SPHERAL_HOST_DEVICE static const_reference  at(const ContainerType& values, const size_t k) { REQUIRE(stride*k < values.size()); return &values[stride*k]; }
};

// When we only have one element per index, return it by reference
template<typename T>
struct ElementAccessor<T, 1u> {
  using ContainerType   = typename T::ContainerType;
  using value_type      = typename ContainerTraits<ContainerType>::value_type&;
  using reference       = typename ContainerTraits<ContainerType>::value_type&;
  using const_reference = const typename ContainerTraits<ContainerType>::value_type&;
  SPHERAL_HOST_DEVICE static reference        at(      ContainerType& values, const size_t k) { REQUIRE(k < values.size()); return values[k]; }
  SPHERAL_HOST_DEVICE static const_reference  at(const ContainerType& values, const size_t k) { REQUIRE(k < values.size()); return values[k]; }
};

}
}

#endif
