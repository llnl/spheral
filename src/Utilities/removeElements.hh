//---------------------------------Spheral++----------------------------------//
// removeElements
//
// Removes the specified elements by index from a std::vector.
// This is needed because there doesn't seem to be a good solution for this in
// the STL.  The std::vector::erase method involves N^2 operations when you're
// removing many elements, while the std::remove and std::remove_if do not work
// for removing elements by index/iterator.
//
// Created by JMO, Mon Dec 14 13:30:19 PST 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_removeElements__
#define __Spheral_removeElements__

#include <vector>
#include <algorithm>
#include "DBC.hh"

#ifdef USE_UVM
#include "../Field/uvm_allocator.hh"
template<typename DataType>
using DataAllocator = typename uvm_allocator::UVMAllocator<DataType>;
#else
template<typename DataType>
using DataAllocator = typename std::allocator<DataType>;
#endif

namespace Spheral {

template<typename Value, typename Allocator, typename index_t>
inline
void
removeElements(std::vector<Value, Allocator>& vec,
               const std::vector<index_t>& elements) {

  // Is there anything to do?
  if (not elements.empty()) {
    const index_t originalSize = vec.size();
    const index_t newSize = originalSize - elements.size();

    // Pre-conditions.
    BEGIN_CONTRACT_SCOPE
    {
      // We require the input IDs be sorted and unique.
      for (typename std::vector<index_t>::const_iterator itr = elements.begin();
           itr + 1 < elements.end();
           ++itr) {
        REQUIRE(*itr < *(itr + 1));
      }
      if (elements.size() > 0) {
        REQUIRE(elements[0] >= 0 && elements.back() < originalSize);
      }
    }
    END_CONTRACT_SCOPE

    // A single value is trivial.
    if (elements.size() == 1) {
      vec.erase(vec.begin() + elements[0]);

    } else {

      // Remove the elements.
      // We prefer not to use the vector::erase here 'cause if we're removing
      // many elements the copy and move behaviour of erase can make this
      // an N^2 thing.  Yuck!
      auto i = elements[0];
      for (auto k = 1; k < elements.size(); ++k) {
        std::copy(vec.begin() + i + 1, vec.begin() + elements[k], vec.begin() + i);
        i = elements[k];
      }

      // Resize vec to it's new size.
      vec.erase(vec.begin() + newSize, vec.end());

    }

    // Post-conditions.
    ENSURE(vec.size() == newSize);
  }
}


template<typename Value, typename index_t>
inline
void
removeElements(std::vector<Value,DataAllocator<Value>>& vec,
               const std::vector<index_t,DataAllocator<index_t>>& elements) {

  // Is there anything to do?
  if (elements.size() > 0) {

    const index_t originalSize = vec.size();
    const index_t newSize = originalSize - elements.size();

    // Pre-conditions.
    BEGIN_CONTRACT_SCOPE
    {
      // We require the input IDs be sorted and unique.
      for (typename std::vector<index_t,DataAllocator<index_t>>::const_iterator itr = elements.begin();
           itr + 1 < elements.end();
           ++itr) {
        REQUIRE(*itr < *(itr + 1));
      }
      if (elements.size() > 0) {
        REQUIRE(elements.back() < originalSize);
      }
    }
    END_CONTRACT_SCOPE
    // Remove the elements.
    // We prefer not to use the vector::erase here 'cause if we're removing
    // many elements the copy and move behaviour of erase can make this
    // an N^2 thing.  Yuck!
    typename std::vector<index_t,DataAllocator<index_t>>::const_iterator delItr = elements.begin();
    typename std::vector<index_t,DataAllocator<index_t>>::const_iterator endItr = elements.end();
    index_t i = *delItr;
    index_t j = i + 1;
    ++delItr;
    while (j != originalSize and delItr != endItr) {
      if (j == *delItr) {
        ++delItr;
        ++j;
      } else {
        vec[i] = vec[j];
        ++i;
        ++j;
      }
    }
    if (j != originalSize) copy(vec.begin() + j, vec.end(), vec.begin() + i);

    // Resize vec to it's new size.
    vec.erase(vec.begin() + newSize, vec.end());

    // Post-conditions.
    ENSURE((int)vec.size() == (int)newSize);
  }
}


} // end Spheral namespace

#endif
