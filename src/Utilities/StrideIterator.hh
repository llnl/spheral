//---------------------------------Spheral++----------------------------------//
// StrideIterator
//
// Provide a basic STL compliant iterator that allows us to specify a stride
//
// Created by J. Michael Owen, Fri Dec 20 10:57:27 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrideIterator__
#define __Spheral_StrideIterator__

#include <cstddef>
#include <iterator>
#include <type_traits>

namespace Spheral {

template<typename T, size_t stride>
class StrideIterator {
public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = std::remove_cv_t<T>;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using reference = T&;
  SPHERAL_HOST_DEVICE StrideIterator(T* ptr): mptr(ptr)                  {}

  SPHERAL_HOST_DEVICE T& operator*()                               const { return *mptr; }
  SPHERAL_HOST_DEVICE StrideIterator& operator++()                       { mptr += stride; return *this; }
  SPHERAL_HOST_DEVICE StrideIterator operator++(int)                     { StrideIterator tmp = *this; ++(*this); return tmp; }
  SPHERAL_HOST_DEVICE StrideIterator& operator--()                       { mptr -= stride; return *this; }
  SPHERAL_HOST_DEVICE StrideIterator operator--(int)                     { StrideIterator tmp = *this; --(*this); return tmp; }

  SPHERAL_HOST_DEVICE bool operator==(const StrideIterator& other) const { return mptr == other.mptr; }
  SPHERAL_HOST_DEVICE bool operator!=(const StrideIterator& other) const { return !(*this == other); }

private:
  template<typename, size_t> friend class StrideIterator;
  pointer mptr;
};

}

#endif
