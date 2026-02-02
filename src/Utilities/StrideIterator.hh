//---------------------------------Spheral++----------------------------------//
// StrideIterator
//
// Provide a basic STL compliant iterator that allows us to specify a stride
//
// Created by J. Michael Owen, Fri Dec 20 10:57:27 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrideIterator__
#define __Spheral_StrideIterator__

#include <iterator>

namespace Spheral {

template<typename T, size_t stride>
class StrideIterator: public std::iterator<std::random_access_iterator_tag, T> {
public:
  StrideIterator(T* ptr): mptr(ptr)                  {}

  T& operator*()                               const { return *mptr; }
  StrideIterator& operator++()                       { mptr += stride; return *this; }
  StrideIterator operator++(int)                     { StrideIterator tmp = *this; ++(*this); return tmp; }
  StrideIterator& operator--()                       { mptr -= stride; return *this; }
  StrideIterator operator--(int)                     { StrideIterator tmp = *this; --(*this); return tmp; }

  bool operator==(const StrideIterator& other) const { return mptr == other.mptr; }
  bool operator!=(const StrideIterator& other) const { return !(*this == other); }

private:
  T* mptr;
};

}

#endif
