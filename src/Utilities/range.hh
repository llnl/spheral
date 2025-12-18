//---------------------------------Spheral++----------------------------------//
// range_iterations
//
// A few useful adapters to make iterating over ranges defined by pairs of
// iterators nicer.
//----------------------------------------------------------------------------//
#ifndef __Spheral_range_iterations__
#define __Spheral_range_iterations__

#include <boost/range/iterator_range.hpp>
#include <utility>
#include <tuple>

namespace Spheral {

//------------------------------------------------------------------------------
// Thin wrapper around boost::make_iterator_range to construct range-based loops
// from iterators.  Only real reason for this wrapper is so we can easily try
// something else out eventually without changing all the loops in Spheral.
//------------------------------------------------------------------------------
template<typename Iter>
constexpr auto range(Iter&& begin, Iter&& end) {
  return boost::make_iterator_range(begin, end);
}

//------------------------------------------------------------------------------
// Nice little C++ method to emulate the Python "enumerate" function, which
// returns an index as well as the values for each value in a container.
// Original version found at
// https://www.reedbeta.com/blog/python-like-enumerate-in-cpp17/
//------------------------------------------------------------------------------
template <typename T,
          typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable) {
  struct iterator {
    size_t i;
    TIter iter;
    bool operator != (const iterator & other) const { return iter != other.iter; }
    void operator ++ () { ++i; ++iter; }
    auto operator * () const { return std::tie(i, *iter); }
  };
  struct iterable_wrapper {
    T iterable;
    auto begin() { return iterator{ 0, std::begin(iterable) }; }
    auto end() { return iterator{ 0, std::end(iterable) }; }
  };
  return iterable_wrapper{ std::forward<T>(iterable) };
}

//------------------------------------------------------------------------------
// Iterator based version of enumerate
//------------------------------------------------------------------------------
template<typename Iter>
constexpr auto enumerate(Iter&& begin, Iter&& end) {
  return enumerate(range(begin, end));
}

}

#endif
