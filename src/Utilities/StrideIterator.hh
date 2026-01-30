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

  constexpr StrideIterator() noexcept: mptr(nullptr) {}
  constexpr explicit StrideIterator(pointer ptr) noexcept: mptr(ptr) {}

  template<typename U, typename = std::enable_if_t<std::is_convertible_v<U*, T*>>>
  constexpr StrideIterator(const StrideIterator<U, stride>& rhs) noexcept:
    mptr(rhs.mptr) {}

  constexpr reference operator*() const noexcept            { return *mptr; }
  constexpr pointer operator->() const noexcept             { return mptr; }

  constexpr StrideIterator& operator++() noexcept           { mptr += stride; return *this; }
  constexpr StrideIterator operator++(int) noexcept         { auto tmp = *this; ++(*this); return tmp; }
  constexpr StrideIterator& operator--() noexcept           { mptr -= stride; return *this; }
  constexpr StrideIterator operator--(int) noexcept         { auto tmp = *this; --(*this); return tmp; }

  constexpr StrideIterator& operator+=(difference_type n) noexcept { mptr += n*stride; return *this; }
  constexpr StrideIterator& operator-=(difference_type n) noexcept { mptr -= n*stride; return *this; }
  constexpr StrideIterator operator+(difference_type n) const noexcept { auto tmp = *this; return tmp += n; }
  constexpr StrideIterator operator-(difference_type n) const noexcept { auto tmp = *this; return tmp -= n; }

  constexpr reference operator[](difference_type n) const noexcept  { return *(*this + n); }

  constexpr bool operator==(const StrideIterator& other) const noexcept { return mptr == other.mptr; }
  constexpr bool operator!=(const StrideIterator& other) const noexcept { return mptr != other.mptr; }
  constexpr bool operator<(const StrideIterator& other) const noexcept  { return mptr < other.mptr; }
  constexpr bool operator>(const StrideIterator& other) const noexcept  { return mptr > other.mptr; }
  constexpr bool operator<=(const StrideIterator& other) const noexcept { return mptr <= other.mptr; }
  constexpr bool operator>=(const StrideIterator& other) const noexcept { return mptr >= other.mptr; }

  constexpr difference_type operator-(const StrideIterator& other) const noexcept {
    return (mptr - other.mptr)/static_cast<difference_type>(stride);
  }

  friend constexpr StrideIterator operator+(difference_type n, StrideIterator it) noexcept { return it += n; }

private:
  template<typename, size_t> friend class StrideIterator;
  pointer mptr;
};

}

#endif
