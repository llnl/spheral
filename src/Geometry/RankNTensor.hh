//---------------------------------Spheral++----------------------------------//
// RankNTensor -- Arbitrary rank tensor.  Assumes each nth index spans 
//                Dimension::nDim.
//
// Created by JMO, Sun Oct 11 10:38:26 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_RankNTensor_hh__
#define __Spheral_RankNTensor_hh__

#include "Utilities/FastMath.hh"

#include <iostream>
#include <array>

namespace Spheral {

template<int nDim, int rank, typename Descendant>
class RankNTensor {

public:
  //--------------------------- Public Interface ---------------------------//
  using const_iterator = const double*;
  using iterator = double*;
  using size_type = unsigned;

  // Useful static member data.
  SPHERAL_HOST_DEVICE static constexpr size_type nrank()                      { return rank; }
  SPHERAL_HOST_DEVICE static constexpr size_type nDimensions()                { return nDim; }
  SPHERAL_HOST_DEVICE static constexpr size_type numElements()                { return FastMath::calcPower(nDim, rank); }

  // Constructors.
  SPHERAL_HOST_DEVICE RankNTensor()                                           = default;
  SPHERAL_HOST_DEVICE explicit RankNTensor(const double val): mElements()     { mElements.fill(val); }
  SPHERAL_HOST_DEVICE RankNTensor(const RankNTensor& rhs)                     = default;

  // Assignment.
  SPHERAL_HOST_DEVICE RankNTensor& operator=(const RankNTensor& rhs)          = default;
  SPHERAL_HOST_DEVICE RankNTensor& operator=(const double rhs)                { for (auto i = 0u; i < numElements(); ++i) mElements[i] = rhs; }

  // More C++ style indexing.
  SPHERAL_HOST_DEVICE double  operator[](size_type index) const               { REQUIRE(index < numElements()); return mElements[index]; }
  SPHERAL_HOST_DEVICE double& operator[](size_type index)                     { REQUIRE(index < numElements()); return mElements[index]; }

  // Iterator access to the raw data.
  SPHERAL_HOST_DEVICE iterator begin()                                        { return mElements.begin(); }
  SPHERAL_HOST_DEVICE iterator end()                                          { return mElements.end(); }

  SPHERAL_HOST_DEVICE const_iterator begin() const                            { return mElements.begin(); }
  SPHERAL_HOST_DEVICE const_iterator end() const                              { return mElements.end(); }

  // Zero out the tensor.
  SPHERAL_HOST_DEVICE constexpr void Zero()                                   { for (auto i = 0u; i < numElements(); ++i) mElements[i] = 0.0; }

  // Assorted operations.
  SPHERAL_HOST_DEVICE Descendant operator-() const;

  SPHERAL_HOST_DEVICE Descendant& operator+=(const RankNTensor& rhs);
  SPHERAL_HOST_DEVICE Descendant& operator-=(const RankNTensor& rhs);

  SPHERAL_HOST_DEVICE Descendant operator+(const RankNTensor& rhs) const;
  SPHERAL_HOST_DEVICE Descendant operator-(const RankNTensor& rhs) const;

  SPHERAL_HOST_DEVICE Descendant& operator*=(const double rhs);
  SPHERAL_HOST_DEVICE Descendant& operator/=(const double rhs);

  SPHERAL_HOST_DEVICE Descendant operator*(const double rhs) const;
  SPHERAL_HOST_DEVICE Descendant operator/(const double rhs) const;

  SPHERAL_HOST_DEVICE bool operator==(const RankNTensor& rhs) const          { return mElements == rhs.mElements; }
  SPHERAL_HOST_DEVICE bool operator!=(const RankNTensor& rhs) const          { return !(*this == rhs); }

  SPHERAL_HOST_DEVICE bool operator==(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const double rhs) const                { return !(*this == rhs); }

  SPHERAL_HOST_DEVICE double doubledot(const RankNTensor<nDim, rank, Descendant>& rhs) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  SPHERAL_HOST_DEVICE Descendant squareElements() const;

  // Return the max absolute element.
  SPHERAL_HOST_DEVICE double maxAbsElement() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  std::array<double, numElements()>  mElements = {};
};

// Forward declare the global functions.
template<int nDim, int rank, typename Descendant> Descendant operator*(const double lhs, const RankNTensor<nDim, rank, Descendant>& rhs);

template<int nDim, int rank, typename Descendant> ::std::istream& operator>>(std::istream& is, RankNTensor<nDim, rank, Descendant>& ten);
template<int nDim, int rank, typename Descendant> ::std::ostream& operator<<(std::ostream& os, const RankNTensor<nDim, rank, Descendant>& ten);

}

#include "RankNTensorInline.hh"

#endif
