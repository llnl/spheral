//---------------------------------Spheral++----------------------------------//
// RankNTensor -- Arbitrary rank tensor.  Assumes each nth index spans 
//                Dimension::nDim.
//
// Created by JMO, Sun Oct 11 10:38:26 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_RankNTensor_hh__
#define __Spheral_RankNTensor_hh__

#include <iostream>
#include "Utilities/FastMath.hh"

namespace Spheral {

template<int nDim, int rank, typename Descendant>
class RankNTensor {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const size_type nrank;
  static const size_type nDimensions;
  static constexpr size_type numElements = FastMath::calcPower(nDim, rank);

  // Constructors.
  SPHERAL_HOST_DEVICE RankNTensor() = default;
  SPHERAL_HOST_DEVICE explicit RankNTensor(const double val);
  SPHERAL_HOST_DEVICE RankNTensor(const RankNTensor& rhs) = default;

  // Assignment.
  SPHERAL_HOST_DEVICE RankNTensor& operator=(const RankNTensor& rhs) = default;
  SPHERAL_HOST_DEVICE RankNTensor& operator=(const double rhs);

  // More C++ style indexing.
  SPHERAL_HOST_DEVICE double operator[](size_type index) const;
  SPHERAL_HOST_DEVICE double& operator[](size_type index);

  // Iterator access to the raw data.
  SPHERAL_HOST_DEVICE iterator begin();
  SPHERAL_HOST_DEVICE iterator end();

  SPHERAL_HOST_DEVICE const_iterator begin() const;
  SPHERAL_HOST_DEVICE const_iterator end() const;

  // Zero out the tensor.
  SPHERAL_HOST_DEVICE void Zero();

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

  SPHERAL_HOST_DEVICE bool operator==(const RankNTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const RankNTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const RankNTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const RankNTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const RankNTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const RankNTensor& rhs) const;

  SPHERAL_HOST_DEVICE bool operator==(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const double rhs) const;

  SPHERAL_HOST_DEVICE double doubledot(const RankNTensor<nDim, rank, Descendant>& rhs) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  SPHERAL_HOST_DEVICE Descendant squareElements() const;

  // Return the max absolute element.
  SPHERAL_HOST_DEVICE double maxAbsElement() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  double mElements[numElements] = {};
};

// Forward declare the global functions.
template<int nDim, int rank, typename Descendant> Descendant operator*(const double lhs, const RankNTensor<nDim, rank, Descendant>& rhs);

template<int nDim, int rank, typename Descendant> ::std::istream& operator>>(std::istream& is, RankNTensor<nDim, rank, Descendant>& ten);
template<int nDim, int rank, typename Descendant> ::std::ostream& operator<<(std::ostream& os, const RankNTensor<nDim, rank, Descendant>& ten);

// Initialize our static variables.
template<int nDim, int rank, typename Descendant> const typename RankNTensor<nDim, rank, Descendant>::size_type RankNTensor<nDim, rank, Descendant>::nrank = rank;
template<int nDim, int rank, typename Descendant> const typename RankNTensor<nDim, rank, Descendant>::size_type RankNTensor<nDim, rank, Descendant>::nDimensions = nDim;

}

#include "RankNTensorInline.hh"

#endif
