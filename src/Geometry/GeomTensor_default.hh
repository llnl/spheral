//---------------------------------Spheral++----------------------------------//
// GeomTensor -- Geometric Tensor Class.
//
// Created by JMO, Thu Apr 22 20:48:23 PDT 1999
// Modified by:
//   Aug 1, 99: JMO, removing MTL inheritance due to compilation problems,
//              and adding symmetric tensor type.
//   2004-08-19: JMO, converting to use boost::ublas library.
//   2004-08-23:  JMO, ublas is still too slow, so going to primitive C 
//                internal data types in accordance with suggestions from
//                Brian White
//   2004-11-10: JMO, experimenting with replacing the double array with individual
//               double elements, again for speed.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomTensor_default_hh__
#define __Spheral_GeomTensor_default_hh__

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"
#include "Geometry/GeomTensorBase.hh"

#include <iostream>
#include "Eigen/Dense"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace Spheral {

template<int nDim>
class GeomTensor: public GeomTensorBase<nDim> {

public:
  //--------------------------- Public Interface ---------------------------//
  using const_iterator = const double*;
  using iterator = double*;
  using size_type = unsigned;
  using EigenType = Eigen::Matrix<double, nDim, nDim>;

  // Useful stuff known at compile time
  SPHERAL_HOST_DEVICE static constexpr size_type nDimensions()            { return nDim; }
  SPHERAL_HOST_DEVICE static constexpr size_type numElements()            { return nDim * nDim; }
  SPHERAL_HOST_DEVICE static constexpr GeomTensor zero()                  { return GeomTensor(); }
  SPHERAL_HOST_DEVICE static constexpr GeomTensor one()                   { return GeomTensor(1.0); }

  // Constructors.
  SPHERAL_HOST_DEVICE GeomTensor() = default;
  SPHERAL_HOST_DEVICE explicit GeomTensor(const double a11);
  SPHERAL_HOST_DEVICE GeomTensor(const double a11, const double a12,
             const double a21, const double a22);
  SPHERAL_HOST_DEVICE GeomTensor(const double a11, const double a12, const double a13,
             const double a21, const double a22, const double a23,
             const double a31, const double a32, const double a33);
  SPHERAL_HOST_DEVICE explicit GeomTensor(const GeomSymmetricTensor<nDim>& ten);
  GeomTensor(const EigenType& ten);
  template<typename Derived> GeomTensor(const Eigen::MatrixBase<Derived>& ten);

  // Assignment.
  SPHERAL_HOST_DEVICE GeomTensor& operator=(const GeomSymmetricTensor<nDim>& rhs);
  template<typename Derived> GeomTensor& operator=(const Eigen::MatrixBase<Derived>& rhs);

  // Access the elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(const size_type row, const size_type column) const;
  SPHERAL_HOST_DEVICE double& operator()(const size_type row, const size_type column);

  // More C++ style indexing.
  SPHERAL_HOST_DEVICE double operator[](size_type index) const;
  SPHERAL_HOST_DEVICE double& operator[](size_type index);

  // Access the individual elements by (x,y,z) notation.
  SPHERAL_HOST_DEVICE double xx() const;
  SPHERAL_HOST_DEVICE double xy() const;
  SPHERAL_HOST_DEVICE double xz() const;
  SPHERAL_HOST_DEVICE double yx() const;
  SPHERAL_HOST_DEVICE double yy() const;
  SPHERAL_HOST_DEVICE double yz() const;
  SPHERAL_HOST_DEVICE double zx() const;
  SPHERAL_HOST_DEVICE double zy() const;
  SPHERAL_HOST_DEVICE double zz() const;
  SPHERAL_HOST_DEVICE void xx(double val);
  SPHERAL_HOST_DEVICE void xy(double val);
  SPHERAL_HOST_DEVICE void xz(double val);
  SPHERAL_HOST_DEVICE void yx(double val);
  SPHERAL_HOST_DEVICE void yy(double val);
  SPHERAL_HOST_DEVICE void yz(double val);
  SPHERAL_HOST_DEVICE void zx(double val);
  SPHERAL_HOST_DEVICE void zy(double val);
  SPHERAL_HOST_DEVICE void zz(double val);

  // Get/set rows and columns.
  SPHERAL_HOST_DEVICE GeomVector<nDim> getRow(size_type index) const;
  SPHERAL_HOST_DEVICE GeomVector<nDim> getColumn(size_type index) const;
  SPHERAL_HOST_DEVICE void setRow(size_type index, const GeomVector<nDim>& vec);
  SPHERAL_HOST_DEVICE void setColumn(size_type index, const GeomVector<nDim>& vec);

  // Iterator access to the raw data.
  SPHERAL_HOST_DEVICE iterator begin();
  SPHERAL_HOST_DEVICE iterator end();

  SPHERAL_HOST_DEVICE const_iterator begin() const;
  SPHERAL_HOST_DEVICE const_iterator end() const;

  // The zero and identity matrices.
  SPHERAL_HOST_DEVICE void Zero();
  SPHERAL_HOST_DEVICE void Identity();

  SPHERAL_HOST_DEVICE GeomTensor operator-() const;

  SPHERAL_HOST_DEVICE GeomTensor operator+(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor operator-(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor operator*(const GeomTensor& rhs) const;

  SPHERAL_HOST_DEVICE GeomTensor operator+(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor operator-(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor operator*(const GeomSymmetricTensor<nDim>& rhs) const;

  SPHERAL_HOST_DEVICE GeomVector<nDim> operator*(const GeomVector<nDim>& rhs) const;
  // GeomTensor operator+(const double val) const;
  // GeomTensor operator-(const double val) const;
  SPHERAL_HOST_DEVICE GeomTensor operator*(const double val) const;
  SPHERAL_HOST_DEVICE GeomTensor operator/(const double val) const;

  SPHERAL_HOST_DEVICE GeomTensor& operator+=(const GeomTensor& rhs);
  SPHERAL_HOST_DEVICE GeomTensor& operator-=(const GeomTensor& rhs);
  SPHERAL_HOST_DEVICE GeomTensor& operator*=(const GeomTensor& rhs);

  SPHERAL_HOST_DEVICE GeomTensor& operator+=(const GeomSymmetricTensor<nDim>& rhs);
  SPHERAL_HOST_DEVICE GeomTensor& operator-=(const GeomSymmetricTensor<nDim>& rhs);
  SPHERAL_HOST_DEVICE GeomTensor& operator*=(const GeomSymmetricTensor<nDim>& rhs);

  template<typename Derived> GeomTensor& operator+=(const Eigen::MatrixBase<Derived>& rhs);
  template<typename Derived> GeomTensor& operator-=(const Eigen::MatrixBase<Derived>& rhs);

  SPHERAL_HOST_DEVICE GeomTensor& operator*=(const double rhs);
  SPHERAL_HOST_DEVICE GeomTensor& operator/=(const double rhs);

  SPHERAL_HOST_DEVICE bool operator==(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const GeomTensor& rhs) const;

  SPHERAL_HOST_DEVICE bool operator==(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const GeomSymmetricTensor<nDim>& rhs) const;

  // bool operator==(const double rhs) const;
  // bool operator!=(const double rhs) const;
  // bool operator<(const double rhs) const;
  // bool operator>(const double rhs) const;
  // bool operator<=(const double rhs) const;
  // bool operator>=(const double rhs) const;

  SPHERAL_HOST_DEVICE GeomSymmetricTensor<nDim> Symmetric() const;
  SPHERAL_HOST_DEVICE GeomTensor SkewSymmetric() const;
  SPHERAL_HOST_DEVICE GeomTensor Transpose() const;
  SPHERAL_HOST_DEVICE GeomTensor Inverse() const;
  SPHERAL_HOST_DEVICE GeomVector<nDim> diagonalElements() const;
  SPHERAL_HOST_DEVICE double Trace() const;
  SPHERAL_HOST_DEVICE double Determinant() const;
  SPHERAL_HOST_DEVICE GeomVector<nDim> dot(const GeomVector<nDim>& vec) const;
  SPHERAL_HOST_DEVICE GeomTensor dot(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor dot(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot(const GeomTensor& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot(const GeomSymmetricTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE double selfDoubledot() const;

  // Return the square of this tensor (using matrix multiplication).
  SPHERAL_HOST_DEVICE GeomTensor square() const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  SPHERAL_HOST_DEVICE GeomTensor squareElements() const;

  // A simple method for returning the eigenvalues of a tensor.
  GeomVector<nDim> eigenValues() const;
  
  // Apply a rotational transform to this tensor ( R^-1 * (*this) * R ).
  SPHERAL_HOST_DEVICE void rotationalTransform(const GeomTensor& R);

  // Return the max absolute element.
  SPHERAL_HOST_DEVICE double maxAbsElement() const;

  //  Convert to an Eigen Vector
  EigenType eigen() const;

private:
  //--------------------------- Private Interface ---------------------------//
  SPHERAL_HOST_DEVICE size_type elementIndex(const size_type row, const size_type column) const;
};

// Declare specializations.
template<> SPHERAL_HOST_DEVICE GeomTensor<2>::GeomTensor(const double, const double,
                                                         const double, const double);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>::GeomTensor(const double, const double, const double,
                                                         const double, const double, const double,
                                                         const double, const double, const double);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator=(const GeomSymmetricTensor<1>& rhs);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator=(const GeomSymmetricTensor<2>& rhs);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator=(const GeomSymmetricTensor<3>& rhs);

template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::xy() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::xz() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::yx() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::yy() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::yz() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::zx() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::zy() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::zz() const;

template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::xz() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::yz() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::zx() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::zy() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::zz() const;

template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::xy(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::xz(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::yx(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::yy(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::yz(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::zx(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::zy(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::zz(const double);

template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::xz(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::yz(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::zx(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::zy(const double);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::zz(const double);

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomTensor<1>::getRow(const GeomTensor<1>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomTensor<2>::getRow(const GeomTensor<2>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomTensor<3>::getRow(const GeomTensor<3>::size_type) const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomTensor<1>::getColumn(const GeomTensor<1>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomTensor<2>::getColumn(const GeomTensor<2>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomTensor<3>::getColumn(const GeomTensor<3>::size_type) const;

template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::setRow(const GeomTensor<1>::size_type, const GeomVector<1>&);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::setRow(const GeomTensor<2>::size_type, const GeomVector<2>&);
template<> SPHERAL_HOST_DEVICE void GeomTensor<3>::setRow(const GeomTensor<3>::size_type, const GeomVector<3>&);

template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::setColumn(const GeomTensor<1>::size_type, const GeomVector<1>&);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::setColumn(const GeomTensor<2>::size_type, const GeomVector<2>&);
template<> SPHERAL_HOST_DEVICE void GeomTensor<3>::setColumn(const GeomTensor<3>::size_type, const GeomVector<3>&);

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::operator-() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::operator-() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::operator-() const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::operator*(const double) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::operator*(const double) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::operator*(const double) const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::operator/(const double) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::operator/(const double) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::operator/(const double) const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator+=(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator+=(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator+=(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator+=(const GeomSymmetricTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator+=(const GeomSymmetricTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator+=(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator-=(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator-=(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator-=(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator-=(const GeomSymmetricTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator-=(const GeomSymmetricTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator-=(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator*=(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator*=(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator*=(const GeomTensor<3>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator*=(const GeomSymmetricTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator*=(const GeomSymmetricTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator*=(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator*=(const double);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator*=(const double);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator*=(const double);

template<> SPHERAL_HOST_DEVICE GeomTensor<1>& GeomTensor<1>::operator/=(const double);
template<> SPHERAL_HOST_DEVICE GeomTensor<2>& GeomTensor<2>::operator/=(const double);
template<> SPHERAL_HOST_DEVICE GeomTensor<3>& GeomTensor<3>::operator/=(const double);

template<> SPHERAL_HOST_DEVICE bool GeomTensor<1>::operator==(const GeomTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomTensor<2>::operator==(const GeomTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomTensor<3>::operator==(const GeomTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE bool GeomTensor<1>::operator==(const GeomSymmetricTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomTensor<2>::operator==(const GeomSymmetricTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomTensor<3>::operator==(const GeomSymmetricTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomTensor<1>::Symmetric() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomTensor<2>::Symmetric() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomTensor<3>::Symmetric() const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::SkewSymmetric() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::SkewSymmetric() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::SkewSymmetric() const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::Transpose() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::Transpose() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::Transpose() const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::Inverse() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::Inverse() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::Inverse() const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomTensor<1>::diagonalElements() const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomTensor<2>::diagonalElements() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomTensor<3>::diagonalElements() const;

template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::Trace() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::Trace() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<3>::Trace() const;

template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::Determinant() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::Determinant() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<3>::Determinant() const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomTensor<1>::dot(const GeomVector<1>&) const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomTensor<2>::dot(const GeomVector<2>&) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomTensor<3>::dot(const GeomVector<3>&) const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::dot(const GeomTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::dot(const GeomTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::dot(const GeomTensor<3>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::dot(const GeomSymmetricTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::dot(const GeomSymmetricTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::dot(const GeomSymmetricTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::doubledot(const GeomTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::doubledot(const GeomTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<3>::doubledot(const GeomTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::doubledot(const GeomSymmetricTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::doubledot(const GeomSymmetricTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<3>::doubledot(const GeomSymmetricTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::square() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::square() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::square() const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomTensor<1>::squareElements() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomTensor<2>::squareElements() const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomTensor<3>::squareElements() const;

template<> SPHERAL_HOST_DEVICE void GeomTensor<1>::rotationalTransform(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE void GeomTensor<2>::rotationalTransform(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE void GeomTensor<3>::rotationalTransform(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE double GeomTensor<1>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<2>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomTensor<3>::maxAbsElement() const;

// Forward declare the global functions.
template<int nDim> GeomTensor<nDim> operator*(double lhs, const GeomTensor<nDim>& rhs);
template<int nDim> ::std::istream& operator>>(std::istream& is, GeomTensor<nDim>& ten);
template<int nDim> std::ostream& operator<<(std::ostream& os, const GeomTensor<nDim>& ten);

#if defined(_OPENMP) && _OPENMP >= 201107
#pragma omp declare reduction(tensadd : GeomTensor<1> : omp_out += omp_in ) initializer( omp_priv = GeomTensor<1>(0.0) )
#pragma omp declare reduction(tensdif : GeomTensor<1> : omp_out -= omp_in ) initializer( omp_priv = GeomTensor<1>(0.0) )
#pragma omp declare reduction(tensadd : GeomTensor<2> : omp_out += omp_in ) initializer( omp_priv = GeomTensor<2>(0.0,0.0,0.0,0.0)) 
#pragma omp declare reduction(tensdif : GeomTensor<2> : omp_out -= omp_in ) initializer( omp_priv = GeomTensor<2>(0.0,0.0,0.0,0.0))
#pragma omp declare reduction(tensadd : GeomTensor<3> : omp_out += omp_in ) initializer( omp_priv = GeomTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(tensdif : GeomTensor<3> : omp_out -= omp_in ) initializer( omp_priv = GeomTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) ) 
#endif

}

#include "GeomTensorInline_default.hh"

#endif
