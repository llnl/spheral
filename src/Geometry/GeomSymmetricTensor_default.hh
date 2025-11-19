//---------------------------------Spheral++----------------------------------//
// GeomSymmetricTensor
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
#ifndef __Spheral_GeomSymmetricTensor_hh_
#define __Spheral_GeomSymmetricTensor_hh_

#include "config.hh"

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"
#include "Geometry/EigenStruct_fwd.hh"
#include "Geometry/GeomSymmetricTensorBase.hh"

#include <iostream>
#include "Eigen/Dense"

namespace Spheral {

template<int nDim>
class GeomSymmetricTensor: public GeomSymmetricTensorBase<nDim> {

public:
  //--------------------------- Public Interface ---------------------------//
  using const_iterator = const double*;
  using iterator = double*;
  using size_type = unsigned;
  using EigenType = Eigen::Matrix<double, nDim, nDim>;
  using EigenStructType = EigenStruct<nDim>;

  // Useful stuff known at compile time
  SPHERAL_HOST_DEVICE static constexpr size_type nDimensions()            { return nDim; }
  SPHERAL_HOST_DEVICE static constexpr size_type numElements()            { return (nDim * (nDim+1)) / 2; }
  SPHERAL_HOST_DEVICE static constexpr GeomSymmetricTensor zero()         { return GeomSymmetricTensor(); }
  SPHERAL_HOST_DEVICE static constexpr GeomSymmetricTensor one();
  SPHERAL_HOST_DEVICE static constexpr double onethird()                  { return 1.0/3.0; }
  SPHERAL_HOST_DEVICE static constexpr double sqrt3()                     { return std::sqrt(3.0); }

  // Constructors.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor() = default;
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const double a11);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const double /*a11*/, const double /*a12*/,
                                          const double /*a21*/, const double /*a22*/);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const double /*a11*/, const double /*a12*/, const double /*a13*/,
                                          const double /*a21*/, const double /*a22*/, const double /*a23*/,
                                          const double /*a31*/, const double /*a32*/, const double /*a33*/);
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const GeomTensor<nDim>& ten);
  template<typename Derived> GeomSymmetricTensor(const Eigen::MatrixBase<Derived>& ten);

  // Assignment.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator=(const GeomTensor<nDim>& rhs);
  template<typename Derived> GeomSymmetricTensor& operator=(const Eigen::MatrixBase<Derived>& rhs);

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

  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator-() const;

  SPHERAL_HOST_DEVICE GeomTensor<nDim> operator+(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> operator-(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> operator*(const GeomTensor<nDim>& rhs) const;

  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator+(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator-(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim>    operator*(const GeomSymmetricTensor& rhs) const;

  SPHERAL_HOST_DEVICE GeomVector<nDim> operator*(const GeomVector<nDim>& rhs) const;
  // GeomSymmetricTensor operator+(const double rhs) const;
  // GeomSymmetricTensor operator-(const double rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator*(const double rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator/(const double rhs) const;

  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator+=(const GeomSymmetricTensor& rhs);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator-=(const GeomSymmetricTensor& rhs);

  template<typename Derived> GeomSymmetricTensor& operator+=(const Eigen::MatrixBase<Derived>& rhs);
  template<typename Derived> GeomSymmetricTensor& operator-=(const Eigen::MatrixBase<Derived>& rhs);

  // GeomSymmetricTensor& operator+=(const double rhs);
  // GeomSymmetricTensor& operator-=(const double rhs);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator*=(const double rhs);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator/=(const double rhs);

  SPHERAL_HOST_DEVICE bool operator==(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const GeomTensor<nDim>& rhs) const;

  SPHERAL_HOST_DEVICE bool operator==(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const GeomSymmetricTensor& rhs) const;

  SPHERAL_HOST_DEVICE bool operator==(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const double rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const double rhs) const;

  SPHERAL_HOST_DEVICE GeomSymmetricTensor Symmetric() const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> SkewSymmetric() const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor Transpose() const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor Inverse() const;
  SPHERAL_HOST_DEVICE GeomVector<nDim> diagonalElements() const;
  SPHERAL_HOST_DEVICE double Trace() const;
  SPHERAL_HOST_DEVICE double Determinant() const;
  SPHERAL_HOST_DEVICE GeomVector<nDim> dot(const GeomVector<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> dot(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> dot(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot(const GeomTensor<nDim>& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE double selfDoubledot() const;

  // Return the square of this tensor (using matrix multiplication).  Note that
  // for a symmetric tensor this is guaranteed to return a symmetric product.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor square() const;

  // Same idea for the cube.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor cube() const;

  // Compute the "square root" of the tensor: the tensor that, 
  // when squared, equals this tensor.
  GeomSymmetricTensor sqrt() const;

  // Similarly, compute the cube root.
  GeomSymmetricTensor cuberoot() const;

  // The general version, raise to an arbitrary power.
  GeomSymmetricTensor pow(const double p) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor squareElements() const;

  // A simple method for returning the eigenvalues of a tensor.
  SPHERAL_HOST_DEVICE GeomVector<nDim> eigenValues() const;
  
  // Apply a rotational transform to this tensor ( R^-1 * (*this) * R ).
  SPHERAL_HOST_DEVICE void rotationalTransform(const GeomTensor<nDim>& R);

  // We also provide a method to retrieve the eigenvectors as an EigenStruct.
  // Note that the eigen vectors are the columns of the full tensor in the resulting
  // struct.
  EigenStructType eigenVectors() const;

  // Return the max absolute element.
  SPHERAL_HOST_DEVICE double maxAbsElement() const;

  // Enforce min/max eigenvalues
  SPHERAL_HOST_DEVICE GeomSymmetricTensor enforceMinEigenValue(const double& x) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor enforceMaxEigenValue(const double& x) const;

  //  Convert to an Eigen Vector
  EigenType eigen() const;

private:
  //--------------------------- Private Interface ---------------------------//
  SPHERAL_HOST_DEVICE size_type elementIndex(const size_type row, const size_type column) const;
};

template<> SPHERAL_HOST_DEVICE constexpr GeomSymmetricTensor<1> GeomSymmetricTensor<1>::one() { GeomSymmetricTensor result; result.mxx = 1.0; return result; }
template<> SPHERAL_HOST_DEVICE constexpr GeomSymmetricTensor<2> GeomSymmetricTensor<2>::one() { GeomSymmetricTensor result; result.mxx = 1.0; result.myy = 1.0; return result; }
template<> SPHERAL_HOST_DEVICE constexpr GeomSymmetricTensor<3> GeomSymmetricTensor<3>::one() { GeomSymmetricTensor result; result.mxx = 1.0; result.myy = 1.0; result.mzz = 1.0; return result; }

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomSymmetricTensor<1>::eigenValues() const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomSymmetricTensor<2>::eigenValues() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::eigenValues() const;

template<> EigenStruct<1> GeomSymmetricTensor<1>::eigenVectors() const;
template<> EigenStruct<2> GeomSymmetricTensor<2>::eigenVectors() const;
template<> EigenStruct<3> GeomSymmetricTensor<3>::eigenVectors() const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>::GeomSymmetricTensor(const double, const double,
                                                                           const double, const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>::GeomSymmetricTensor(const double, const double, const double,
                                                                           const double, const double, const double,
                                                                           const double, const double, const double);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>::GeomSymmetricTensor(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>::GeomSymmetricTensor(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>::GeomSymmetricTensor(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator=(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator=(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator=(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::xy() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::xz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::yx() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::yy() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::yz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::zx() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::zy() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::zz() const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::xz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::yz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::zx() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::zy() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::zz() const;

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::xy(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::xz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::yx(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::yy(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::yz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::zx(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::zy(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::zz(const double);

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::xz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::yz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::zx(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::zy(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::zz(const double);

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomSymmetricTensor<1>::getRow(const GeomSymmetricTensor<1>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomSymmetricTensor<2>::getRow(const GeomSymmetricTensor<2>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::getRow(const GeomSymmetricTensor<3>::size_type) const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomSymmetricTensor<1>::getColumn(const GeomSymmetricTensor<1>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomSymmetricTensor<2>::getColumn(const GeomSymmetricTensor<2>::size_type) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::getColumn(const GeomSymmetricTensor<3>::size_type) const;

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::setRow(const GeomSymmetricTensor<1>::size_type, const GeomVector<1>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::setRow(const GeomSymmetricTensor<2>::size_type, const GeomVector<2>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<3>::setRow(const GeomSymmetricTensor<3>::size_type, const GeomVector<3>&);

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::setColumn(const GeomSymmetricTensor<1>::size_type, const GeomVector<1>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::setColumn(const GeomSymmetricTensor<2>::size_type, const GeomVector<2>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<3>::setColumn(const GeomSymmetricTensor<3>::size_type, const GeomVector<3>&);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::operator-() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::operator-() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::operator-() const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::operator*(const double) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::operator*(const double) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::operator*(const double) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::operator/(const double) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::operator/(const double) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::operator/(const double) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator+=(const GeomSymmetricTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator+=(const GeomSymmetricTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator+=(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator-=(const GeomSymmetricTensor<1>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator-=(const GeomSymmetricTensor<2>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator-=(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator*=(const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator*=(const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator*=(const double);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator/=(const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator/=(const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator/=(const double);

template<> SPHERAL_HOST_DEVICE bool GeomSymmetricTensor<1>::operator==(const GeomTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomSymmetricTensor<2>::operator==(const GeomTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomSymmetricTensor<3>::operator==(const GeomTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE bool GeomSymmetricTensor<1>::operator==(const GeomSymmetricTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomSymmetricTensor<2>::operator==(const GeomSymmetricTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE bool GeomSymmetricTensor<3>::operator==(const GeomSymmetricTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomSymmetricTensor<1>::diagonalElements() const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomSymmetricTensor<2>::diagonalElements() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::diagonalElements() const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::Trace() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::Trace() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::Trace() const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::Determinant() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::Determinant() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::Determinant() const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomSymmetricTensor<1>::dot(const GeomVector<1>&) const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomSymmetricTensor<2>::dot(const GeomVector<2>&) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::dot(const GeomVector<3>&) const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomSymmetricTensor<1>::dot(const GeomTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomSymmetricTensor<2>::dot(const GeomTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomSymmetricTensor<3>::dot(const GeomTensor<3>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomSymmetricTensor<1>::dot(const GeomSymmetricTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomSymmetricTensor<2>::dot(const GeomSymmetricTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomSymmetricTensor<3>::dot(const GeomSymmetricTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::doubledot(const GeomTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::doubledot(const GeomTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::doubledot(const GeomTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::doubledot(const GeomSymmetricTensor<1>&) const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::doubledot(const GeomSymmetricTensor<2>&) const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::doubledot(const GeomSymmetricTensor<3>&) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::square() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::square() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::square() const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::cube() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::cube() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::cube() const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::squareElements() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::squareElements() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::squareElements() const;

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::rotationalTransform(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::rotationalTransform(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<3>::rotationalTransform(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::maxAbsElement() const;

// Forward declare the global functions.
template<int nDim> GeomSymmetricTensor<nDim> operator*(double lhs, const GeomSymmetricTensor<nDim>& rhs);
template<int nDim> ::std::istream& operator>>(::std::istream& is, GeomSymmetricTensor<nDim>& ten);
template<int nDim> ::std::ostream& operator<<(::std::ostream& os, const GeomSymmetricTensor<nDim>& ten);

#if defined(_OPENMP) && _OPENMP >= 201107
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<1> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<1>(0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<1> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<1>(0.0) )
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<2> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<2>(0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<2> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<2>(0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<3> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<3> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#endif

}

#include "GeomSymmetricTensorInline_default.hh"

#endif
