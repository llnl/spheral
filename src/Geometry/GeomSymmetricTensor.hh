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
#include "Threading/GPUUtils.hh"

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
  using VectorType = GeomVector<nDim>;
  using SymTensorType = GeomSymmetricTensor<nDim>;
  using TensorType = GeomTensor<nDim>;

  // Useful stuff known at compile time
  SPHERAL_HOST_DEVICE static constexpr size_type nDimensions()            { return nDim; }
  SPHERAL_HOST_DEVICE static constexpr size_type numElements()            { return 3u + nDim*(nDim - 1u)/2u; }
  SPHERAL_HOST_DEVICE static constexpr GeomSymmetricTensor zero()         { return GeomSymmetricTensor(); }
  SPHERAL_HOST_DEVICE static constexpr GeomSymmetricTensor one()          { GeomSymmetricTensor result; result.mxx = 1.0; result.myy = 1.0; result.mzz = 1.0; return result; }
  SPHERAL_HOST_DEVICE static constexpr double onethird()                  { return 1.0/3.0; }
  SPHERAL_HOST_DEVICE static constexpr double sqrt3()                     { return std::sqrt(3.0); }

  // Constructors.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor() = default;
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const double a);                                          // Any dimension
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const double a11,                                         // 1D
                                                   const double a22,
                                                   const double a33) requires (nDim == 1u);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const double a11, const double a12,                                // 2D
                                          const double a21, const double a22,
                                          const double a33 = 0.0) requires (nDim == 2u);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const double a11, const double a12, const double a13,              // 3D
                                          const double a21, const double a22, const double a23,
                                          const double a31, const double a32, const double a33) requires (nDim == 3u);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const GeomSymmetricTensor& rhs) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(GeomSymmetricTensor&& rhs) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const TensorType& rhs);
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const GeomSymmetricTensor<3>& rhs) requires (nDim < 3u);
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const GeomTensor<3>& rhs) requires (nDim < 3u);
  template<typename Derived> GeomSymmetricTensor(const Eigen::MatrixBase<Derived>& ten);

  // Construct with the given diagonal elements
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const VectorType& diagonal);
  SPHERAL_HOST_DEVICE explicit GeomSymmetricTensor(const GeomVector<3>& diagonal) requires (nDim < 3u);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const VectorType& diagonal,    const double minval, const double maxval);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor(const GeomVector<3>& diagonal, const double minval, const double maxval) requires (nDim < 3u);

  // Assignment.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator=(const GeomSymmetricTensor& rhs) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator=(GeomSymmetricTensor&& rhs) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator=(const TensorType& rhs);
  template<typename Derived> GeomSymmetricTensor& operator=(const Eigen::MatrixBase<Derived>& rhs);

  // Access the elements by indicies.
  SPHERAL_HOST_DEVICE double  operator()(const size_type row, const size_type column) const;
  SPHERAL_HOST_DEVICE double& operator()(const size_type row, const size_type column);

  // More C++ style indexing.
  SPHERAL_HOST_DEVICE double  operator[](size_type index) const;
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
  SPHERAL_HOST_DEVICE VectorType getRow(size_type index) const;
  SPHERAL_HOST_DEVICE VectorType getColumn(size_type index) const;
  SPHERAL_HOST_DEVICE void setRow(size_type index, const VectorType& vec);
  SPHERAL_HOST_DEVICE void setColumn(size_type index, const VectorType& vec);

  // Iterator access to the raw data.
  SPHERAL_HOST_DEVICE iterator begin();
  SPHERAL_HOST_DEVICE iterator end();

  SPHERAL_HOST_DEVICE const_iterator begin() const;
  SPHERAL_HOST_DEVICE const_iterator end() const;

  // The zero and identity matrices.
  SPHERAL_HOST_DEVICE void Zero();
  SPHERAL_HOST_DEVICE void Identity();

  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator-() const;

  SPHERAL_HOST_DEVICE TensorType operator+(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE TensorType operator-(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE TensorType operator*(const TensorType& rhs) const;

  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator+(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator-(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE TensorType          operator*(const GeomSymmetricTensor& rhs) const;

  SPHERAL_HOST_DEVICE VectorType          operator*(const VectorType& rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator*(const double rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor operator/(const double rhs) const;

  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator+=(const GeomSymmetricTensor& rhs);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator-=(const GeomSymmetricTensor& rhs);

  template<typename Derived> GeomSymmetricTensor& operator+=(const Eigen::MatrixBase<Derived>& rhs);
  template<typename Derived> GeomSymmetricTensor& operator-=(const Eigen::MatrixBase<Derived>& rhs);

  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator*=(const double rhs);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator/=(const double rhs);

  SPHERAL_HOST_DEVICE bool operator==(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<=(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE bool operator>=(const TensorType& rhs) const;

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
  SPHERAL_HOST_DEVICE TensorType SkewSymmetric() const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor Transpose() const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor Inverse() const;
  SPHERAL_HOST_DEVICE VectorType diagonalElements() const;
  SPHERAL_HOST_DEVICE double Trace() const;
  SPHERAL_HOST_DEVICE double Determinant() const;
  SPHERAL_HOST_DEVICE VectorType dot(const VectorType& rhs) const;
  SPHERAL_HOST_DEVICE TensorType dot(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE TensorType dot(const GeomSymmetricTensor& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot(const SymTensorType& rhs) const;
  SPHERAL_HOST_DEVICE double selfDoubledot() const;

  // Return the square of this tensor (using matrix multiplication).  Note that
  // for a symmetric tensor this is guaranteed to return a symmetric product.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor square() const;

  // Compute the "square root" of the tensor: the tensor that, 
  // when squared, equals this tensor.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor sqrt() const;

  // The general version, raise to an arbitrary power.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor pow(const double p) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  SPHERAL_HOST_DEVICE GeomSymmetricTensor squareElements() const;

  // A simple method for returning the eigenvalues of a tensor.
  SPHERAL_HOST_DEVICE VectorType eigenValues() const;
  
  // We also provide a method to retrieve the eigenvectors as an EigenStruct.
  // Note that the eigen vectors are the columns of the full tensor in the resulting
  // struct.
  SPHERAL_HOST_DEVICE EigenStructType eigenVectors() const;

  // Apply a rotational transform to this tensor ( R^-1 * (*this) * R ).
  SPHERAL_HOST_DEVICE void rotationalTransform(const TensorType& R);
  SPHERAL_HOST_DEVICE void rotationalTransform(const GeomTensor<3>& R) requires (nDim < 3u);

  // Return the max absolute element.
  SPHERAL_HOST_DEVICE double maxAbsElement() const;

  // Enforce min/max eigenvalues
  SPHERAL_HOST_DEVICE GeomSymmetricTensor enforceMinEigenValue(const double& minvalue) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor enforceMaxEigenValue(const double& maxvalue) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor clampEigenValues(const double& minvalue, const double& maxvalue) const;

  //  Convert to an Eigen Vector
  EigenType eigen() const;

  // A subset of operations are specially defined for full 3D results
  SPHERAL_HOST_DEVICE GeomVector<3> diagonalElements3D() const;
  SPHERAL_HOST_DEVICE double Trace3D() const;
  SPHERAL_HOST_DEVICE double Determinant3D() const;
  SPHERAL_HOST_DEVICE double doubledot3D(const TensorType& rhs) const;
  SPHERAL_HOST_DEVICE double doubledot3D(const SymTensorType& rhs) const;
  SPHERAL_HOST_DEVICE double selfDoubledot3D() const;
  SPHERAL_HOST_DEVICE double maxAbsElement3D() const;
  SPHERAL_HOST_DEVICE GeomVector<3> eigenValues3D() const;
  SPHERAL_HOST_DEVICE EigenStruct<3> eigenVectors3D() const;
  SPHERAL_HOST_DEVICE GeomVector<3> dot(const GeomVector<3>& rhs) const                   requires (nDim < 3u);  // Screen out as redundant with normal vector dot above in 3d
  SPHERAL_HOST_DEVICE GeomVector<3> operator*(const GeomVector<3>& rhs) const             requires (nDim < 3u);
  SPHERAL_HOST_DEVICE GeomSymmetricTensor& operator+=(const GeomSymmetricTensor<3>& rhs)  requires (nDim < 3u);

  // Support for atomic operations on device
  template<typename Op> SPHERAL_HOST_DEVICE void atomicOp(const TensorType& rhs);
  template<typename Op> SPHERAL_HOST_DEVICE void atomicOp(const SymTensorType& rhs);

  SPHERAL_HOST_DEVICE void atomicAdd(const TensorType& rhs)    { atomicOp<GPUUtils::AtomicAddOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicAdd(const SymTensorType& rhs) { atomicOp<GPUUtils::AtomicAddOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicSub(const TensorType& rhs)    { atomicOp<GPUUtils::AtomicSubOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicSub(const SymTensorType& rhs) { atomicOp<GPUUtils::AtomicSubOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicMax(const TensorType& rhs)    { atomicOp<GPUUtils::AtomicMaxOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicMax(const SymTensorType& rhs) { atomicOp<GPUUtils::AtomicMaxOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicMin(const TensorType& rhs)    { atomicOp<GPUUtils::AtomicMinOp>(rhs); }
  SPHERAL_HOST_DEVICE void atomicMin(const SymTensorType& rhs) { atomicOp<GPUUtils::AtomicMinOp>(rhs); }

private:
  //--------------------------- Private Interface ---------------------------//
  SPHERAL_HOST_DEVICE size_type elementIndex(const size_type row, const size_type column) const;
};

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>::GeomSymmetricTensor(const double, const double, const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>::GeomSymmetricTensor(const double, const double,
                                                                           const double, const double,
                                                                           const double);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>::GeomSymmetricTensor(const double, const double, const double,
                                                                           const double, const double, const double,
                                                                           const double, const double, const double);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>::GeomSymmetricTensor(const GeomSymmetricTensor<3>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>::GeomSymmetricTensor(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator=(const GeomTensor<1>& rhs);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator=(const GeomTensor<2>& rhs);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator=(const GeomTensor<3>& rhs);

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::xy() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::xz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::yx() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::yz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::zx() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::zy() const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::xz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::yz() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::zx() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::zy() const;

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::xy(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::xz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::yx(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::yz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::zx(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::zy(const double);

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::xz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::yz(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::zx(const double);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::zy(const double);

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

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::sqrt() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::sqrt() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::sqrt() const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::pow(const double) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::pow(const double) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::pow(const double) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::squareElements() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::squareElements() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::squareElements() const;

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomSymmetricTensor<1>::eigenValues() const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomSymmetricTensor<2>::eigenValues() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::eigenValues() const;

template<> SPHERAL_HOST_DEVICE EigenStruct<1> GeomSymmetricTensor<1>::eigenVectors() const;
template<> SPHERAL_HOST_DEVICE EigenStruct<2> GeomSymmetricTensor<2>::eigenVectors() const;
template<> SPHERAL_HOST_DEVICE EigenStruct<3> GeomSymmetricTensor<3>::eigenVectors() const;

template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<1>::rotationalTransform(const GeomTensor<1>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<2>::rotationalTransform(const GeomTensor<2>&);
template<> SPHERAL_HOST_DEVICE void GeomSymmetricTensor<3>::rotationalTransform(const GeomTensor<3>&);

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::maxAbsElement() const;

template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<1>::maxAbsElement3D() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<2>::maxAbsElement3D() const;
template<> SPHERAL_HOST_DEVICE double GeomSymmetricTensor<3>::maxAbsElement3D() const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator+=(const GeomSymmetricTensor<3>&);
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator+=(const GeomSymmetricTensor<3>&);

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::enforceMinEigenValue(const double&) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::enforceMinEigenValue(const double&) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::enforceMinEigenValue(const double&) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::enforceMaxEigenValue(const double&) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::enforceMaxEigenValue(const double&) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::enforceMaxEigenValue(const double&) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomSymmetricTensor<1>::clampEigenValues(const double&, const double&) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomSymmetricTensor<2>::clampEigenValues(const double&, const double&) const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomSymmetricTensor<3>::clampEigenValues(const double&, const double&) const;

template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<1>::eigenValues3D() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<2>::eigenValues3D() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomSymmetricTensor<3>::eigenValues3D() const;

template<> SPHERAL_HOST_DEVICE EigenStruct<3> GeomSymmetricTensor<1>::eigenVectors3D() const;
template<> SPHERAL_HOST_DEVICE EigenStruct<3> GeomSymmetricTensor<2>::eigenVectors3D() const;
template<> SPHERAL_HOST_DEVICE EigenStruct<3> GeomSymmetricTensor<3>::eigenVectors3D() const;

// Forward declare the global functions.
template<int nDim> SPHERAL_HOST_DEVICE GeomSymmetricTensor<nDim> operator*(double lhs, const GeomSymmetricTensor<nDim>& rhs);
template<int nDim> ::std::istream& operator>>(::std::istream& is, GeomSymmetricTensor<nDim>& ten);
template<int nDim> ::std::ostream& operator<<(::std::ostream& os, const GeomSymmetricTensor<nDim>& ten);

#if defined(_OPENMP) && _OPENMP >= 201107
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<1> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<1>(0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<1> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<1>(0.0,0.0,0.0) )
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<2> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<2>(0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<2> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<2>(0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<3> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<3> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#endif

}

#include "GeomSymmetricTensorInline.hh"

#endif
