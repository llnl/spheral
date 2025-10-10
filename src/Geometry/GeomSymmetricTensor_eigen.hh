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
//   2017-10-09: JMO, new Eigen based version.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomSymmetricTensor_eigen_hh_
#define __Spheral_GeomSymmetricTensor_eigen_hh_

#include <iostream>

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"
#include "Geometry/EigenStruct_fwd.hh"
#include "Eigen/Dense"

namespace Spheral {

template<int nDim>
class GeomSymmetricTensor {

public:
  //--------------------------- Public Interface ---------------------------//
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef Eigen::Matrix<double, nDim, 1>                     VectorStorage;
  typedef Eigen::Matrix<double, nDim, nDim, Eigen::RowMajor> TensorStorage;
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;
  typedef EigenStruct<nDim> EigenStructType;

  // Useful static memeber data.
  static const size_type nDimensions;
  static const size_type numElements;
  static const GeomSymmetricTensor zero;
  static const GeomSymmetricTensor one;
  static const double onethird;
  static const double sqrt3;

  // Constructors.
  GeomSymmetricTensor();
  explicit GeomSymmetricTensor(const double a11);
  GeomSymmetricTensor(const double a11, const double a12,
                      const double a21, const double a22);
  GeomSymmetricTensor(const double a11, const double a12, const double a13,
                      const double a21, const double a22, const double a23,
                      const double a31, const double a32, const double a33);
  GeomSymmetricTensor(const GeomSymmetricTensor& ten);
  explicit GeomSymmetricTensor(const GeomTensor<nDim>& ten);
  GeomSymmetricTensor(const TensorStorage& ten);

  // Destructor.
  ~GeomSymmetricTensor();

  // Assignment.
  GeomSymmetricTensor& operator=(const GeomTensor<nDim>& rhs);
  GeomSymmetricTensor& operator=(const GeomSymmetricTensor& rhs);
  GeomSymmetricTensor& operator=(const TensorStorage& rhs);

  // Access the elements by indicies.
  double operator()(const size_type row, const size_type column) const;
  double& operator()(const size_type row, const size_type column);

  // More C++ style indexing.
  double operator[](size_type index) const;
  double& operator[](size_type index);

  // Access the individual elements by (x,y,z) notation.
  double xx() const;
  double xy() const;
  double xz() const;
  double yx() const;
  double yy() const;
  double yz() const;
  double zx() const;
  double zy() const;
  double zz() const;
  void xx(double val);
  void xy(double val);
  void xz(double val);
  void yx(double val);
  void yy(double val);
  void yz(double val);
  void zx(double val);
  void zy(double val);
  void zz(double val);

  // Get/set rows and columns.
  GeomVector<nDim> getRow(size_type index) const;
  GeomVector<nDim> getColumn(size_type index) const;
  void setRow(size_type index, const GeomVector<nDim>& vec);
  void setColumn(size_type index, const GeomVector<nDim>& vec);

  // Iterator access to the raw data.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // The zero and identity matrices.
  void Zero();
  void Identity();

  GeomSymmetricTensor operator-() const;

  GeomTensor<nDim> operator+(const GeomTensor<nDim>& rhs) const;
  GeomTensor<nDim> operator-(const GeomTensor<nDim>& rhs) const;
  GeomTensor<nDim> operator*(const GeomTensor<nDim>& rhs) const;

  GeomSymmetricTensor operator+(const GeomSymmetricTensor& rhs) const;
  GeomSymmetricTensor operator-(const GeomSymmetricTensor& rhs) const;
  GeomTensor<nDim>    operator*(const GeomSymmetricTensor& rhs) const;

  GeomVector<nDim> operator*(const GeomVector<nDim>& rhs) const;
  GeomSymmetricTensor operator*(const double rhs) const;
  GeomSymmetricTensor operator/(const double rhs) const;

  GeomSymmetricTensor& operator+=(const GeomSymmetricTensor& rhs);
  GeomSymmetricTensor& operator-=(const GeomSymmetricTensor& rhs);

  GeomSymmetricTensor& operator*=(const double rhs);
  GeomSymmetricTensor& operator/=(const double rhs);

  bool operator==(const GeomTensor<nDim>& rhs) const;
  bool operator!=(const GeomTensor<nDim>& rhs) const;
  bool operator<(const GeomTensor<nDim>& rhs) const;
  bool operator>(const GeomTensor<nDim>& rhs) const;
  bool operator<=(const GeomTensor<nDim>& rhs) const;
  bool operator>=(const GeomTensor<nDim>& rhs) const;

  bool operator==(const GeomSymmetricTensor& rhs) const;
  bool operator!=(const GeomSymmetricTensor& rhs) const;
  bool operator<(const GeomSymmetricTensor& rhs) const;
  bool operator>(const GeomSymmetricTensor& rhs) const;
  bool operator<=(const GeomSymmetricTensor& rhs) const;
  bool operator>=(const GeomSymmetricTensor& rhs) const;

  GeomSymmetricTensor Symmetric() const;
  GeomTensor<nDim> SkewSymmetric() const;
  GeomSymmetricTensor Transpose() const;
  GeomSymmetricTensor Inverse() const;
  GeomVector<nDim> diagonalElements() const;
  double Trace() const;
  double Determinant() const;
  GeomVector<nDim> dot(const GeomVector<nDim>& rhs) const;
  GeomTensor<nDim> dot(const GeomTensor<nDim>& rhs) const;
  GeomTensor<nDim> dot(const GeomSymmetricTensor& rhs) const;
  double doubledot(const GeomTensor<nDim>& rhs) const;
  double doubledot(const GeomSymmetricTensor& rhs) const;
  double selfDoubledot() const;

  // Return the square of this tensor (using matrix multiplication).  Note that
  // for a symmetric tensor this is guaranteed to return a symmetric product.
  GeomSymmetricTensor square() const;

  // Same idea for the cube.
  GeomSymmetricTensor cube() const;

  // Compute the "square root" of the tensor: the tensor that, 
  // when squared, equals this tensor.
  GeomSymmetricTensor sqrt() const;

  // Similarly, compute the cube root.
  GeomSymmetricTensor cuberoot() const;

  // The general version, raise to an arbitrary power.
  GeomSymmetricTensor pow(const double p) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  GeomSymmetricTensor squareElements() const;

  // Apply a rotational transform to this tensor ( R^-1 * (*this) * R ).
  void rotationalTransform(const GeomTensor<nDim>& R);

  // Return the max absolute element.
  double maxAbsElement() const;

  // A simple method for returning the eigenvalues of a tensor.
  GeomVector<nDim> eigenValues() const;
  
  // We also provide a method to retrieve the eigenvectors as an EigenStruct.
  // Note that the eigen vectors are the columns of the full tensor in the resulting
  // struct.
  EigenStructType eigenVectors() const;

  //  Access the internal Eigen type.
  TensorStorage& native();
  const TensorStorage& native() const;

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStorage mTensorData;
};

// Declare specializations.
template<> GeomSymmetricTensor<2>::GeomSymmetricTensor(const double, const double,
                                                       const double, const double);
template<> GeomSymmetricTensor<3>::GeomSymmetricTensor(const double, const double, const double,
                                                       const double, const double, const double,
                                                       const double, const double, const double);

template<> double GeomSymmetricTensor<1>::xy() const;
template<> double GeomSymmetricTensor<1>::xz() const;
template<> double GeomSymmetricTensor<1>::yx() const;
template<> double GeomSymmetricTensor<1>::yy() const;
template<> double GeomSymmetricTensor<1>::yz() const;
template<> double GeomSymmetricTensor<1>::zx() const;
template<> double GeomSymmetricTensor<1>::zy() const;
template<> double GeomSymmetricTensor<1>::zz() const;

template<> double GeomSymmetricTensor<2>::xz() const;
template<> double GeomSymmetricTensor<2>::yz() const;
template<> double GeomSymmetricTensor<2>::zx() const;
template<> double GeomSymmetricTensor<2>::zy() const;
template<> double GeomSymmetricTensor<2>::zz() const;

template<> void GeomSymmetricTensor<1>::xy(const double);
template<> void GeomSymmetricTensor<1>::xz(const double);
template<> void GeomSymmetricTensor<1>::yx(const double);
template<> void GeomSymmetricTensor<1>::yy(const double);
template<> void GeomSymmetricTensor<1>::yz(const double);
template<> void GeomSymmetricTensor<1>::zx(const double);
template<> void GeomSymmetricTensor<1>::zy(const double);
template<> void GeomSymmetricTensor<1>::zz(const double);

template<> void GeomSymmetricTensor<2>::xz(const double);
template<> void GeomSymmetricTensor<2>::yz(const double);
template<> void GeomSymmetricTensor<2>::zx(const double);
template<> void GeomSymmetricTensor<2>::zy(const double);
template<> void GeomSymmetricTensor<2>::zz(const double);

// Forward declare the global functions.
template<int nDim> GeomSymmetricTensor<nDim> operator*(double lhs, const GeomSymmetricTensor<nDim>& rhs);
template<int nDim> ::std::istream& operator>>(::std::istream& is, GeomSymmetricTensor<nDim>& ten);
template<int nDim> ::std::ostream& operator<<(::std::ostream& os, const GeomSymmetricTensor<nDim>& ten);

#ifdef _OPENMP
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<1> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<1>(0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<1> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<1>(0.0) )
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<2> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<2>(0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<2> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<2>(0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensadd : GeomSymmetricTensor<3> : omp_out += omp_in ) initializer( omp_priv = GeomSymmetricTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(symtensdif : GeomSymmetricTensor<3> : omp_out -= omp_in ) initializer( omp_priv = GeomSymmetricTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#endif
}

#include "GeomSymmetricTensorInline_eigen.hh"

#endif
