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
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;
  typedef Eigen::Matrix<double, nDim, nDim> EigenType;
  typedef EigenStruct<nDim> EigenStructType;

  // Useful static memeber data.
  static constexpr size_type nDimensions = nDim;
  static constexpr size_type numElements = (nDim * (nDim+1)) / 2;
  static const GeomSymmetricTensor zero;
  static const GeomSymmetricTensor one;
  static constexpr double onethird = 1.0/3.0;
  static const double sqrt3;

  // Constructors.
  GeomSymmetricTensor();
  explicit GeomSymmetricTensor(const double a11);
  GeomSymmetricTensor(const double /*a11*/, const double /*a12*/,
                      const double /*a21*/, const double /*a22*/);
  GeomSymmetricTensor(const double /*a11*/, const double /*a12*/, const double /*a13*/,
                      const double /*a21*/, const double /*a22*/, const double /*a23*/,
                      const double /*a31*/, const double /*a32*/, const double /*a33*/);
  GeomSymmetricTensor(const GeomSymmetricTensor& ten);
  explicit GeomSymmetricTensor(const GeomTensor<nDim>& ten);
  template<typename Derived> GeomSymmetricTensor(const Eigen::MatrixBase<Derived>& ten);

  // Destructor.
  ~GeomSymmetricTensor();

  // Assignment.
  GeomSymmetricTensor& operator=(const GeomTensor<nDim>& rhs);
  GeomSymmetricTensor& operator=(const GeomSymmetricTensor& rhs);
  template<typename Derived> GeomSymmetricTensor& operator=(const Eigen::MatrixBase<Derived>& rhs);

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
  // GeomSymmetricTensor operator+(const double rhs) const;
  // GeomSymmetricTensor operator-(const double rhs) const;
  GeomSymmetricTensor operator*(const double rhs) const;
  GeomSymmetricTensor operator/(const double rhs) const;

  GeomSymmetricTensor& operator+=(const GeomSymmetricTensor& rhs);
  GeomSymmetricTensor& operator-=(const GeomSymmetricTensor& rhs);

  template<typename Derived> GeomSymmetricTensor& operator+=(const Eigen::MatrixBase<Derived>& rhs);
  template<typename Derived> GeomSymmetricTensor& operator-=(const Eigen::MatrixBase<Derived>& rhs);

  // GeomSymmetricTensor& operator+=(const double rhs);
  // GeomSymmetricTensor& operator-=(const double rhs);
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

  // bool operator==(const double rhs) const;
  // bool operator!=(const double rhs) const;
  // bool operator<(const double rhs) const;
  // bool operator>(const double rhs) const;
  // bool operator<=(const double rhs) const;
  // bool operator>=(const double rhs) const;

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

  // A simple method for returning the eigenvalues of a tensor.
  GeomVector<nDim> eigenValues() const;
  
  // Apply a rotational transform to this tensor ( R^-1 * (*this) * R ).
  void rotationalTransform(const GeomTensor<nDim>& R);

  // We also provide a method to retrieve the eigenvectors as an EigenStruct.
  // Note that the eigen vectors are the columns of the full tensor in the resulting
  // struct.
  EigenStructType eigenVectors() const;

  // Return the max absolute element.
  double maxAbsElement() const;

  //  Convert to an Eigen Vector
  EigenType eigen() const;

private:
  //--------------------------- Private Interface ---------------------------//
  size_type elementIndex(const size_type row, const size_type column) const;
};

// Declare specializations.
template<> GeomVector<1> GeomSymmetricTensor<1>::eigenValues() const;
template<> GeomVector<2> GeomSymmetricTensor<2>::eigenValues() const;
template<> GeomVector<3> GeomSymmetricTensor<3>::eigenValues() const;

template<> EigenStruct<1> GeomSymmetricTensor<1>::eigenVectors() const;
template<> EigenStruct<2> GeomSymmetricTensor<2>::eigenVectors() const;
template<> EigenStruct<3> GeomSymmetricTensor<3>::eigenVectors() const;

template<> GeomSymmetricTensor<2>::GeomSymmetricTensor(const double, const double,
                                                       const double, const double);
template<> GeomSymmetricTensor<3>::GeomSymmetricTensor(const double, const double, const double,
                                                       const double, const double, const double,
                                                       const double, const double, const double);

template<> GeomSymmetricTensor<1>::GeomSymmetricTensor(const GeomTensor<1>&);
template<> GeomSymmetricTensor<2>::GeomSymmetricTensor(const GeomTensor<2>&);
template<> GeomSymmetricTensor<3>::GeomSymmetricTensor(const GeomTensor<3>&);

template<> GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator=(const GeomTensor<1>&);
template<> GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator=(const GeomTensor<2>&);
template<> GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator=(const GeomTensor<3>&);

template<> GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator=(const GeomSymmetricTensor<1>&);
template<> GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator=(const GeomSymmetricTensor<2>&);
template<> GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator=(const GeomSymmetricTensor<3>&);

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

template<> GeomVector<1> GeomSymmetricTensor<1>::getRow(const GeomSymmetricTensor<1>::size_type) const;
template<> GeomVector<2> GeomSymmetricTensor<2>::getRow(const GeomSymmetricTensor<2>::size_type) const;
template<> GeomVector<3> GeomSymmetricTensor<3>::getRow(const GeomSymmetricTensor<3>::size_type) const;

template<> GeomVector<1> GeomSymmetricTensor<1>::getColumn(const GeomSymmetricTensor<1>::size_type) const;
template<> GeomVector<2> GeomSymmetricTensor<2>::getColumn(const GeomSymmetricTensor<2>::size_type) const;
template<> GeomVector<3> GeomSymmetricTensor<3>::getColumn(const GeomSymmetricTensor<3>::size_type) const;

template<> void GeomSymmetricTensor<1>::setRow(const GeomSymmetricTensor<1>::size_type, const GeomVector<1>&);
template<> void GeomSymmetricTensor<2>::setRow(const GeomSymmetricTensor<2>::size_type, const GeomVector<2>&);
template<> void GeomSymmetricTensor<3>::setRow(const GeomSymmetricTensor<3>::size_type, const GeomVector<3>&);

template<> void GeomSymmetricTensor<1>::setColumn(const GeomSymmetricTensor<1>::size_type, const GeomVector<1>&);
template<> void GeomSymmetricTensor<2>::setColumn(const GeomSymmetricTensor<2>::size_type, const GeomVector<2>&);
template<> void GeomSymmetricTensor<3>::setColumn(const GeomSymmetricTensor<3>::size_type, const GeomVector<3>&);

template<> GeomSymmetricTensor<1> GeomSymmetricTensor<1>::operator-() const;
template<> GeomSymmetricTensor<2> GeomSymmetricTensor<2>::operator-() const;
template<> GeomSymmetricTensor<3> GeomSymmetricTensor<3>::operator-() const;

template<> GeomSymmetricTensor<1> GeomSymmetricTensor<1>::operator*(const double) const;
template<> GeomSymmetricTensor<2> GeomSymmetricTensor<2>::operator*(const double) const;
template<> GeomSymmetricTensor<3> GeomSymmetricTensor<3>::operator*(const double) const;

template<> GeomSymmetricTensor<1> GeomSymmetricTensor<1>::operator/(const double) const;
template<> GeomSymmetricTensor<2> GeomSymmetricTensor<2>::operator/(const double) const;
template<> GeomSymmetricTensor<3> GeomSymmetricTensor<3>::operator/(const double) const;

template<> GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator+=(const GeomSymmetricTensor<1>&);
template<> GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator+=(const GeomSymmetricTensor<2>&);
template<> GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator+=(const GeomSymmetricTensor<3>&);

template<> GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator-=(const GeomSymmetricTensor<1>&);
template<> GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator-=(const GeomSymmetricTensor<2>&);
template<> GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator-=(const GeomSymmetricTensor<3>&);

template<> GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator*=(const double);
template<> GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator*=(const double);
template<> GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator*=(const double);

template<> GeomSymmetricTensor<1>& GeomSymmetricTensor<1>::operator/=(const double);
template<> GeomSymmetricTensor<2>& GeomSymmetricTensor<2>::operator/=(const double);
template<> GeomSymmetricTensor<3>& GeomSymmetricTensor<3>::operator/=(const double);

template<> bool GeomSymmetricTensor<1>::operator==(const GeomTensor<1>&) const;
template<> bool GeomSymmetricTensor<2>::operator==(const GeomTensor<2>&) const;
template<> bool GeomSymmetricTensor<3>::operator==(const GeomTensor<3>&) const;

template<> bool GeomSymmetricTensor<1>::operator==(const GeomSymmetricTensor<1>&) const;
template<> bool GeomSymmetricTensor<2>::operator==(const GeomSymmetricTensor<2>&) const;
template<> bool GeomSymmetricTensor<3>::operator==(const GeomSymmetricTensor<3>&) const;

template<> GeomVector<1> GeomSymmetricTensor<1>::diagonalElements() const;
template<> GeomVector<2> GeomSymmetricTensor<2>::diagonalElements() const;
template<> GeomVector<3> GeomSymmetricTensor<3>::diagonalElements() const;

template<> double GeomSymmetricTensor<1>::Trace() const;
template<> double GeomSymmetricTensor<2>::Trace() const;
template<> double GeomSymmetricTensor<3>::Trace() const;

template<> double GeomSymmetricTensor<1>::Determinant() const;
template<> double GeomSymmetricTensor<2>::Determinant() const;
template<> double GeomSymmetricTensor<3>::Determinant() const;

template<> GeomVector<1> GeomSymmetricTensor<1>::dot(const GeomVector<1>&) const;
template<> GeomVector<2> GeomSymmetricTensor<2>::dot(const GeomVector<2>&) const;
template<> GeomVector<3> GeomSymmetricTensor<3>::dot(const GeomVector<3>&) const;

template<> GeomTensor<1> GeomSymmetricTensor<1>::dot(const GeomTensor<1>&) const;
template<> GeomTensor<2> GeomSymmetricTensor<2>::dot(const GeomTensor<2>&) const;
template<> GeomTensor<3> GeomSymmetricTensor<3>::dot(const GeomTensor<3>&) const;
template<> GeomTensor<1> GeomSymmetricTensor<1>::dot(const GeomSymmetricTensor<1>&) const;
template<> GeomTensor<2> GeomSymmetricTensor<2>::dot(const GeomSymmetricTensor<2>&) const;
template<> GeomTensor<3> GeomSymmetricTensor<3>::dot(const GeomSymmetricTensor<3>&) const;

template<> double GeomSymmetricTensor<1>::doubledot(const GeomTensor<1>&) const;
template<> double GeomSymmetricTensor<2>::doubledot(const GeomTensor<2>&) const;
template<> double GeomSymmetricTensor<3>::doubledot(const GeomTensor<3>&) const;

template<> double GeomSymmetricTensor<1>::doubledot(const GeomSymmetricTensor<1>&) const;
template<> double GeomSymmetricTensor<2>::doubledot(const GeomSymmetricTensor<2>&) const;
template<> double GeomSymmetricTensor<3>::doubledot(const GeomSymmetricTensor<3>&) const;

template<> GeomSymmetricTensor<1> GeomSymmetricTensor<1>::square() const;
template<> GeomSymmetricTensor<2> GeomSymmetricTensor<2>::square() const;
template<> GeomSymmetricTensor<3> GeomSymmetricTensor<3>::square() const;

template<> GeomSymmetricTensor<1> GeomSymmetricTensor<1>::cube() const;
template<> GeomSymmetricTensor<2> GeomSymmetricTensor<2>::cube() const;
template<> GeomSymmetricTensor<3> GeomSymmetricTensor<3>::cube() const;

template<> GeomSymmetricTensor<1> GeomSymmetricTensor<1>::squareElements() const;
template<> GeomSymmetricTensor<2> GeomSymmetricTensor<2>::squareElements() const;
template<> GeomSymmetricTensor<3> GeomSymmetricTensor<3>::squareElements() const;

template<> void GeomSymmetricTensor<1>::rotationalTransform(const GeomTensor<1>&);
template<> void GeomSymmetricTensor<2>::rotationalTransform(const GeomTensor<2>&);
template<> void GeomSymmetricTensor<3>::rotationalTransform(const GeomTensor<3>&);

template<> double GeomSymmetricTensor<1>::maxAbsElement() const;
template<> double GeomSymmetricTensor<2>::maxAbsElement() const;
template<> double GeomSymmetricTensor<3>::maxAbsElement() const;

template<> const GeomSymmetricTensor<1> GeomSymmetricTensor<1>::zero;
template<> const GeomSymmetricTensor<2> GeomSymmetricTensor<2>::zero;
template<> const GeomSymmetricTensor<3> GeomSymmetricTensor<3>::zero;

template<> const GeomSymmetricTensor<1> GeomSymmetricTensor<1>::one;
template<> const GeomSymmetricTensor<2> GeomSymmetricTensor<2>::one;
template<> const GeomSymmetricTensor<3> GeomSymmetricTensor<3>::one;

template<> const double GeomSymmetricTensor<1>::sqrt3;
template<> const double GeomSymmetricTensor<2>::sqrt3;
template<> const double GeomSymmetricTensor<3>::sqrt3;

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
