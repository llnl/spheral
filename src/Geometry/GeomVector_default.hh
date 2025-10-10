//---------------------------------Spheral++----------------------------------//
// GeomVector -- Geometric Vector Class.
//
// Created by JMO, Thu Apr 15 17:39:47 PDT 1999
// Modified by:
//   July 31, 99: JMO, removing MTL based version due to problems compiling
//                MTL code with gcc2.95.
//   2004-08-19:  JMO, switch to using boost::ublas, and our vector class is
//                just an interface.
//   2004-08-23:  JMO, ublas is still too slow, so going to primitive C 
//                internal data types in accordance with suggestions from
//                Brian White
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomVector_default_hh__
#define __Spheral_GeomVector_default_hh__

#include "config.hh"

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"
#include "GeomVectorBase_default.hh"

#include <iostream>
#include "Eigen/Dense"

namespace Spheral {

template<int nDim>
class GeomVector: public GeomVectorBase<nDim> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;
  typedef Eigen::Matrix<double, nDim, 1> EigenType;

  // Useful static member data.
  static const size_type nDimensions;
  static const size_type numElements;
  static const GeomVector zero;
  static const GeomVector one;

  // Constructors.
  SPHERAL_HOST_DEVICE
  GeomVector(const double x = 0.0,
             const double y = 0.0,
             const double z = 0.0);
  template<typename Derived> GeomVector(const Eigen::MatrixBase<Derived>& vec);

  // Assignment.
  SPHERAL_HOST_DEVICE GeomVector& operator=(const double val);
  template<typename Derived> GeomVector& operator=(const Eigen::MatrixBase<Derived>& vec);

  // Allow the elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(size_type index) const;
  SPHERAL_HOST_DEVICE double& operator()(size_type index);

  // More C++ style indexing.
  SPHERAL_HOST_DEVICE double operator[](size_type index) const;
  SPHERAL_HOST_DEVICE double& operator[](size_type index);

  // Access the individual elements by (x, y, z) notation.
  SPHERAL_HOST_DEVICE double x() const;
  SPHERAL_HOST_DEVICE double y() const;
  SPHERAL_HOST_DEVICE double z() const;
  SPHERAL_HOST_DEVICE void x(const double val);
  SPHERAL_HOST_DEVICE void y(const double val);
  SPHERAL_HOST_DEVICE void z(const double val);

  // Iterator access to the raw data.
  SPHERAL_HOST_DEVICE iterator begin();
  SPHERAL_HOST_DEVICE iterator end();

  SPHERAL_HOST_DEVICE const_iterator begin() const;
  SPHERAL_HOST_DEVICE const_iterator end() const;

  // Zero the vector.
  SPHERAL_HOST_DEVICE void Zero();

  // Mathematical operators.
  SPHERAL_HOST_DEVICE GeomVector operator-() const;

  SPHERAL_HOST_DEVICE GeomVector operator+(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE GeomVector operator-(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE GeomVector operator*(const double val) const;
  SPHERAL_HOST_DEVICE GeomVector operator/(const double val) const;

  SPHERAL_HOST_DEVICE GeomVector& operator+=(const GeomVector& vec);
  SPHERAL_HOST_DEVICE GeomVector& operator-=(const GeomVector& vec);

  template<typename Derived> GeomVector& operator+=(const Eigen::MatrixBase<Derived>& vec);
  template<typename Derived> GeomVector& operator-=(const Eigen::MatrixBase<Derived>& vec);

  SPHERAL_HOST_DEVICE GeomVector& operator*=(const double val);
  SPHERAL_HOST_DEVICE GeomVector& operator/=(const double val);

  SPHERAL_HOST_DEVICE int compare(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE int compare(const double val) const;

  SPHERAL_HOST_DEVICE bool operator==(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE bool operator!=(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE bool operator<(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE bool operator>(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE bool operator<=(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE bool operator>=(const GeomVector& vec) const;

  SPHERAL_HOST_DEVICE bool operator==(const double val) const;
  SPHERAL_HOST_DEVICE bool operator!=(const double val) const;
  SPHERAL_HOST_DEVICE bool operator<(const double val) const;
  SPHERAL_HOST_DEVICE bool operator>(const double val) const;
  SPHERAL_HOST_DEVICE bool operator<=(const double val) const;
  SPHERAL_HOST_DEVICE bool operator>=(const double val) const;

  SPHERAL_HOST_DEVICE double dot(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE GeomVector<3> cross(const GeomVector& vec) const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> dyad(const GeomVector& rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor<nDim> selfdyad() const;
  SPHERAL_HOST_DEVICE GeomTensor<nDim> operator*(const GeomVector& vec) const;

  SPHERAL_HOST_DEVICE GeomVector unitVector() const;

  SPHERAL_HOST_DEVICE double magnitude() const;
  SPHERAL_HOST_DEVICE double magnitude2() const;
  SPHERAL_HOST_DEVICE double minElement() const;
  SPHERAL_HOST_DEVICE double maxElement() const;
  SPHERAL_HOST_DEVICE double maxAbsElement() const;
  SPHERAL_HOST_DEVICE double sumElements() const;
  
  //  Convert to an Eigen Vector
  EigenType eigen() const;
};

// Declare explicit specializations.
template<> SPHERAL_HOST_DEVICE GeomVector<1>::GeomVector(const double, const double, const double);
template<> SPHERAL_HOST_DEVICE GeomVector<2>::GeomVector(const double, const double, const double);
template<> SPHERAL_HOST_DEVICE GeomVector<3>::GeomVector(const double, const double, const double);

template<> SPHERAL_HOST_DEVICE GeomVector<1>& GeomVector<1>::operator=(const double val);
template<> SPHERAL_HOST_DEVICE GeomVector<2>& GeomVector<2>::operator=(const double val);
template<> SPHERAL_HOST_DEVICE GeomVector<3>& GeomVector<3>::operator=(const double val);

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::y() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<1>::z() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::z() const;

template<> SPHERAL_HOST_DEVICE void GeomVector<1>::y(const double val);
template<> SPHERAL_HOST_DEVICE void GeomVector<1>::z(const double val);
template<> SPHERAL_HOST_DEVICE void GeomVector<2>::z(const double val);

template<> SPHERAL_HOST_DEVICE void GeomVector<1>::Zero();
template<> SPHERAL_HOST_DEVICE void GeomVector<2>::Zero();
template<> SPHERAL_HOST_DEVICE void GeomVector<3>::Zero();

template<> SPHERAL_HOST_DEVICE GeomVector<1> GeomVector<1>::operator-() const;
template<> SPHERAL_HOST_DEVICE GeomVector<2> GeomVector<2>::operator-() const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomVector<3>::operator-() const;

template<> SPHERAL_HOST_DEVICE GeomVector<1>& GeomVector<1>::operator+=(const GeomVector<1>& vec);
template<> SPHERAL_HOST_DEVICE GeomVector<2>& GeomVector<2>::operator+=(const GeomVector<2>& vec);
template<> SPHERAL_HOST_DEVICE GeomVector<3>& GeomVector<3>::operator+=(const GeomVector<3>& vec);

template<> SPHERAL_HOST_DEVICE GeomVector<1>& GeomVector<1>::operator-=(const GeomVector<1>& vec);
template<> SPHERAL_HOST_DEVICE GeomVector<2>& GeomVector<2>::operator-=(const GeomVector<2>& vec);
template<> SPHERAL_HOST_DEVICE GeomVector<3>& GeomVector<3>::operator-=(const GeomVector<3>& vec);

#if defined(_OPENMP) && _OPENMP >= 201107
#pragma omp declare reduction(vecadd : GeomVector<1> : omp_out += omp_in ) initializer( omp_priv = GeomVector<1>(0.0,0.0,0.0) )
#pragma omp declare reduction(vecdif : GeomVector<1> : omp_out -= omp_in ) initializer( omp_priv = GeomVector<1>(0.0,0.0,0.0) )
#pragma omp declare reduction(vecadd : GeomVector<2> : omp_out += omp_in ) initializer( omp_priv = GeomVector<2>(0.0,0.0,0.0) )
#pragma omp declare reduction(vecdif : GeomVector<2> : omp_out -= omp_in ) initializer( omp_priv = GeomVector<2>(0.0,0.0,0.0) )
#pragma omp declare reduction(vecadd : GeomVector<3> : omp_out += omp_in ) initializer( omp_priv = GeomVector<3>(0.0,0.0,0.0) )
#pragma omp declare reduction(vecdif : GeomVector<3> : omp_out -= omp_in ) initializer( omp_priv = GeomVector<3>(0.0,0.0,0.0) )
#endif

template<> SPHERAL_HOST_DEVICE GeomVector<1>& GeomVector<1>::operator*=(const double val);
template<> SPHERAL_HOST_DEVICE GeomVector<2>& GeomVector<2>::operator*=(const double val);
template<> SPHERAL_HOST_DEVICE GeomVector<3>& GeomVector<3>::operator*=(const double val);

template<> SPHERAL_HOST_DEVICE GeomVector<1>& GeomVector<1>::operator/=(const double val);
template<> SPHERAL_HOST_DEVICE GeomVector<2>& GeomVector<2>::operator/=(const double val);
template<> SPHERAL_HOST_DEVICE GeomVector<3>& GeomVector<3>::operator/=(const double val);

template<> SPHERAL_HOST_DEVICE int GeomVector<1>::compare(const GeomVector<1>& vec) const;
template<> SPHERAL_HOST_DEVICE int GeomVector<2>::compare(const GeomVector<2>& vec) const;
template<> SPHERAL_HOST_DEVICE int GeomVector<3>::compare(const GeomVector<3>& vec) const;

template<> SPHERAL_HOST_DEVICE int GeomVector<1>::compare(const double val) const;
template<> SPHERAL_HOST_DEVICE int GeomVector<2>::compare(const double val) const;
template<> SPHERAL_HOST_DEVICE int GeomVector<3>::compare(const double val) const;

template<> SPHERAL_HOST_DEVICE bool GeomVector<1>::operator==(const GeomVector<1>& vec) const;
template<> SPHERAL_HOST_DEVICE bool GeomVector<2>::operator==(const GeomVector<2>& vec) const;
template<> SPHERAL_HOST_DEVICE bool GeomVector<3>::operator==(const GeomVector<3>& vec) const;

template<> SPHERAL_HOST_DEVICE bool GeomVector<1>::operator==(const double val) const;
template<> SPHERAL_HOST_DEVICE bool GeomVector<2>::operator==(const double val) const;
template<> SPHERAL_HOST_DEVICE bool GeomVector<3>::operator==(const double val) const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::dot(const GeomVector<1>& vec) const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::dot(const GeomVector<2>& vec) const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::dot(const GeomVector<3>& vec) const;

template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomVector<1>::cross(const GeomVector<1>& vec) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomVector<2>::cross(const GeomVector<2>& vec) const;
template<> SPHERAL_HOST_DEVICE GeomVector<3> GeomVector<3>::cross(const GeomVector<3>& vec) const;

template<> SPHERAL_HOST_DEVICE GeomTensor<1> GeomVector<1>::dyad(const GeomVector<1>& rhs) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<2> GeomVector<2>::dyad(const GeomVector<2>& rhs) const;
template<> SPHERAL_HOST_DEVICE GeomTensor<3> GeomVector<3>::dyad(const GeomVector<3>& rhs) const;

template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<1> GeomVector<1>::selfdyad() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<2> GeomVector<2>::selfdyad() const;
template<> SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> GeomVector<3>::selfdyad() const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::magnitude() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::magnitude() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::magnitude() const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::magnitude2() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::magnitude2() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::magnitude2() const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::minElement() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::minElement() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::minElement() const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::maxElement() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::maxElement() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::maxElement() const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::maxAbsElement() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::maxAbsElement() const;

template<> SPHERAL_HOST_DEVICE double GeomVector<1>::sumElements() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<2>::sumElements() const;
template<> SPHERAL_HOST_DEVICE double GeomVector<3>::sumElements() const;

// Forward declare the global functions.
template<int nDim> SPHERAL_HOST_DEVICE GeomVector<nDim> elementWiseMin(const GeomVector<nDim>& lhs,
                                                   const GeomVector<nDim>& rhs);
template<int nDim> SPHERAL_HOST_DEVICE GeomVector<nDim> elementWiseMax(const GeomVector<nDim>& lhs,
                                                   const GeomVector<nDim>& rhs);

template<> SPHERAL_HOST_DEVICE GeomVector<1> elementWiseMin(const GeomVector<1>& lhs,
                                        const GeomVector<1>& rhs);
template<> SPHERAL_HOST_DEVICE GeomVector<2> elementWiseMin(const GeomVector<2>& lhs,
                                        const GeomVector<2>& rhs);
template<> SPHERAL_HOST_DEVICE GeomVector<3> elementWiseMin(const GeomVector<3>& lhs,
                                        const GeomVector<3>& rhs);

template<> SPHERAL_HOST_DEVICE GeomVector<1> elementWiseMax(const GeomVector<1>& lhs,
                                        const GeomVector<1>& rhs);
template<> SPHERAL_HOST_DEVICE GeomVector<2> elementWiseMax(const GeomVector<2>& lhs,
                                        const GeomVector<2>& rhs);
template<> SPHERAL_HOST_DEVICE GeomVector<3> elementWiseMax(const GeomVector<3>& lhs,
                                        const GeomVector<3>& rhs);

template<int nDim> std::istream& operator>>(std::istream& is, GeomVector<nDim>& vec);
template<int nDim> std::ostream& operator<<(std::ostream& os, const GeomVector<nDim>& vec);

}

#include "GeomVectorInline_default.hh"

#endif

