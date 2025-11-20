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
#ifndef __Spheral_GeomVector_array_hh__
#define __Spheral_GeomVector_array_hh__

#include <iostream>

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"

namespace Spheral {

template<int nDim, bool ownMemory>
class GeomVector {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const size_type nDimensions;
  static const size_type numElements;
  static const GeomVector zero;
  static const GeomVector one;

  // Constructors.
  GeomVector(const double x = 0.0,
             const double y = 0.0,
             const double z = 0.0);
  template<bool otherMemory> GeomVector(const GeomVector<nDim, otherMemory>& vec);

  // Destructor.
  ~GeomVector();

  // Assignment.
  template<bool otherMemory> GeomVector& operator=(const GeomVector<nDim, otherMemory>& vec);
  GeomVector& operator=(const double val);

  // Allow the elements by indicies.
  double operator()(size_type index) const;
  double& operator()(size_type index);

  // Access the individual elements by (x, y, z) notation.
  double x() const;
  double y() const;
  double z() const;
  void x(const double val);
  void y(const double val);
  void z(const double val);

  // Iterator access to the raw data.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Zero the vector.
  void Zero();

  // Mathematical operators.
  GeomVector<nDim, true> operator-() const;

  template<bool otherMemory> GeomVector<nDim, true> operator+(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> GeomVector<nDim, true> operator-(const GeomVector<nDim, otherMemory>& vec) const;
  GeomVector<nDim, true> operator*(const double val) const;
  GeomVector<nDim, true> operator/(const double val) const;

  template<bool otherMemory> GeomVector& operator+=(const GeomVector<nDim, otherMemory>& vec);
  template<bool otherMemory> GeomVector& operator-=(const GeomVector<nDim, otherMemory>& vec);
  GeomVector& operator*=(const double val);
  GeomVector& operator/=(const double val);

  template<bool otherMemory> int compare(const GeomVector<nDim, otherMemory>& vec) const;
  int compare(const double val) const;

  template<bool otherMemory> bool operator==(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator!=(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator< (const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator> (const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator<=(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator>=(const GeomVector<nDim, otherMemory>& vec) const;

  bool operator==(const double val) const;
  bool operator!=(const double val) const;
  bool operator< (const double val) const;
  bool operator> (const double val) const;
  bool operator<=(const double val) const;
  bool operator>=(const double val) const;

  template<bool otherMemory> double dot(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> GeomVector<3, true> cross(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> GeomTensor<nDim> dyad(const GeomVector<nDim, otherMemory>& rhs) const;
  GeomSymmetricTensor<nDim> selfdyad() const;
  template<bool otherMemory> GeomTensor<nDim> operator*(const GeomVector<nDim, otherMemory>& vec) const;

  GeomVector<nDim, true> unitVector() const;

  double magnitude() const;
  double magnitude2() const;
  double minElement() const;
  double maxElement() const;
  double maxAbsElement() const;
  double sumElements() const;

private:
  double* mData;

  //  friend class GeomVector<nDim, !ownMemory>;
};

// Declare explicit specializations.
template<> GeomVector<1, true>::GeomVector(const double, const double, const double);
template<> GeomVector<2, true>::GeomVector(const double, const double, const double);
template<> GeomVector<3, true>::GeomVector(const double, const double, const double);

template<> template<bool otherMemory> GeomVector<1, true>::GeomVector(const GeomVector<1, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<2, true>::GeomVector(const GeomVector<2, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<3, true>::GeomVector(const GeomVector<3, otherMemory>& vec);

template<> GeomVector<1, true>::~GeomVector();
template<> GeomVector<2, true>::~GeomVector();
template<> GeomVector<3, true>::~GeomVector();
template<> GeomVector<1, false>::~GeomVector();
template<> GeomVector<2, false>::~GeomVector();
template<> GeomVector<3, false>::~GeomVector();

template<> GeomVector<1> GeomVector<1, true>::operator-() const;
template<> GeomVector<2> GeomVector<2, true>::operator-() const;
template<> GeomVector<3> GeomVector<3, true>::operator-() const;
template<> GeomVector<1> GeomVector<1, false>::operator-() const;
template<> GeomVector<2> GeomVector<2, false>::operator-() const;
template<> GeomVector<3> GeomVector<3, false>::operator-() const;

template<> template<bool otherMemory> GeomVector<1, true>& GeomVector<1, true>::operator+=(const GeomVector<1, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<2, true>& GeomVector<2, true>::operator+=(const GeomVector<2, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<3, true>& GeomVector<3, true>::operator+=(const GeomVector<3, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<1, false>& GeomVector<1, false>::operator+=(const GeomVector<1, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<2, false>& GeomVector<2, false>::operator+=(const GeomVector<2, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<3, false>& GeomVector<3, false>::operator+=(const GeomVector<3, otherMemory>& vec);

template<> template<bool otherMemory> GeomVector<1, true>& GeomVector<1, true>::operator-=(const GeomVector<1, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<2, true>& GeomVector<2, true>::operator-=(const GeomVector<2, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<3, true>& GeomVector<3, true>::operator-=(const GeomVector<3, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<1, false>& GeomVector<1, false>::operator-=(const GeomVector<1, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<2, false>& GeomVector<2, false>::operator-=(const GeomVector<2, otherMemory>& vec);
template<> template<bool otherMemory> GeomVector<3, false>& GeomVector<3, false>::operator-=(const GeomVector<3, otherMemory>& vec);

template<> GeomVector<1, true>& GeomVector<1, true>::operator*=(const double val);
template<> GeomVector<2, true>& GeomVector<2, true>::operator*=(const double val);
template<> GeomVector<3, true>& GeomVector<3, true>::operator*=(const double val);
template<> GeomVector<1, false>& GeomVector<1, false>::operator*=(const double val);
template<> GeomVector<2, false>& GeomVector<2, false>::operator*=(const double val);
template<> GeomVector<3, false>& GeomVector<3, false>::operator*=(const double val);

template<> GeomVector<1, true>& GeomVector<1, true>::operator/=(const double val);
template<> GeomVector<2, true>& GeomVector<2, true>::operator/=(const double val);
template<> GeomVector<3, true>& GeomVector<3, true>::operator/=(const double val);
template<> GeomVector<1, false>& GeomVector<1, false>::operator/=(const double val);
template<> GeomVector<2, false>& GeomVector<2, false>::operator/=(const double val);
template<> GeomVector<3, false>& GeomVector<3, false>::operator/=(const double val);

template<> template<bool otherMemory> int GeomVector<1, true>::compare(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> int GeomVector<2, true>::compare(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> int GeomVector<3, true>::compare(const GeomVector<3, otherMemory>& vec) const;
template<> template<bool otherMemory> int GeomVector<1, false>::compare(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> int GeomVector<2, false>::compare(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> int GeomVector<3, false>::compare(const GeomVector<3, otherMemory>& vec) const;

template<> int GeomVector<1, true>::compare(const double val) const;
template<> int GeomVector<2, true>::compare(const double val) const;
template<> int GeomVector<3, true>::compare(const double val) const;
template<> int GeomVector<1, false>::compare(const double val) const;
template<> int GeomVector<2, false>::compare(const double val) const;
template<> int GeomVector<3, false>::compare(const double val) const;

template<> template<bool otherMemory> bool GeomVector<1, true>::operator==(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> bool GeomVector<2, true>::operator==(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> bool GeomVector<3, true>::operator==(const GeomVector<3, otherMemory>& vec) const;
template<> template<bool otherMemory> bool GeomVector<1, false>::operator==(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> bool GeomVector<2, false>::operator==(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> bool GeomVector<3, false>::operator==(const GeomVector<3, otherMemory>& vec) const;

template<> bool GeomVector<1, true>::operator==(const double val) const;
template<> bool GeomVector<2, true>::operator==(const double val) const;
template<> bool GeomVector<3, true>::operator==(const double val) const;
template<> bool GeomVector<1, false>::operator==(const double val) const;
template<> bool GeomVector<2, false>::operator==(const double val) const;
template<> bool GeomVector<3, false>::operator==(const double val) const;

template<> template<bool otherMemory> double GeomVector<1, true>::dot(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> double GeomVector<2, true>::dot(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> double GeomVector<3, true>::dot(const GeomVector<3, otherMemory>& vec) const;
template<> template<bool otherMemory> double GeomVector<1, false>::dot(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> double GeomVector<2, false>::dot(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> double GeomVector<3, false>::dot(const GeomVector<3, otherMemory>& vec) const;

template<> template<bool otherMemory> GeomVector<3> GeomVector<1, true>::cross(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> GeomVector<3> GeomVector<2, true>::cross(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> GeomVector<3> GeomVector<3, true>::cross(const GeomVector<3, otherMemory>& vec) const;
template<> template<bool otherMemory> GeomVector<3> GeomVector<1, false>::cross(const GeomVector<1, otherMemory>& vec) const;
template<> template<bool otherMemory> GeomVector<3> GeomVector<2, false>::cross(const GeomVector<2, otherMemory>& vec) const;
template<> template<bool otherMemory> GeomVector<3> GeomVector<3, false>::cross(const GeomVector<3, otherMemory>& vec) const;

template<> template<bool otherMemory> GeomTensor<1> GeomVector<1, true>::dyad(const GeomVector<1, otherMemory>& rhs) const;
template<> template<bool otherMemory> GeomTensor<2> GeomVector<2, true>::dyad(const GeomVector<2, otherMemory>& rhs) const;
template<> template<bool otherMemory> GeomTensor<3> GeomVector<3, true>::dyad(const GeomVector<3, otherMemory>& rhs) const;
template<> template<bool otherMemory> GeomTensor<1> GeomVector<1, false>::dyad(const GeomVector<1, otherMemory>& rhs) const;
template<> template<bool otherMemory> GeomTensor<2> GeomVector<2, false>::dyad(const GeomVector<2, otherMemory>& rhs) const;
template<> template<bool otherMemory> GeomTensor<3> GeomVector<3, false>::dyad(const GeomVector<3, otherMemory>& rhs) const;

template<> GeomSymmetricTensor<1> GeomVector<1, true>::selfdyad() const;
template<> GeomSymmetricTensor<2> GeomVector<2, true>::selfdyad() const;
template<> GeomSymmetricTensor<3> GeomVector<3, true>::selfdyad() const;
template<> GeomSymmetricTensor<1> GeomVector<1, false>::selfdyad() const;
template<> GeomSymmetricTensor<2> GeomVector<2, false>::selfdyad() const;
template<> GeomSymmetricTensor<3> GeomVector<3, false>::selfdyad() const;

template<> double GeomVector<1, true>::maxAbsElement() const;
template<> double GeomVector<2, true>::maxAbsElement() const;
template<> double GeomVector<3, true>::maxAbsElement() const;
template<> double GeomVector<1, false>::maxAbsElement() const;
template<> double GeomVector<2, false>::maxAbsElement() const;
template<> double GeomVector<3, false>::maxAbsElement() const;

// Forward declare the global functions.
template<bool ownMemory, bool otherMemory> GeomVector<1> elementWiseMin(const GeomVector<1, ownMemory>& lhs,
                                                                        const GeomVector<1, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<2> elementWiseMin(const GeomVector<2, ownMemory>& lhs,
                                                                        const GeomVector<2, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<3> elementWiseMin(const GeomVector<3, ownMemory>& lhs,
                                                                        const GeomVector<3, otherMemory>& rhs);

template<bool ownMemory, bool otherMemory> GeomVector<1> elementWiseMax(const GeomVector<1, ownMemory>& lhs,
                                                                        const GeomVector<1, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<2> elementWiseMax(const GeomVector<2, ownMemory>& lhs,
                                                                        const GeomVector<2, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<3> elementWiseMax(const GeomVector<3, ownMemory>& lhs,
                                                                        const GeomVector<3, otherMemory>& rhs);

template<int nDim, bool ownMemory> std::istream& operator>>(std::istream& is, GeomVector<nDim, ownMemory>& vec);
template<int nDim, bool ownMemory> std::ostream& operator<<(std::ostream& os, const GeomVector<nDim, ownMemory>& vec);

}

#include "GeomVectorInline_array.hh"

#endif

