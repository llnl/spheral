//---------------------------------Spheral++----------------------------------//
// Geom3Vector -- Explicitly 3-D Geometric Vector Class.
//
// Created by JNJ, Wed Dec 21 13:13:00 PST 2005
//----------------------------------------------------------------------------//

#ifndef __Spheral__Geom3Vector_hh__
#define __Spheral__Geom3Vector_hh__

#include "Geometry/GeomVector.hh"

namespace Spheral {

// Three-dimensional vectors can exist in any dimension, since they may not 
// interact with or create forces in transverse dimensions.  Magnetic field 
// lines in a conducting fluid are a good example of three-dimensional 
// vectors embedded in a lower-dimensional space.  We need only create a 
// type that is distinct from GeomVector<3>, which we do below.
class Geom3Vector
{
  public:

  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const size_type nDimensions;
  static const size_type numElements;
  static const Geom3Vector zero;

  // Constructors.
  SPHERAL_HOST_DEVICE explicit Geom3Vector(const double x = 0.0,
                       const double y = 0.0,
                       const double z = 0.0);

  // Implicit conversion from vectors that happen to be 3-D.
  SPHERAL_HOST_DEVICE Geom3Vector(const GeomVector<3>& vec);

  // Access elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(size_type index) const;
  SPHERAL_HOST_DEVICE double& operator()(size_type index);

  // Access the individual elements by (x, y, z) notation.
  SPHERAL_HOST_DEVICE double x() const;
  SPHERAL_HOST_DEVICE double y() const;
  SPHERAL_HOST_DEVICE double z() const;
  SPHERAL_HOST_DEVICE void x(const double val);
  SPHERAL_HOST_DEVICE void y(const double val);
  SPHERAL_HOST_DEVICE void z(const double val);

  // Iterator access to the raw data.
  SPHERAL_HOST_DEVICE iterator begin() { return mGeomVector.begin(); }
  SPHERAL_HOST_DEVICE iterator end() { return mGeomVector.end(); }

  SPHERAL_HOST_DEVICE const_iterator begin() const { return mGeomVector.begin(); }
  SPHERAL_HOST_DEVICE const_iterator end() const { return mGeomVector.end(); }

  // Zero the vector.
  SPHERAL_HOST_DEVICE void Zero() { mGeomVector.Zero(); }

  // Mathematical operators.
  SPHERAL_HOST_DEVICE Geom3Vector operator-() const { return Geom3Vector(-mGeomVector); }

  SPHERAL_HOST_DEVICE Geom3Vector operator+(const Geom3Vector& vec) const { return Geom3Vector(mGeomVector + vec.mGeomVector); }
  SPHERAL_HOST_DEVICE Geom3Vector operator-(const Geom3Vector& vec) const { return Geom3Vector(mGeomVector - vec.mGeomVector); }
  SPHERAL_HOST_DEVICE Geom3Vector operator*(const double val) const { return Geom3Vector(mGeomVector * val); }
  SPHERAL_HOST_DEVICE Geom3Vector operator/(const double val) const { return Geom3Vector(mGeomVector / val); }

  SPHERAL_HOST_DEVICE Geom3Vector& operator+=(const Geom3Vector& vec) { mGeomVector += vec.mGeomVector; return *this; }
  SPHERAL_HOST_DEVICE Geom3Vector& operator-=(const Geom3Vector& vec) { mGeomVector -= vec.mGeomVector; return *this; }
  SPHERAL_HOST_DEVICE Geom3Vector& operator*=(const double val) { mGeomVector *= val; return *this; }
  SPHERAL_HOST_DEVICE Geom3Vector& operator/=(const double val) { mGeomVector /= val; return *this; }

  SPHERAL_HOST_DEVICE bool operator==(const Geom3Vector& vec) const { return mGeomVector == vec.mGeomVector; }
  SPHERAL_HOST_DEVICE bool operator!=(const Geom3Vector& vec) const { return mGeomVector != vec.mGeomVector; }
  SPHERAL_HOST_DEVICE bool operator<(const Geom3Vector& vec) const { return mGeomVector < vec.mGeomVector; }
  SPHERAL_HOST_DEVICE bool operator>(const Geom3Vector& vec) const { return mGeomVector > vec.mGeomVector; }
  SPHERAL_HOST_DEVICE bool operator<=(const Geom3Vector& vec) const { return mGeomVector <= vec.mGeomVector; }
  SPHERAL_HOST_DEVICE bool operator>=(const Geom3Vector& vec) const { return mGeomVector >= vec.mGeomVector; }

  SPHERAL_HOST_DEVICE double dot(const Geom3Vector& rhs) const { return mGeomVector.dot(rhs.mGeomVector); }
  SPHERAL_HOST_DEVICE Geom3Vector cross(const Geom3Vector& rhs) const { return Geom3Vector(mGeomVector.cross(rhs.mGeomVector)); }
  SPHERAL_HOST_DEVICE GeomTensor<3> dyad(const Geom3Vector& rhs) const;
  SPHERAL_HOST_DEVICE GeomSymmetricTensor<3> selfdyad() const;
  SPHERAL_HOST_DEVICE GeomTensor<3> operator*(const Geom3Vector& rhs) const;

  SPHERAL_HOST_DEVICE Geom3Vector unitVector() const { return Geom3Vector(mGeomVector.unitVector()); }

  SPHERAL_HOST_DEVICE double magnitude() const { return mGeomVector.magnitude(); }
  SPHERAL_HOST_DEVICE double magnitude2() const { return mGeomVector.magnitude2(); }
  SPHERAL_HOST_DEVICE double minElement() const { return mGeomVector.minElement(); }
  SPHERAL_HOST_DEVICE double maxElement() const { return mGeomVector.maxElement(); }
  SPHERAL_HOST_DEVICE double sumElements() const { return mGeomVector.sumElements(); }

private:
  GeomVector<3> mGeomVector;

};

std::istream& operator>>(std::istream& is, Geom3Vector& vec);
std::ostream& operator<<(std::ostream& os, const Geom3Vector& vec);

}


#include "Geom3VectorInline.hh"

#endif

