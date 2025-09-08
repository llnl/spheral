#include "GeomTensor.hh"
#include "GeomSymmetricTensor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Geom3Vector::
Geom3Vector(const double x, const double y, const double z):
  mGeomVector(x, y, z) {
}

//------------------------------------------------------------------------------
// Implicit conversion.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Geom3Vector::
Geom3Vector(const GeomVector<3>& vec):
  mGeomVector(vec) {
}

//------------------------------------------------------------------------------
// Access elements by indicies.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
Geom3Vector::
operator()(size_type index) const {
  return mGeomVector(index);
}

SPHERAL_HOST_DEVICE inline
double&
Geom3Vector::
operator()(size_type index) {
  return mGeomVector(index);
}

//------------------------------------------------------------------------------
// Access elements by (x, y, z) notation.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
Geom3Vector::
x() const {
  return mGeomVector.x();
}

SPHERAL_HOST_DEVICE inline
double
Geom3Vector::
y() const {
  return mGeomVector.y();
}

SPHERAL_HOST_DEVICE inline
double
Geom3Vector::
z() const {
  return mGeomVector.z();
}

//------------------------------------------------------------------------------
// Set elements by (x, y, z) notation.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
void
Geom3Vector::
x(const double val) {
  mGeomVector.x(val);
}

SPHERAL_HOST_DEVICE inline
void
Geom3Vector::
y(const double val) {
  mGeomVector.y(val);
}

SPHERAL_HOST_DEVICE inline
void
Geom3Vector::
z(const double val) {
  mGeomVector.z(val);
}

//------------------------------------------------------------------------------
// dyad
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
GeomTensor<3>
Geom3Vector::
dyad(const Geom3Vector& rhs) const {
  return mGeomVector.dyad(rhs.mGeomVector);
}

//------------------------------------------------------------------------------
// selfdyad
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
GeomSymmetricTensor<3>
Geom3Vector::
selfdyad() const {
  return mGeomVector.selfdyad();
}

//------------------------------------------------------------------------------
// Multiplication by another Vector3.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
GeomTensor<3>
Geom3Vector::
operator*(const Geom3Vector& rhs) const {
  return mGeomVector*rhs.mGeomVector;
}

//------------------------------------------------------------------------------
// Multiply a scalar by a vector.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Geom3Vector
operator*(const double val, const Geom3Vector& vec) {
  return vec*val;
}

//------------------------------------------------------------------------------
// Multiplication with tensors.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Geom3Vector
operator*(const GeomTensor<3>& lhs, const Geom3Vector& rhs) {
  return Geom3Vector(lhs*GeomVector<3>(rhs.x(), rhs.y(), rhs.z()));
}

SPHERAL_HOST_DEVICE inline
Geom3Vector
operator*(const GeomSymmetricTensor<3>& lhs, const Geom3Vector& rhs) {
  return Geom3Vector(lhs*GeomVector<3>(rhs.x(), rhs.y(), rhs.z()));
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
inline
std::istream&
operator>>(std::istream& is, Geom3Vector& vec) {
  std::string parenthesis;
  is >> parenthesis;
  for (Geom3Vector::iterator elementItr = vec.begin();
       elementItr < vec.end();
       ++elementItr) {
    is >> *elementItr;
  }
  is >> parenthesis;
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
inline
std::ostream&
operator<<(std::ostream& os, const Geom3Vector& vec) {
  os << "( ";
  for (Geom3Vector::const_iterator elementItr = vec.begin();
       elementItr < vec.end();
       ++elementItr) {
    os << *elementItr << " ";
  }
  os << ")";
  return os;
}

}
