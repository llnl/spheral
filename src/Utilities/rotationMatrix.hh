//------------------------------------------------------------------------------
// Determine the rotation matrix to rotate into the frame such that the given
// vector is aligned with the x axis.
//------------------------------------------------------------------------------
#ifndef rotationMatrix_HH
#define rotationMatrix_HH

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

SPHERAL_HOST_DEVICE
inline
Dim<1>::Tensor
rotationMatrix(const Dim<1>::Vector& runit) {
  REQUIRE(fuzzyEqual(runit.magnitude2(), 1.0, 1.0e-10));
  if (runit.x() > 0.0) {
    return Dim<1>::Tensor(1.0);
  } else {
    return Dim<1>::Tensor(-1.0);
  }
}

SPHERAL_HOST_DEVICE
inline
Dim<2>::Tensor
rotationMatrix(const Dim<2>::Vector& runit) {
  REQUIRE(fuzzyEqual(runit.magnitude2(), 1.0, 1.0e-10));
  typedef Dim<2>::Tensor Tensor;
  const double x = runit.x();
  const double y = runit.y();
  return Tensor( x, y,
                -y, x);
}

SPHERAL_HOST_DEVICE
inline
Dim<3>::Tensor
rotationMatrix(const Dim<3>::Vector& runit) {
  REQUIRE(fuzzyEqual(runit.magnitude2(), 1.0, 1.0e-10));
  typedef Dim<3>::Tensor Tensor;
  const double rxy = sqrt(runit.x()*runit.x() + 
                          runit.y()*runit.y());
  const double x = runit.x();
  const double y = runit.y();
  const double z = runit.z();
  if (fuzzyEqual(rxy, 0.0)) {
    // The runit vector is aligned with the z-axis, so we only want the
    // second rotation operation (R2).
    return Tensor( 0.0, 0.0, z,
                   0.0, 1.0, 0.0,
                  -z,   0.0, 0.0);
  } else {
    return Tensor( x,        y,       z,
                  -y/rxy,    x/rxy,   0.0,
                  -x*z/rxy, -y*z/rxy, rxy);
  }
}

}

#endif
