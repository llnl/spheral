//---------------------------------Spheral++----------------------------------//
// GeomSymmetricTensorBase -- provides the dimension dependent storage for
// GeomSymmetricTensor.
//
// Created by JMO, Wed Nov 10 21:33:17 PST 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral__GeomSymmetricTensorBase__
#define __Spheral__GeomSymmetricTensorBase__

#include "config.hh"

namespace Spheral {

template<int nDim> class GeomSymmetricTensorBase {};

template<>
class GeomSymmetricTensorBase<1> {
 public:
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const double a):
    mxx(a),
    myy(a),
    mzz(a) {}
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const double xx,
                                              const double yy = 0.0,
                                              const double zz = 0.0):
    mxx(xx),
    myy(yy),
    mzz(zz) {}
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase() = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const GeomSymmetricTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(GeomSymmetricTensorBase&&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase& operator=(const GeomSymmetricTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase& operator=(GeomSymmetricTensorBase&&) = default;
 protected:
  double mxx = 0.0;
  double myy = 0.0;
  double mzz = 0.0;
};

template<>
class GeomSymmetricTensorBase<2> {
 public:
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const double a):
    mxx(a),
    mxy(0.0),
    myy(a),
    mzz(a) {}
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const double xx, const double xy,
                                                               const double yy,
                                              const double zz = 0.0):
    mxx(xx),
    mxy(xy),
    myy(yy),
    mzz(zz) {}
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase() = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const GeomSymmetricTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(GeomSymmetricTensorBase&&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase& operator=(const GeomSymmetricTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase& operator=(GeomSymmetricTensorBase&&) = default;
 protected:
  double mxx = 0.0;
  double mxy = 0.0;
  double myy = 0.0;
  double mzz = 0.0;
};

template<>
class GeomSymmetricTensorBase<3> {
 public:
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const double a):
    mxx(a),
    mxy(0.0),
    mxz(0.0),
    myy(a),
    myz(0.0),
    mzz(a) {}
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const double xx, const double xy, const double xz,
                                                               const double yy, const double yz,
                                                                                const double zz):
    mxx(xx),
    mxy(xy),
    mxz(xz),
    myy(yy),
    myz(yz),
    mzz(zz) {}
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase() = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(const GeomSymmetricTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase(GeomSymmetricTensorBase&&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase& operator=(const GeomSymmetricTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomSymmetricTensorBase& operator=(GeomSymmetricTensorBase&&) = default;
 protected:
  double mxx = 0.0;
  double mxy = 0.0;
  double mxz = 0.0;
  double myy = 0.0;
  double myz = 0.0;
  double mzz = 0.0;
};

}

#endif
