//---------------------------------Spheral++----------------------------------//
// GeomTensorBase -- provides the dimension dependent storage for GeomTensor.
//
// Created by JMO, Wed Nov 10 16:40:29 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral__GeomTensorBase__
#define __Spheral__GeomTensorBase__

#include "config.hh"

namespace Spheral {

template<int nDim> class GeomTensorBase {};

template<>
class GeomTensorBase<1> {
 public:
  SPHERAL_HOST_DEVICE GeomTensorBase(const double a):
    mxx(a),
    myy(a),
    mzz(a) {}
  SPHERAL_HOST_DEVICE GeomTensorBase(const double xx,
                                     const double yy = 0.0,
                                     const double zz = 0.0):
    mxx(xx),
    myy(yy),
    mzz(zz) {}
  SPHERAL_HOST_DEVICE GeomTensorBase() = default;
  SPHERAL_HOST_DEVICE GeomTensorBase(const GeomTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase(GeomTensorBase&&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase& operator=(const GeomTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase& operator=(GeomTensorBase&&) = default;
 protected:
  double mxx = 0.0;
  double myy = 0.0;
  double mzz = 0.0;
};

template<>
class GeomTensorBase<2> {
 public:
  SPHERAL_HOST_DEVICE GeomTensorBase(const double a):
    mxx(a),
    mxy(0.0),
    myx(0.0),
    myy(a),
    mzz(a) {}
  SPHERAL_HOST_DEVICE GeomTensorBase(const double xx, const double xy,
                                     const double yx, const double yy,
                                                                       const double zz = 0.0):
    mxx(xx),
    mxy(xy),
    myx(yx),
    myy(yy),
    mzz(zz) {}
  SPHERAL_HOST_DEVICE GeomTensorBase() = default;
  SPHERAL_HOST_DEVICE GeomTensorBase(const GeomTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase(GeomTensorBase&&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase& operator=(const GeomTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase& operator=(GeomTensorBase&&) = default;
 protected:
  double mxx = 0.0;
  double mxy = 0.0;
  double myx = 0.0;
  double myy = 0.0;
  double mzz = 0.0;
};

template<>
class GeomTensorBase<3> {
 public:
  SPHERAL_HOST_DEVICE GeomTensorBase(const double a):
    mxx(a),
    mxy(0.0),
    mxz(0.0),
    myx(0.0),
    myy(a),
    myz(0.0),
    mzx(0.0),
    mzy(0.0), 
    mzz(a) {}
  SPHERAL_HOST_DEVICE GeomTensorBase(
                 const double xx, const double xy, const double xz,
                 const double yx, const double yy, const double yz,
                 const double zx, const double zy, const double zz):
    mxx(xx),
    mxy(xy),
    mxz(xz),
    myx(yx),
    myy(yy),
    myz(yz),
    mzx(zx),
    mzy(zy), 
    mzz(zz) {}
  SPHERAL_HOST_DEVICE GeomTensorBase() = default;
  SPHERAL_HOST_DEVICE GeomTensorBase(const GeomTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase(GeomTensorBase&&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase& operator=(const GeomTensorBase&) = default;
  SPHERAL_HOST_DEVICE GeomTensorBase& operator=(GeomTensorBase&&) = default;
 protected:
  double mxx = 0.0;
  double mxy = 0.0;
  double mxz = 0.0;
  double myx = 0.0;
  double myy = 0.0;
  double myz = 0.0;
  double mzx = 0.0;
  double mzy = 0.0;
  double mzz = 0.0;
};

}

#endif
