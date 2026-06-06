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
  SPHERAL_HOST_DEVICE GeomTensorBase(const double xx):
    mxx(xx) {}
 protected:
  double mxx = 0.0;
  SPHERAL_HOST_DEVICE GeomTensorBase() = default;
};

template<>
class GeomTensorBase<2> {
 public:
  SPHERAL_HOST_DEVICE GeomTensorBase(const double a):
    mxx(a),
    mxy(a),
    myx(a),
    myy(a) {}
  SPHERAL_HOST_DEVICE GeomTensorBase(
                 const double xx, const double xy,
                 const double yx, const double yy):
    mxx(xx),
    mxy(xy),
    myx(yx),
    myy(yy) {}
 protected:
  double mxx = 0.0;
  double mxy = 0.0;
  double myx = 0.0;
  double myy = 0.0;
  SPHERAL_HOST_DEVICE GeomTensorBase() = default;
};

template<>
class GeomTensorBase<3> {
 public:
  SPHERAL_HOST_DEVICE GeomTensorBase(const double a):
    mxx(a),
    mxy(a),
    mxz(a),
    myx(a),
    myy(a),
    myz(a),
    mzx(a),
    mzy(a), 
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
  SPHERAL_HOST_DEVICE GeomTensorBase() = default;
};

}

#endif
