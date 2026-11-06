//---------------------------------Spheral++----------------------------------//
// GridCellIndexBase -- provides the dimension dependent storage for 
//                      GridCellIndex.
//
// Created by JMO, Tue Aug 16 11:13:22 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral__GridCellIndexBase__
#define __Spheral__GridCellIndexBase__

#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension> class GridCellIndexBase {};

template<>
class GridCellIndexBase<Dim<1> > {
 public:
  GridCellIndexBase():
    mx(0) {}
  GridCellIndexBase(const int x):
    mx(x) {}
  GridCellIndexBase(const GridCellIndexBase& rhs):
    mx(rhs.mx) {}
  GridCellIndexBase& operator=(const GridCellIndexBase& rhs) {
    if (this != &rhs) {
      mx = rhs.mx;
    }
    return *this;
  }
  virtual ~GridCellIndexBase() {}
 protected:
  int mx;
};

template<>
class GridCellIndexBase<Dim<2> > {
 public:
  GridCellIndexBase():
    mx(0),
    my(0) {}
  GridCellIndexBase(const int x,
                    const int y):
    mx(x),
    my(y) {}
  GridCellIndexBase(const GridCellIndexBase& rhs):
    mx(rhs.mx),
    my(rhs.my) {}
  GridCellIndexBase& operator=(const GridCellIndexBase& rhs) {
    if (this != &rhs) {
      mx = rhs.mx;
      my = rhs.my;
    }
    return *this;
  }
  virtual ~GridCellIndexBase() {}
 protected:
  int mx;
  int my;
  int mz;
};

template<>
class GridCellIndexBase<Dim<3> > {
 public:
  GridCellIndexBase():
    mx(0),
    my(0),
    mz(0) {}
  GridCellIndexBase(const int x,
                    const int y,
                    const int z):
    mx(x),
    my(y),
    mz(z) {}
  GridCellIndexBase(const GridCellIndexBase& rhs):
    mx(rhs.mx),
    my(rhs.my),
    mz(rhs.mz) {}
  GridCellIndexBase& operator=(const GridCellIndexBase& rhs) {
    if (this != &rhs) {
      mx = rhs.mx;
      my = rhs.my;
      mz = rhs.mz;
    }
    return *this;
  }
  virtual ~GridCellIndexBase() {}
 protected:
  int mx;
  int my;
  int mz;
};

}

#endif
