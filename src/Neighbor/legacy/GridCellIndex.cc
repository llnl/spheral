//---------------------------------Spheral++----------------------------------//
// GridCellIndex -- a simple little class to hold the gridcell indices for
// a gridcell in each dimension.  Basically an integer version of the
// GeomVector<Dimension> class.
//
// Created by:  JMO, Mon Dec 27 10:47:34 PST 1999
//----------------------------------------------------------------------------//
#include "GridCellIndex.hh"
#include "Geometry/Dimension.hh"

#include <limits.h>
#include <cmath>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Define a global function which returns a vector of GridCellIndex objects
// within a given range of grid cells.
//------------------------------------------------------------------------------
inline
vector<GridCellIndex<Dim<1> > >
GridCellIndexRange(const GridCellIndex<Dim<1> >& gridCellMin,
                   const GridCellIndex<Dim<1> >& gridCellMax) {
  int numCells = gridCellMax.xIndex() - gridCellMin.xIndex() + 1;
  CHECK(numCells > 0);

  vector<GridCellIndex<Dim<1> > > result;
  result.reserve(numCells);
  for (int ix = gridCellMin.xIndex(); ix <= gridCellMax.xIndex(); ++ix) {
    result.push_back(GridCellIndex<Dim<1> >(ix));
    CHECK((int)result.size() <= numCells);
  }

  return result;  
}

inline
vector<GridCellIndex<Dim<2> > >
GridCellIndexRange(const GridCellIndex<Dim<2> >& gridCellMin,
                   const GridCellIndex<Dim<2> >& gridCellMax) {
  int dx = gridCellMax.xIndex() - gridCellMin.xIndex();
  int dy = gridCellMax.yIndex() - gridCellMin.yIndex();
  CHECK(dx >= 0 && dy >= 0);
  int numCells = (dx + 1)*(dy + 1);
  CHECK(numCells >= 0);

  vector<GridCellIndex<Dim<2> > > result;
  result.reserve(numCells);
  for (int ix = gridCellMin.xIndex(); ix <= gridCellMax.xIndex(); ++ix) {
    for (int iy = gridCellMin.yIndex(); iy <= gridCellMax.yIndex(); ++iy) {
      result.push_back(GridCellIndex<Dim<2> >(ix, iy));
      CHECK((int)result.size() <= numCells);
    }
  }

  return result;  
}

inline
vector<GridCellIndex<Dim<3> > >
GridCellIndexRange(const GridCellIndex<Dim<3> >& gridCellMin,
                   const GridCellIndex<Dim<3> >& gridCellMax) {
  int dx = gridCellMax.xIndex() - gridCellMin.xIndex();
  int dy = gridCellMax.yIndex() - gridCellMin.yIndex();
  int dz = gridCellMax.zIndex() - gridCellMin.zIndex();
  CHECK(dx >= 0 && dy >= 0 && dz >= 0);
  int numCells = (dx + 1)*(dy + 1)*(dz + 1);
  CHECK(numCells >= 0);

  vector<GridCellIndex<Dim<3> > > result;
  result.reserve(numCells);
  for (int ix = gridCellMin.xIndex(); ix <= gridCellMax.xIndex(); ++ix) {
    for (int iy = gridCellMin.yIndex(); iy <= gridCellMax.yIndex(); ++iy) {
      for (int iz = gridCellMin.zIndex(); iz <= gridCellMax.zIndex(); ++iz) {
        result.push_back(GridCellIndex<Dim<3> >(ix, iy, iz));
        CHECK((int)result.size() <= numCells);
      }
    }
  }

  return result;  
}

}
