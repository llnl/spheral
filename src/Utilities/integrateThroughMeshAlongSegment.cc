//---------------------------------Spheral++----------------------------------//
// integrateThroughMeshAlongSegment
//
// Return the result of integrating a quantity along a line segment.
// The quantity here is assumed to be represented a values in a vector<Value>,
// where the vector<Value> is the value of the quantity in a series of cartesian
// cells whose box is defined by by xmin, xmax, and ncells.
//
// Created by JMO, Wed Feb  3 16:03:46 PST 2010
//----------------------------------------------------------------------------//
#include "integrateThroughMeshAlongSegment.hh"
#include "lineSegmentIntersections.hh"
#include "safeInv.hh"
#include "testBoxIntersection.hh"
#include "DataTypeTraits.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/binFieldList2Lattice.hh"

#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Find the index of the cell boundary just left or right (1-D) of the given
// coordinate.
//------------------------------------------------------------------------------
inline
unsigned
leftCellBoundary(const double xi,
                 const double xmin,
                 const double xmax,
                 const unsigned ncells) {
  return unsigned(max(0.0, min(1.0, (xi - xmin)/(xmax - xmin)))*ncells);
}

inline
unsigned
rightCellBoundary(const double xi,
                  const double xmin,
                  const double xmax,
                  const unsigned ncells) {
  const double f = max(0.0, min(1.0, (xi - xmin)/(xmax - xmin)))*ncells;
  return f - unsigned(f) > 1.0e-10 ? unsigned(f) + 1U : unsigned(f);
}
}
  
//------------------------------------------------------------------------------
// Find the points of intersection with the mesh planes for the given segment.
//------------------------------------------------------------------------------
#ifdef SPHERAL1D
// 1-D.
vector<Dim<1>::Vector>
findIntersections(const Dim<1>::Vector& xmin,
		  const Dim<1>::Vector& xmax,
		  const vector<unsigned>& ncells,
		  const Dim<1>::Vector& s0,
		  const Dim<1>::Vector& s1) {
  REQUIRE(xmin.x() < xmax.x());
  REQUIRE(ncells.size() == 1);
  REQUIRE(ncells[0] > 0);

  // Find the min and max bounding mesh planes indicies.
  typedef Dim<1>::Vector Vector;
  const Vector smin = elementWiseMin(s0, s1);
  const Vector smax = elementWiseMax(s0, s1);
  const unsigned ixmin = rightCellBoundary(smin(0), xmin(0), xmax(0), ncells[0]);
  const unsigned ixmax =  leftCellBoundary(smax(0), xmin(0), xmax(0), ncells[0]);
  CHECK(ixmin <= ncells[0]);
  CHECK(ixmax <= ncells[0]);

  // The intersections are just the intermediate planes.
  vector<Vector> result;
  const double xstep = (xmax.x() - xmin.x())/ncells[0];
  for (unsigned iplane = ixmin; iplane < ixmax; ++iplane) {
    result.push_back(xmin + Dim<1>::Vector(iplane*xstep));
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    if (result.size() > 0) {
      ENSURE(result.front().x() >= smin.x());
      ENSURE(result.back().x() <= smax.x());
    }
  }
  END_CONTRACT_SCOPE

  return result;
}
#endif
  
//------------------------------------------------------------------------------
// 2-D.
#ifdef SPHERAL2D
vector<Dim<2>::Vector> 
findIntersections(const Dim<2>::Vector& xmin,
		  const Dim<2>::Vector& xmax,
		  const vector<unsigned>& ncells,
		  const Dim<2>::Vector& s0,
		  const Dim<2>::Vector& s1) {
  REQUIRE(xmin.x() < xmax.x() and xmin.y() < xmax.y());
  REQUIRE(ncells.size() == 2);
  REQUIRE(ncells[0] > 0 and ncells[1] > 0);

  // Find the min and max bounding mesh planes indicies.
  typedef Dim<2>::Vector Vector;
  const Vector smin = elementWiseMin(s0, s1);
  const Vector smax = elementWiseMax(s0, s1);
  vector<unsigned> ixmin, ixmax;
  for (size_t idim = 0; idim != 2; ++idim) {
    ixmin.push_back(rightCellBoundary(smin(idim), xmin(idim), xmax(idim), ncells[idim]));
    ixmax.push_back( leftCellBoundary(smax(idim), xmin(idim), xmax(idim), ncells[idim]));
    CHECK(ixmin.back() <= ncells[idim]);
    CHECK(ixmax.back() <= ncells[idim]);
  }
  CHECK(ixmin.size() == 2);
  CHECK(ixmax.size() == 2);

  // The intersections are just the intermediate planes.
  vector<Vector> result;
  const double xstep = (xmax.x() - xmin.x())/ncells[0];
  const double ystep = (xmax.y() - xmin.y())/ncells[1];
  for (unsigned iplane = ixmin[0]; iplane < ixmax[0]; ++iplane) {
    const double xseg = xmin.x() + iplane*xstep;
    const Vector meshSeg0(xseg, -1.0e10);
    const Vector meshSeg1(xseg,  1.0e10);
    Vector intersect1, intersect2;
    const char test = segmentSegmentIntersection(s0, s1, meshSeg0, meshSeg1, intersect1, intersect2);
    CONTRACT_VAR(test);
    CHECK(test != '0');
    result.push_back(intersect1);
  }
  for (unsigned iplane = ixmin[1]; iplane < ixmax[1]; ++iplane) {
    const double yseg = xmin.y() + iplane*ystep;
    const Vector meshSeg0(-1.0e10, yseg);
    const Vector meshSeg1( 1.0e10, yseg);
    Vector intersect1, intersect2;
    const char test = segmentSegmentIntersection(s0, s1, meshSeg0, meshSeg1, intersect1, intersect2);
    CONTRACT_VAR(test);
    CHECK(test != '0');
    result.push_back(intersect1);
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    const Vector shat = (s1 - s0).unitVector();
    CONTRACT_VAR(shat);
    const double segLen = (s1 - s0).magnitude();
    for (vector<Vector>::const_iterator itr = result.begin();
	 itr != result.end();
	 ++itr) {
      CONTRACT_VAR(segLen);
      ENSURE((*itr - s0).dot(shat) >= 0.0);
      ENSURE(fuzzyEqual(abs((*itr - s0).dot(shat)), (*itr - s0).magnitude(), 1.0e-10));
      ENSURE((*itr - s0).dot(shat)*safeInv(segLen) <= (1.0 + 1.0e-8));
    }
  }
  END_CONTRACT_SCOPE

  return result;
}
#endif
  
//------------------------------------------------------------------------------
// 3-D.
#ifdef SPHERAL3D
vector<Dim<3>::Vector> 
findIntersections(const Dim<3>::Vector& xmin,
		  const Dim<3>::Vector& xmax,
		  const vector<unsigned>& ncells,
		  const Dim<3>::Vector& s0,
		  const Dim<3>::Vector& s1) {
  REQUIRE(xmin.x() < xmax.x() and xmin.y() < xmax.y() and xmin.z() < xmax.z());
  REQUIRE(ncells.size() == 3);
  REQUIRE(ncells[0] > 0 and ncells[1] > 0 and ncells[2] > 0);

  // Find the min and max bounding mesh planes indicies.
  typedef Dim<3>::Vector Vector;
  const Vector smin = elementWiseMin(s0, s1);
  const Vector smax = elementWiseMax(s0, s1);
  vector<unsigned> ixmin, ixmax;
  for (size_t idim = 0; idim != 3; ++idim) {
    ixmin.push_back(rightCellBoundary(smin(idim), xmin(idim), xmax(idim), ncells[idim]));
    ixmax.push_back( leftCellBoundary(smax(idim), xmin(idim), xmax(idim), ncells[idim]));
    CHECK(ixmin.back() <= ncells[idim]);
    CHECK(ixmax.back() <= ncells[idim]);
  }
  CHECK(ixmin.size() == 3);
  CHECK(ixmax.size() == 3);

  // The intersections are just the intermediate planes.
  vector<Vector> result;
  const double xstep = (xmax.x() - xmin.x())/ncells[0];
  const double ystep = (xmax.y() - xmin.y())/ncells[1];
  const double zstep = (xmax.z() - xmin.z())/ncells[2];
  for (unsigned iplane = ixmin[0]; iplane < ixmax[0]; ++iplane) {
    const double xplane = xmin.x() + iplane*xstep;
    const Vector meshPlane0(xplane, 0.0, 0.0);
    const Vector meshPlane1(xplane, 0.0, 1.0);
    const Vector meshPlane2(xplane, 1.0, 1.0);
    Vector intersect;
    const char test = segmentPlaneIntersection(s0, s1, meshPlane0, meshPlane1, meshPlane2, intersect);
    CONTRACT_VAR(test);
    CHECK(test != '0');
    result.push_back(intersect);
  }
  for (unsigned iplane = ixmin[1]; iplane < ixmax[1]; ++iplane) {
    const double yplane = xmin.y() + iplane*ystep;
    const Vector meshPlane0(0.0, yplane, 0.0);
    const Vector meshPlane1(0.0, yplane, 1.0);
    const Vector meshPlane2(1.0, yplane, 1.0);
    Vector intersect;
    const char test = segmentPlaneIntersection(s0, s1, meshPlane0, meshPlane1, meshPlane2, intersect);
    CONTRACT_VAR(test);
    CHECK(test != '0');
    result.push_back(intersect);
  }
  for (unsigned iplane = ixmin[2]; iplane < ixmax[2]; ++iplane) {
    const double zplane = xmin.z() + iplane*zstep;
    const Vector meshPlane0(0.0, 0.0, zplane);
    const Vector meshPlane1(0.0, 1.0, zplane);
    const Vector meshPlane2(1.0, 1.0, zplane);
    Vector intersect;
    const char test = segmentPlaneIntersection(s0, s1, meshPlane0, meshPlane1, meshPlane2, intersect);
    CONTRACT_VAR(test);
    CHECK(test != '0');
    result.push_back(intersect);
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    const Vector shat = (s1 - s0).unitVector();
    CONTRACT_VAR(shat);
    const double segLen = (s1 - s0).magnitude();
    for (vector<Vector>::const_iterator itr = result.begin();
	 itr != result.end();
	 ++itr) {
      CONTRACT_VAR(segLen);
      ENSURE((*itr - s0).dot(shat) >= 0.0);
      ENSURE(fuzzyEqual(abs((*itr - s0).dot(shat)), (*itr - s0).magnitude(), 1.0e-10));
      ENSURE((*itr - s0).dot(shat)*safeInv(segLen) <= 1.0);
    }
  }
  END_CONTRACT_SCOPE

  return result;
}
#endif

//------------------------------------------------------------------------------
// Find the finest non-zero value in the level set of values at the give point.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
Value
finestNonZeroValue(const vector<vector<Value> >& values,
                   const typename Dimension::Vector& xmin,
                   const typename Dimension::Vector& xmax,
                   const vector<unsigned>& ncells,
                   const typename Dimension::Vector& point) {
  int level = -1;
  Value result = DataTypeTraits<Value>::zero();
  vector<unsigned> ncellsLevel(Dimension::nDim);
  for (unsigned idim = 0; idim != Dimension::nDim; ++idim) ncellsLevel[idim] = 2*ncells[idim];
  while ((result == DataTypeTraits<Value>::zero()) and level < int(values.size() - 1)) {
    ++level;
    for (unsigned idim = 0; idim != Dimension::nDim; ++idim) ncellsLevel[idim] /= 2;
    const size_t index = latticeIndex(point, xmin, xmax, ncellsLevel);
    CHECK(index < values[level].size());
    result = values[level][index];
  }
  return result;
}

//------------------------------------------------------------------------------
// A helpful functor for sorting points according to their distance from a 
// given point.
//------------------------------------------------------------------------------
template<typename Vector>
struct DistanceFromPoint {
  DistanceFromPoint(const Vector& point1, const Vector& point2): 
    mPoint(point1),
    mDelta(point2 - point1) {}
  bool operator()(const Vector& lhs, const Vector& rhs) const {
    return (lhs - mPoint).dot(mDelta) < (rhs - mPoint).dot(mDelta);
  }
  Vector mPoint, mDelta;
};

//------------------------------------------------------------------------------
// integrateThroughMeshAlongSegment
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
Value
integrateThroughMeshAlongSegment(const vector<vector<Value> >& values,
                                 const typename Dimension::Vector& xmin,
                                 const typename Dimension::Vector& xmax,
                                 const vector<unsigned>& ncells,
                                 const typename Dimension::Vector& s0,
                                 const typename Dimension::Vector& s1) {

  typedef typename Dimension::Vector Vector;

  // Preconditions.
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(ncells.size() == Dimension::nDim);
    for (unsigned level = 0; level != values.size(); ++level) {
      unsigned ncellsTotal = 1;
      CONTRACT_VAR(ncellsTotal);
      for (int i = 0; i != Dimension::nDim; ++i) ncellsTotal *= ncells[i]/(1U << level);
      REQUIRE(values[level].size() == ncellsTotal);
    }
  }
  END_CONTRACT_SCOPE

  // Find the points of intersection with the cartesian planes.
  vector<Vector> intersections = findIntersections(xmin, xmax, ncells, s0, s1);

  // Sort the intersection points in order along the line from s0 -> s1.
  sort(intersections.begin(), intersections.end(), DistanceFromPoint<Vector>(s0, s1));

  // Iterate through the intersection points, the interval between each of which 
  // represents a path segment through a cell.
  Value result = DataTypeTraits<Value>::zero();
  Vector lastPoint = s0;
  double cumulativeLength = 0.0;
  CONTRACT_VAR(cumulativeLength);
  for (typename vector<Vector>::const_iterator itr = intersections.begin();
       itr != intersections.end();
       ++itr) {
    const Vector point = 0.5*(lastPoint + *itr);
    const double dl = (*itr - lastPoint).magnitude();
    result += dl*finestNonZeroValue<Dimension, Value>(values, xmin, xmax, ncells, point);
    cumulativeLength += dl;
    lastPoint = *itr;
  }

  // Add the last bit from the last intersection to the end point.
  const Vector point = 0.5*(lastPoint + s1);
  const double dl = (s1 - lastPoint).magnitude();
  cumulativeLength += dl;
  result += dl*finestNonZeroValue<Dimension, Value>(values, xmin, xmax, ncells, point);

  // That's it.
  ENSURE(fuzzyEqual(cumulativeLength, (s1 - s0).magnitude(), 1.0e-10));
  return result;
}

}

