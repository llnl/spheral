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

