//---------------------------------Spheral++----------------------------------//
// GridCellPlane -- define an integer version of a plane in terms of 
// GridCellIndices.
//
// Created by:  JMO, Wed Mar  8 21:43:23 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_GridCellPlane_hh__
#define __Spheral_GridCellPlane_hh__

namespace Spheral {

template<typename Dimension> class GridCellIndex;

template<typename Dimension>
class GridCellPlane {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  GridCellPlane();
  GridCellPlane(const GridCellPlane& rhs);
  GridCellPlane(const GridCellIndex<Dimension>& point,
                const GridCellIndex<Dimension>& normal);

  // Destructor.
  ~GridCellPlane();

  // Assignment.
  GridCellPlane& operator=(const GridCellPlane& rhs);

  // Access the point and normal which define the plane.
  const GridCellIndex<Dimension>& point() const;
  void setPoint(const GridCellIndex<Dimension>& point);

  const GridCellIndex<Dimension>& normal() const;
  void setNormal(const GridCellIndex<Dimension>& normal);

  // Calculate the minimum distance between a given point and the plane.
  double minimumDistance(const GridCellIndex<Dimension>& point) const;

  // Test if the given point is in the plane.
  bool coplanar(const GridCellIndex<Dimension>& point) const;

  // Various tests which can be applied between planes.
  bool parallel(const GridCellPlane& rhs) const;
  bool operator==(const GridCellPlane& rhs) const;
  bool operator!=(const GridCellPlane& rhs) const;

  // Some tests which can applied between a plane and a point.
  // These are only meant as tests for whether a point lies "above" or "below" the
  // plane, as represented by the direction of the normal.
  bool operator>(const GridCellIndex<Dimension>& point) const;
  bool operator<(const GridCellIndex<Dimension>& point) const;
  bool operator>=(const GridCellIndex<Dimension>& point) const;
  bool operator<=(const GridCellIndex<Dimension>& point) const;

  // The usual validity test.
  bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  GridCellIndex<Dimension> mPoint;
  GridCellIndex<Dimension> mNormal;
};

}

#include "GridCellPlaneInline.hh"

#endif
