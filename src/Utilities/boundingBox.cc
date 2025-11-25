//---------------------------------Spheral++----------------------------------//
// boundingBox
//
// Compute the minimum bounding box for a set of points.
//
// Created by JMO, Wed Oct 19 09:58:58 PDT 2011
//----------------------------------------------------------------------------//
#include "boundingBox.hh"
#include "Geometry/Dimension.hh"
#include "Field/FieldList.hh"

#include <vector>
#include <algorithm>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the minimum volume box containing all the points in a vector of
// positions.
// Note no global operations here!
//------------------------------------------------------------------------------
template<typename Vector>
void
boundingBox(const vector<Vector>& positions,
            Vector& xmin,
            Vector& xmax) {

  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  const unsigned n = positions.size();
  for (unsigned i = 0; i != n; ++i) {
    xmin = elementWiseMin(xmin, positions[i]);
    xmax = elementWiseMax(xmax, positions[i]);
  }
}

//------------------------------------------------------------------------------
// Compute the minimum volume box containing all the points in a FieldList of
// positions.
// Note no global operations here!
//------------------------------------------------------------------------------
template<typename Dimension>
void
boundingBox(const FieldList<Dimension, typename Dimension::Vector>& positions,
            typename Dimension::Vector& xmin,
            typename Dimension::Vector& xmax,
            const bool useGhosts) {

  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  const unsigned numFields = positions.numFields();
  for (unsigned ifield = 0; ifield != numFields; ++ifield) {
    const unsigned n = (useGhosts ? positions[ifield]->numElements() : positions[ifield]->numInternalElements());
    for (unsigned i = 0; i != n; ++i) {
      xmin = elementWiseMin(xmin, positions(ifield, i));
      xmax = elementWiseMax(xmax, positions(ifield, i));
    }
  }
}

}
