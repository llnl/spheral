//------------------------------------------------------------------------------
// A collection of little standalone functions that we use as part of building
// the meshes.
//------------------------------------------------------------------------------
#ifndef __Spheral_MeshConstructionUtilities__
#define __Spheral_MeshConstructionUtilities__

#include "Utilities/boundPointWithinBox.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/packElement.hh"
#include "Distributed/Communicator.hh"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdint.h>
#include <iostream>
#include <vector>
#include <set>
#include <tuple>

namespace Spheral {

//------------------------------------------------------------------------------
// A trait class to hold some specializations per mesh dimension.
//------------------------------------------------------------------------------
template<typename Vector>
struct MeshTraits {
  static const uint64_t ncells, ncoarse, ncellsperbin;
};

//------------------------------------------------------------------------------
// Return the quantized cell size in each dimension.
// It is critical that ncells be the same in hashPosition and quantizedPosition!
//------------------------------------------------------------------------------
template<typename Vector>
inline
Vector
cellSize(const Vector& xmin,
         const Vector& xmax) {
  return (xmax - xmin)/MeshTraits<Vector>::ncells;
}

//------------------------------------------------------------------------------
// hashPosition
// It is critical that ncells be the same in hashPosition and quantizedPosition!
//------------------------------------------------------------------------------
template<typename Int>
inline
Int
hashX(const double x,
      const double x0, 
      const double x1,
      const double boxInv,
      const Int ncells,
      const Int ncoarse,
      const Int ncellsperbin) {
  CONTRACT_VAR(x1);
  REQUIRE2((uint64_t(ncells)) % ncoarse == 0U, ncells << " " << ncoarse << " : " << (ncells % ncoarse));
  REQUIRE2(ncellsperbin == (uint64_t(ncells))/ncoarse, ncellsperbin << " " << ncells/ncoarse);
  double f = std::max(0.0, std::min(1.0 - std::numeric_limits<double>::epsilon(), 
                                    std::max(0.0, std::min(1.0, (x - x0)*boxInv))));
  const unsigned icoarse = unsigned(f*ncoarse);
  f = std::max(0.0, f - double(icoarse)/ncoarse);
  return icoarse * ncellsperbin + Int(f * ncells);
}

inline
std::tuple<uint64_t, uint64_t, uint64_t>
hashPosition(const Dim<1>::Vector& position,
             const Dim<1>::Vector& xmin,
             const Dim<1>::Vector& xmax,
             const Dim<1>::Vector& boxInv) {
  typedef Dim<1>::Vector Vector;
  const uint64_t ix0 = hashX(position.x(), xmin.x(), xmax.x(), boxInv.x(), MeshTraits<Vector>::ncells, MeshTraits<Vector>::ncoarse, MeshTraits<Vector>::ncellsperbin);
  return std::make_tuple(ix0, 0U, 0U);
}

inline
std::tuple<uint64_t, uint64_t, uint64_t>
hashPosition(const Dim<2>::Vector& position,
             const Dim<2>::Vector& xmin,
             const Dim<2>::Vector& xmax,
             const Dim<2>::Vector& boxInv) {
  typedef Dim<2>::Vector Vector;
  const uint64_t ix0 = hashX(position.x(), xmin.x(), xmax.x(), boxInv.x(), MeshTraits<Vector>::ncells, MeshTraits<Vector>::ncoarse, MeshTraits<Vector>::ncellsperbin);
  const uint64_t iy0 = hashX(position.y(), xmin.y(), xmax.y(), boxInv.y(), MeshTraits<Vector>::ncells, MeshTraits<Vector>::ncoarse, MeshTraits<Vector>::ncellsperbin);
  return std::make_tuple(ix0, iy0, 0U);
}

inline
std::tuple<uint64_t, uint64_t, uint64_t>
hashPosition(const Dim<3>::Vector& position,
             const Dim<3>::Vector& xmin,
             const Dim<3>::Vector& xmax,
             const Dim<3>::Vector& boxInv) {
  typedef Dim<3>::Vector Vector;
  const uint64_t ix0 = hashX(position.x(), xmin.x(), xmax.x(), boxInv.x(), MeshTraits<Vector>::ncells, MeshTraits<Vector>::ncoarse, MeshTraits<Vector>::ncellsperbin);
  const uint64_t iy0 = hashX(position.y(), xmin.y(), xmax.y(), boxInv.y(), MeshTraits<Vector>::ncells, MeshTraits<Vector>::ncoarse, MeshTraits<Vector>::ncellsperbin);
  const uint64_t iz0 = hashX(position.z(), xmin.z(), xmax.z(), boxInv.z(), MeshTraits<Vector>::ncells, MeshTraits<Vector>::ncoarse, MeshTraits<Vector>::ncellsperbin);
  return std::make_tuple(ix0, iy0, iz0);
}

//------------------------------------------------------------------------------
// Return the center of cell position for the given tuple hash.
// It is critical that ncells be the same in hashPosition and quantizedPosition!
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
quantizedPosition(const std::tuple<uint64_t, uint64_t, uint64_t>& hash,
                  const Dim<1>::Vector& xmin,
                  const Dim<1>::Vector& xmax) {
  const Dim<1>::Vector dx = cellSize(xmin, xmax);
  return boundPointWithinBox(Dim<1>::Vector(xmin.x() + std::get<0>(hash)*dx.x()) +
                             0.5*dx, xmin, xmax);
}

inline
Dim<2>::Vector
quantizedPosition(const std::tuple<uint64_t, uint64_t, uint64_t>& hash,
                  const Dim<2>::Vector& xmin,
                  const Dim<2>::Vector& xmax) {
  const Dim<2>::Vector dx = cellSize(xmin, xmax);
  return boundPointWithinBox(Dim<2>::Vector(xmin.x() + std::get<0>(hash)*dx.x(),
                                            xmin.y() + std::get<1>(hash)*dx.y()) +
                             0.5*dx, xmin, xmax);
}

inline
Dim<3>::Vector
quantizedPosition(const std::tuple<uint64_t, uint64_t, uint64_t>& hash,
                  const Dim<3>::Vector& xmin,
                  const Dim<3>::Vector& xmax) {
  const Dim<3>::Vector dx = cellSize(xmin, xmax);
  return boundPointWithinBox(Dim<3>::Vector(xmin.x() + std::get<0>(hash)*dx.x(),
                                            xmin.y() + std::get<1>(hash)*dx.y(),
                                            xmin.z() + std::get<2>(hash)*dx.z()) +
                             0.5*dx, xmin, xmax);
}

//------------------------------------------------------------------------------
// Specialized comparison for a pair of hashed positions to see if they're equal.
//------------------------------------------------------------------------------
inline
int
compare(const std::tuple<uint64_t, uint64_t, uint64_t>& lhs,
        const std::tuple<uint64_t, uint64_t, uint64_t>& rhs) {
  using namespace std;
  return (get<2>(lhs) < get<2>(rhs)                     ? -1 :
          get<2>(lhs) > get<2>(rhs)                     ?  1 :
          get<1>(lhs) < get<1>(rhs)                     ? -1 :
          get<1>(lhs) > get<1>(rhs)                     ?  1 :
          get<0>(lhs) < get<0>(rhs)                     ? -1 : 
          get<0>(lhs) > get<0>(rhs)                     ?  1 : 
                                                           0);
}

inline
int
fuzzyCompare(const std::tuple<uint64_t, uint64_t, uint64_t>& lhs,
             const std::tuple<uint64_t, uint64_t, uint64_t>& rhs,
             const uint64_t delta) {
  const unsigned dx = std::max(std::get<0>(lhs), std::get<0>(rhs)) - std::min(std::get<0>(lhs), std::get<0>(rhs));
  const unsigned dy = std::max(std::get<1>(lhs), std::get<1>(rhs)) - std::min(std::get<1>(lhs), std::get<1>(rhs));
  const unsigned dz = std::max(std::get<2>(lhs), std::get<2>(rhs)) - std::min(std::get<2>(lhs), std::get<2>(rhs));
  return ((dx <= delta and dy <= delta and dz <= delta) ?  0 :
          std::get<2>(lhs) < std::get<2>(rhs)                     ? -1 :
          std::get<2>(lhs) > std::get<2>(rhs)                     ?  1 :
          std::get<1>(lhs) < std::get<1>(rhs)                     ? -1 :
          std::get<1>(lhs) > std::get<1>(rhs)                     ?  1 :
          std::get<0>(lhs) < std::get<0>(rhs)                     ? -1 : 
                                                           1);
}

//------------------------------------------------------------------------------
// Add two keys.
//------------------------------------------------------------------------------
template<typename T>
inline
void
incrementTuple(std::tuple<T, T, T>& lhs,
               const std::tuple<T, T, T>& rhs) {
  std::get<0>(lhs) += std::get<0>(rhs);
  std::get<1>(lhs) += std::get<1>(rhs);
  std::get<2>(lhs) += std::get<2>(rhs);
}

//------------------------------------------------------------------------------
// Divide a key.
//------------------------------------------------------------------------------
template<typename T>
inline
void
divideTuple(std::tuple<T, T, T>& lhs, const T& rhs) {
  std::get<0>(lhs) /= rhs;
  std::get<1>(lhs) /= rhs;
  std::get<2>(lhs) /= rhs;
}

//------------------------------------------------------------------------------
// Specialized comparison for a pair of hashed edges to see if they're equal.
//------------------------------------------------------------------------------
inline
int
compare(const std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>& lhs,
        const std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>& rhs) {
  using namespace std;
  return (get<5>(lhs) < get<5>(rhs)                     ? -1 :
          get<5>(lhs) > get<5>(rhs)                     ?  1 :
          get<4>(lhs) < get<4>(rhs)                     ? -1 :
          get<4>(lhs) > get<4>(rhs)                     ?  1 :
          get<3>(lhs) < get<3>(rhs)                     ? -1 : 
          get<3>(lhs) > get<3>(rhs)                     ?  1 : 
          get<2>(lhs) < get<2>(rhs)                     ? -1 :
          get<2>(lhs) > get<2>(rhs)                     ?  1 :
          get<1>(lhs) < get<1>(rhs)                     ? -1 :
          get<1>(lhs) > get<1>(rhs)                     ?  1 :
          get<0>(lhs) < get<0>(rhs)                     ? -1 : 
          get<0>(lhs) > get<0>(rhs)                     ?  1 : 
                                                           0);
}

//------------------------------------------------------------------------------
// Take a set of vertices and return their Key values, collapsing degenerate 
// vertices into single hashed values.
//------------------------------------------------------------------------------
template<typename Vector, typename Uint>
inline
std::vector<std::tuple<Uint, Uint, Uint> >
collapseDegenerateVertices(const std::vector<Vector>& vertices,
                           const Vector& xmin,
                           const Vector& xmax,
                           const Vector& boxInv,
                           const Uint tol) {
  typedef std::tuple<Uint, Uint, Uint> Key;
  using std::vector;
  using std::set;

  unsigned i, j;
  const unsigned n = vertices.size();
  Key hashi;

  // Make a first pass through all the points, converting them to their integer
  // hashed positions.
  vector<Key> result;
  result.reserve(n);
  for (i = 0; i != n; ++i) {
    result.push_back(hashPosition(vertices[i], xmin, xmax, boxInv));
  }
  CHECK(result.size() == n);

  // Look for who is linked to whom, friends of friends style.
  if (tol > 0U) {
    Key hashi;
    vector<unsigned> assigned(n, 0U);
    vector<set<unsigned> > neighborsForPoint(n);
    for (i = 0; i <= (n - 1); ++i) {
      neighborsForPoint[i].insert(i);
      for (j = i + 1; j != n; ++j) {
        if (fuzzyCompare(result[j], result[i], tol) == 0) {
          neighborsForPoint[i].insert(j);
          neighborsForPoint[j].insert(i);
        }
      }
    }
    for (i = 0; i != n; ++i) {
      if (assigned[i] == 0U) {
        set<unsigned> fullNeighbors;
        for (typename set<unsigned>::const_iterator itr = neighborsForPoint[i].begin();
             itr != neighborsForPoint[i].end();
             ++itr) {
          for (typename set<unsigned>::const_iterator otherItr = neighborsForPoint[*itr].begin();
               otherItr != neighborsForPoint[*itr].end();
               ++otherItr) {
            fullNeighbors.insert(*otherItr);
          }
        }
        CHECK(fullNeighbors.size() > 0U);
        hashi = Key(0U, 0U, 0U);
        for (typename set<unsigned>::iterator itr = fullNeighbors.begin();
             itr != fullNeighbors.end();
             ++itr) {
          incrementTuple(hashi, result[*itr]);
        }
        divideTuple<Uint>(hashi, fullNeighbors.size());
        for (typename set<unsigned>::iterator itr = fullNeighbors.begin();
             itr != fullNeighbors.end();
             ++itr) {
          result[*itr] = hashi;
          assigned[*itr] = 1U;
        }
      }
    }
  }

  // That's it.
  ENSURE(result.size() == n);
  return result;
}

//------------------------------------------------------------------------------
// Convert a Voro++ vertex position to lab coordinates.
//------------------------------------------------------------------------------
template<typename Vector, typename Float>
inline
Vector
voroVertex(const Vector& centroid,
           const Vector& xmin,
           const Vector& xmax,
           const Float& ptx,
           const Float& pty,
           const Float& ptz) {
  return boundPointWithinBox(centroid + Vector(0.5*ptx, 0.5*pty, 0.5*ptz), xmin, xmax);
}

//------------------------------------------------------------------------------
// Exchange a set of tuple keys.
//------------------------------------------------------------------------------
template<typename T>
inline
void
exchangeTuples(const std::vector<std::tuple<T, T, T> >& localKeys,
               const std::vector<unsigned>& neighborDomains,
               std::vector<std::vector<std::tuple<T, T, T> > >& neighborKeys) {
  CONTRACT_VAR(localKeys);
  CONTRACT_VAR(neighborDomains);
  CONTRACT_VAR(neighborKeys);
#ifdef USE_MPI
  typedef std::tuple<T, T, T> Key;
  using std::vector;

  // Pack up our local keys.
  std::vector<char> localPacked;
  packElement(localKeys, localPacked);
  unsigned localBufferSize = localPacked.size();

  // Exchange the sizes of the hashed keys.
  const unsigned numNeighborDomains = neighborDomains.size();
  if (numNeighborDomains > 0) {
    std::vector<unsigned> neighborBufferSizes(numNeighborDomains);
    {
      std::vector<MPI_Request> requests(2*numNeighborDomains);
      for (unsigned k = 0; k != numNeighborDomains; ++k) {
        const unsigned otherProc = neighborDomains[k];
        MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &requests[k]);
        MPI_Irecv(&neighborBufferSizes[k], 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &requests[numNeighborDomains + k]);
      }
      std::vector<MPI_Status> status(requests.size());
      MPI_Waitall(2*numNeighborDomains, &requests.front(), &status.front());
    }

    // Now exchange the hashed keys themselves.
    std::vector<vector<char> > buffers;
    for (unsigned k = 0; k != numNeighborDomains; ++k) buffers.push_back(std::vector<char>(neighborBufferSizes[k]));
    {
      std::vector<MPI_Request> requests(2*numNeighborDomains);
      for (unsigned k = 0; k != numNeighborDomains; ++k) {
        const unsigned otherProc = neighborDomains[k];
        MPI_Isend(&localPacked.front(), localBufferSize, MPI_CHAR, otherProc, 2, Communicator::communicator(), &requests[k]);
        CHECK(buffers[k].size() == neighborBufferSizes[k]);
        MPI_Irecv(&buffers[k].front(), neighborBufferSizes[k], MPI_CHAR, otherProc, 2, Communicator::communicator(), &requests[numNeighborDomains + k]);
      }
      std::vector<MPI_Status> status(requests.size());
      MPI_Waitall(2*numNeighborDomains, &requests.front(), &status.front());
    }

    // Unpack the neighbor hashes.
    neighborKeys = std::vector<std::vector<Key> >(numNeighborDomains);
    for (unsigned k = 0; k != numNeighborDomains; ++k) {
      std::vector<char>::const_iterator bufItr = buffers[k].begin();
      unpackElement(neighborKeys[k], bufItr, buffers[k].end());
      CHECK(bufItr == buffers[k].end());
    }
  }
#endif
}

//------------------------------------------------------------------------------
// Exchange a set of tuple keys, in this case only a subset.
//------------------------------------------------------------------------------
template<typename T>
inline
void
exchangeTuples(const std::vector<std::tuple<T, T, T> >& localKeys,
               const std::vector<unsigned>& neighborDomains,
               const std::vector<std::vector<unsigned> >& sendIndices,
               std::vector<std::vector<std::tuple<T, T, T> > >& neighborKeys) {
  CONTRACT_VAR(localKeys);
  CONTRACT_VAR(neighborDomains);
  CONTRACT_VAR(sendIndices);
  CONTRACT_VAR(neighborKeys);
#ifdef USE_MPI
  typedef std::tuple<T, T, T> Key;

  const unsigned numNeighborDomains = neighborDomains.size();
  REQUIRE(sendIndices.size() == numNeighborDomains);
  if (numNeighborDomains > 0) {

    // Pack up our local keys.
    std::vector<std::vector<char> > localPacked(numNeighborDomains);
    std::vector<unsigned> localBufferSizes;
    for (unsigned k = 0; k != numNeighborDomains; ++k) {
      std::vector<Key> thpt;
      for (unsigned i = 0; i != sendIndices[k].size(); ++i) thpt.push_back(localKeys[sendIndices[k][i]]);
      packElement(thpt, localPacked[k]);
      localBufferSizes.push_back(localPacked[k].size());
    }
    CHECK(localPacked.size() == numNeighborDomains);
    CHECK(localBufferSizes.size() == numNeighborDomains);

    // Exchange the sizes of the hashed keys.
    std::vector<unsigned> neighborBufferSizes(numNeighborDomains);
    {
      std::vector<MPI_Request> requests(2*numNeighborDomains);
      for (unsigned k = 0; k != numNeighborDomains; ++k) {
        const unsigned otherProc = neighborDomains[k];
        MPI_Isend(&localBufferSizes[k], 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &requests[k]);
        MPI_Irecv(&neighborBufferSizes[k], 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &requests[numNeighborDomains + k]);
      }
      std::vector<MPI_Status> status(requests.size());
      MPI_Waitall(2*numNeighborDomains, &requests.front(), &status.front());
    }

    // Now exchange the hashed keys themselves.
    std::vector<std::vector<char> > buffers;
    for (unsigned k = 0; k != numNeighborDomains; ++k) buffers.push_back(std::vector<char>(neighborBufferSizes[k]));
    {
      std::vector<MPI_Request> requests(2*numNeighborDomains);
      for (unsigned k = 0; k != numNeighborDomains; ++k) {
        const unsigned otherProc = neighborDomains[k];
        MPI_Isend(&localPacked[k].front(), localBufferSizes[k], MPI_CHAR, otherProc, 2, Communicator::communicator(), &requests[k]);
        CHECK(buffers[k].size() == neighborBufferSizes[k]);
        MPI_Irecv(&buffers[k].front(), neighborBufferSizes[k], MPI_CHAR, otherProc, 2, Communicator::communicator(), &requests[numNeighborDomains + k]);
      }
      std::vector<MPI_Status> status(requests.size());
      MPI_Waitall(2*numNeighborDomains, &requests.front(), &status.front());
    }

    // Unpack the neighbor hashes.
    neighborKeys = std::vector<std::vector<Key> >(numNeighborDomains);
    for (unsigned k = 0; k != numNeighborDomains; ++k) {
      std::vector<char>::const_iterator bufItr = buffers[k].begin();
      unpackElement(neighborKeys[k], bufItr, buffers[k].end());
      CHECK(bufItr == buffers[k].end());
    }
  }
#endif
}

//------------------------------------------------------------------------------
// Append the result of intersecting two line segments to the given list of 
// points.
//------------------------------------------------------------------------------
inline
void
appendIntersections(const Dim<2>::Vector& a0,
                    const Dim<2>::Vector& a1,
                    const Dim<2>::Vector& b0,
                    const Dim<2>::Vector& b1,
                    std::vector<Dim<2>::Vector>& result) {
  typedef Dim<2>::Vector Vector;
  Vector intersect1, intersect2;
  const char code = segmentSegmentIntersection(a0, a1, b0, b1, intersect1, intersect2);
  if (code == '1' or code == 'v') {
    // Intersect at one point.
    result.push_back(intersect1);
  } else if (code == 'e') {
    // The edge is colinear and overlapping with the boundary.
    result.push_back(intersect1);
    result.push_back(intersect2);
  }
}

//------------------------------------------------------------------------------
// Compute the set of points that results from bounding the given set within
// the given box (defined by it's corner points).
//------------------------------------------------------------------------------
inline
Dim<2>::ConvexHull
boundPolygonInBox(const Dim<2>::ConvexHull& polygon,
                  const std::vector<Dim<2>::Vector>& boundPoints,
                  const double xtol2) {
  CONTRACT_VAR(xtol2);
  REQUIRE(boundPoints.size() == 4);

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::ConvexHull ConvexHull;

  // Generate a set of flags indicating whether points are inside or outside
  // of the box.
  const std::vector<Vector> points = polygon.vertices();
  std::vector<unsigned> flags(points.size(), 0);
  for (size_t i = 0; i != points.size(); ++i) flags[i] = (testPointInBox(points[i], boundPoints[0], boundPoints[2]) ? 1 : 0);
//   CHECK(*max_element(flags.begin(), flags.end()) == 1);

  // If the polygon is entirely within the boundary there's nothing to do.
  if (*min_element(flags.begin(), flags.end()) == 1) return polygon;

  // The polygon straddles the box boundary, so generate new points on the 
  // polygon where it intersects the boundary and remove any points outside.
  std::vector<Vector> newPoints;
  for (size_t i = 0; i != points.size(); ++i) {
    if (flags[i] == 1) newPoints.push_back(points[i]);

    // Check if this edge intersects the boundary.
    const size_t j = (i + 1) % points.size();
    appendIntersections(points[i], points[j], boundPoints[0], boundPoints[1], newPoints);
    appendIntersections(points[i], points[j], boundPoints[1], boundPoints[2], newPoints);
    appendIntersections(points[i], points[j], boundPoints[2], boundPoints[3], newPoints);
    appendIntersections(points[i], points[j], boundPoints[3], boundPoints[0], newPoints);
  }

  // If the polygon contains any of the boundary corners, add them.
  for (unsigned i = 0; i != boundPoints.size(); ++i) {
    if (polygon.contains(boundPoints[i])) newPoints.push_back(boundPoints[i]);
  }
  CHECK(newPoints.size() >= 3);

  // Make sure all points are contained in the boundaries.
  for (unsigned i = 0; i != newPoints.size(); ++i) newPoints[i] = boundPointWithinBox(newPoints[i], boundPoints[0], boundPoints[2]);

  // That's it.
  return ConvexHull(newPoints);
}

//------------------------------------------------------------------------------
// A trait class to help us specialize constructing a Vector2d from Vector3d.
//------------------------------------------------------------------------------
template<unsigned nDim> struct Vector3dto2d;

// Indexing by 0.
template<> struct Vector3dto2d<0> {
  static Dim<2>::Vector extract(const Dim<3>::Vector& v) { return Dim<2>::Vector(v.y(), v.z()); }
  static Dim<3>::Vector build(const Dim<2>::Vector& v2d,
                              const Dim<3>::Vector& v3d) { return Dim<3>::Vector(v3d.x(), v2d.x(), v2d.y()); }
  static Dim<2>::Vector xmin(const Dim<3>::Vector& v1,
                             const Dim<3>::Vector& v2,
                             const Dim<3>::Vector& v3,
                             const Dim<3>::Vector& v4) {
    return Dim<2>::Vector(std::min(v1.y(), std::min(v2.y(), std::min(v3.y(), v4.y()))),
                          std::min(v1.z(), std::min(v2.z(), std::min(v3.z(), v4.z()))));
  }
  static Dim<2>::Vector xmax(const Dim<3>::Vector& v1,
                             const Dim<3>::Vector& v2,
                             const Dim<3>::Vector& v3,
                             const Dim<3>::Vector& v4) {
    return Dim<2>::Vector(std::max(v1.y(), std::max(v2.y(), std::max(v3.y(), v4.y()))),
                          std::max(v1.z(), std::max(v2.z(), std::max(v3.z(), v4.z()))));
  }
};

// Indexing by 1.
template<> struct Vector3dto2d<1> {
  static Dim<2>::Vector extract(const Dim<3>::Vector& v) { return Dim<2>::Vector(v.x(), v.z()); }
  static Dim<3>::Vector build(const Dim<2>::Vector& v2d,
                              const Dim<3>::Vector& v3d) { return Dim<3>::Vector(v2d.x(), v3d.y(), v2d.y()); }
  static Dim<2>::Vector xmin(const Dim<3>::Vector& v1,
                             const Dim<3>::Vector& v2,
                             const Dim<3>::Vector& v3,
                             const Dim<3>::Vector& v4) {
    return Dim<2>::Vector(std::min(v1.x(), std::min(v2.x(), std::min(v3.x(), v4.x()))),
                          std::min(v1.z(), std::min(v2.z(), std::min(v3.z(), v4.z()))));
  }
  static Dim<2>::Vector xmax(const Dim<3>::Vector& v1,
                             const Dim<3>::Vector& v2,
                             const Dim<3>::Vector& v3,
                             const Dim<3>::Vector& v4) {
    return Dim<2>::Vector(std::max(v1.x(), std::max(v2.x(), std::max(v3.x(), v4.x()))),
                          std::max(v1.z(), std::max(v2.z(), std::max(v3.z(), v4.z()))));
  }
};

// Indexing by 2.
template<> struct Vector3dto2d<2> {
  static Dim<2>::Vector extract(const Dim<3>::Vector& v) { return Dim<2>::Vector(v.x(), v.y()); }
  static Dim<3>::Vector build(const Dim<2>::Vector& v2d,
                              const Dim<3>::Vector& v3d) { return Dim<3>::Vector(v2d.x(), v2d.y(), v3d.z()); }
  static Dim<2>::Vector xmin(const Dim<3>::Vector& v1,
                             const Dim<3>::Vector& v2,
                             const Dim<3>::Vector& v3,
                             const Dim<3>::Vector& v4) {
    return Dim<2>::Vector(std::min(v1.x(), std::min(v2.x(), std::min(v3.x(), v4.x()))),
                          std::min(v1.y(), std::min(v2.y(), std::min(v3.y(), v4.y()))));
  }
  static Dim<2>::Vector xmax(const Dim<3>::Vector& v1,
                             const Dim<3>::Vector& v2,
                             const Dim<3>::Vector& v3,
                             const Dim<3>::Vector& v4) {
    return Dim<2>::Vector(std::max(v1.x(), std::max(v2.x(), std::max(v3.x(), v4.x()))),
                          std::max(v1.y(), std::max(v2.y(), std::max(v3.y(), v4.y()))));
  }
};

//------------------------------------------------------------------------------
// For coplanar points, append the intersections with the planar bounding box.
//------------------------------------------------------------------------------
template<unsigned planarDim>
inline
void
append2dIntersections(const Dim<3>::Vector& a0_3d,
                      const Dim<3>::Vector& a1_3d,
                      const Dim<2>::Vector& a0,
                      const Dim<2>::Vector& a1,
                      const Dim<2>::Vector& b0,
                      const Dim<2>::Vector& b1,
                      std::vector<Dim<3>::Vector>& result) {
  CONTRACT_VAR(a1_3d);
  typedef Dim<2>::Vector Vector2d;
  Vector2d intersect1, intersect2;
  const char code = segmentSegmentIntersection(a0, a1, b0, b1, intersect1, intersect2);
  if (code == '1' or code == 'v') {
    // Intersect at one point.
    result.push_back(Vector3dto2d<planarDim>::build(intersect1, a0_3d));
  } else if (code == 'e') {
    // The edge is colinear and overlapping with the boundary.
    result.push_back(Vector3dto2d<planarDim>::build(intersect1, a0_3d));
    result.push_back(Vector3dto2d<planarDim>::build(intersect2, a0_3d));
  }
}

//------------------------------------------------------------------------------
// Append the result of intersecting a line segment with a square planar 
// section.
//------------------------------------------------------------------------------
template<unsigned planarDim>
inline
void
appendIntersections(const Dim<3>::Vector& a0,
                    const Dim<3>::Vector& a1,
                    const Dim<3>::Vector& p0,
                    const Dim<3>::Vector& p1,
                    const Dim<3>::Vector& p2,
                    const Dim<3>::Vector& p3,
                    std::vector<Dim<3>::Vector>& result) {
  REQUIRE(p1(planarDim) == p0(planarDim) and
          p2(planarDim) == p0(planarDim) and
          p3(planarDim) == p0(planarDim));
  typedef Dim<3>::Vector Vector;
  typedef Dim<2>::Vector Vector2d;
  Vector intersect;
  const char code = segmentPlaneIntersection(a0, a1, p0, p1, p2, intersect);
  CHECK(code != 'd');
  if (code == '0') {
    return;
  } else if (code == '1') {
    const Vector delta = intersect - p1;
    const Vector l1 = p0 - p1;
    const Vector l1hat = l1.unitVector();
    const double t1 = delta.dot(l1hat);
    if (t1 >= 0.0 and t1 <= l1.magnitude()) {
      const Vector l2 = p2 - p1;
      const Vector l2hat = l2.unitVector();
      const double t2 = delta.dot(l2hat);
      if (t2 >= 0.0 and t2 <= (p2 - p1).magnitude()) result.push_back(intersect);
    }
  } else {
    CHECK(code == 'p');
    const Vector2d xmin = Vector3dto2d<planarDim>::xmin(p0, p1, p2, p3);
    const Vector2d xmax = Vector3dto2d<planarDim>::xmax(p0, p1, p2, p3);
    const Vector2d s0 = Vector3dto2d<planarDim>::extract(a0);
    const Vector2d s1 = Vector3dto2d<planarDim>::extract(a1);
    const bool t1 = testPointInBox(s0, xmin, xmax);
    const bool t2 = testPointInBox(s1, xmin, xmax);
    if (t1) result.push_back(a0);
    if (t2) result.push_back(a1);
    if ((not t1) or (not t2)) {
      const Vector2d b0 = Vector3dto2d<planarDim>::extract(p0);
      const Vector2d b1 = Vector3dto2d<planarDim>::extract(p1);
      const Vector2d b2 = Vector3dto2d<planarDim>::extract(p2);
      const Vector2d b3 = Vector3dto2d<planarDim>::extract(p3);
      append2dIntersections<planarDim>(a0, a1, s0, s1, b0, b1, result);
      append2dIntersections<planarDim>(a0, a1, s0, s1, b1, b2, result);
      append2dIntersections<planarDim>(a0, a1, s0, s1, b2, b3, result);
      append2dIntersections<planarDim>(a0, a1, s0, s1, b3, b0, result);
    }
  }
}

//------------------------------------------------------------------------------
// Compute the set of points that results from bounding the given set within
// the given box (defined by it's corner points).
//------------------------------------------------------------------------------
inline
Dim<3>::ConvexHull
boundPolyhedronInBox(const Dim<3>::ConvexHull& polyhedron,
                     const std::vector<Dim<3>::Vector>& boundPoints,
                     const double xtol2) {
  CONTRACT_VAR(xtol2);
  REQUIRE(boundPoints.size() == 8);

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::ConvexHull ConvexHull;

  // Generate a set of flags indicating whether points are inside or outside
  // of the box.
  const std::vector<Vector> points = polyhedron.vertices();
  std::vector<unsigned> flags(points.size(), 0);
  const Vector& xmin = boundPoints[0];
  const Vector& xmax = boundPoints[6];
  for (size_t i = 0; i != points.size(); ++i) flags[i] = (testPointInBox(points[i], xmin, xmax) ? 1 : 0);
//   CHECK(*max_element(flags.begin(), flags.end()) == 1);

  // If the polyhedron is entirely within the boundary there's nothing to do.
  if (*min_element(flags.begin(), flags.end()) == 1) return polyhedron;

  // The polyhedron straddles the box boundary.  First take any vertices of the 
  // polyhedron that are inside the box.
  std::vector<Vector> newPoints;
  for (size_t i = 0; i != points.size(); ++i) {
    if (flags[i] == 1) newPoints.push_back(points[i]);
  }

  // Walk the edge and look for any that intersect the boundary.
  const std::vector<Vector>& vertices = polyhedron.vertices();
  const std::vector<std::pair<unsigned, unsigned> > edges = polyhedron.edges();
  unsigned i, j;
  for (std::vector<std::pair<unsigned, unsigned> >::const_iterator edgeItr = edges.begin();
       edgeItr != edges.end();
       ++edgeItr) {
    i = edgeItr->first;
    j = edgeItr->second;
    if (flags[i] + flags[j] == 1) {
      const Vector& p1 = vertices[i];
      const Vector& p2 = vertices[j];
      appendIntersections<0>(p1, p2, boundPoints[3], boundPoints[0], boundPoints[4], boundPoints[7], newPoints);
      appendIntersections<0>(p1, p2, boundPoints[1], boundPoints[2], boundPoints[6], boundPoints[5], newPoints);
      appendIntersections<1>(p1, p2, boundPoints[2], boundPoints[3], boundPoints[7], boundPoints[6], newPoints);
      appendIntersections<1>(p1, p2, boundPoints[0], boundPoints[1], boundPoints[5], boundPoints[4], newPoints);
      appendIntersections<2>(p1, p2, boundPoints[0], boundPoints[3], boundPoints[2], boundPoints[1], newPoints);
      appendIntersections<2>(p1, p2, boundPoints[5], boundPoints[6], boundPoints[7], boundPoints[4], newPoints);
    }
  }

  // If the polyhedron contains any of the boundary corners, add them.
  for (unsigned k = 0; k != boundPoints.size(); ++k) {
    if (polyhedron.contains(boundPoints[k])) newPoints.push_back(boundPoints[k]);
  }
  CHECK(newPoints.size() > 3);

  // Make sure all points are contained in the boundaries.
  for (unsigned i = 0; i != newPoints.size(); ++i) newPoints[i] = boundPointWithinBox(newPoints[i], xmin, xmax);

  // That's it.
  return ConvexHull(newPoints);
}

//------------------------------------------------------------------------------
// Generic method for adding things to our element hash sets.
//------------------------------------------------------------------------------
template<typename Hash2KeyMap, typename Key2HashMap>
inline
void
extendElementHashSets(const typename Hash2KeyMap::key_type& key,
                      Hash2KeyMap& hash2id,
                      Key2HashMap& id2hash,
                      typename Hash2KeyMap::mapped_type& nextID) {
  if (hash2id.find(key) == hash2id.end()) {
    hash2id[key] = nextID;
    id2hash[nextID] = key;
    ++nextID;
  }
}

//------------------------------------------------------------------------------
// Construct an edge hash as a std::pair<unsigned, unsigned>.
//------------------------------------------------------------------------------
inline
std::pair<unsigned, unsigned>
hashEdge(const unsigned i, const unsigned j) {
  return i < j ? std::make_pair(i, j) : std::make_pair(j, i);
}

template<typename T>
inline
std::tuple<T, T, T, T, T, T>
hashEdge(const std::tuple<T, T, T>& hashi,
         const std::tuple<T, T, T>& hashj) {
  return (hashi < hashj ? 
          std::make_tuple(std::get<0>(hashi), std::get<1>(hashi), std::get<2>(hashi),
                          std::get<0>(hashj), std::get<1>(hashj), std::get<2>(hashj)) :
          std::make_tuple(std::get<0>(hashj), std::get<1>(hashj), std::get<2>(hashj),
                          std::get<0>(hashi), std::get<1>(hashi), std::get<2>(hashi)));
}

}

//------------------------------------------------------------------------------
// Provide comparison and other useful operators for the tuple Key type.
//------------------------------------------------------------------------------
namespace std {
  //----------------------------------------------------------------------------
  // tuple<T, T, T>
  //----------------------------------------------------------------------------
  // operator==
  template<typename T>
  inline
  bool operator==(const ::std::tuple<T, T, T>& lhs,
                  const ::std::tuple<T, T, T>& rhs) {
    return (Spheral::compare(lhs, rhs) == 0);
  }

  // operator!=
  template<typename T>
  inline
  bool operator!=(const ::std::tuple<T, T, T>& lhs,
                  const ::std::tuple<T, T, T>& rhs) {
    return (Spheral::compare(lhs, rhs) != 0);
  }

  // operator<
  template<typename T>
  inline
  bool operator<(const ::std::tuple<T, T, T>& lhs,
                 const ::std::tuple<T, T, T>& rhs) {
    return (Spheral::compare(lhs, rhs) == -1);
  }

  // operator<<
  template<typename T>
  inline
  std::ostream&
  operator<<(std::ostream& os, const ::std::tuple<T, T, T>& x) {
    using namespace std;
    os << "(" << get<0>(x) << " " << get<1>(x) << " " << get<2>(x) << ")";
    return os;
  }

  //----------------------------------------------------------------------------
  // tuple<T, T, T, T, T, T>
  //----------------------------------------------------------------------------
  // operator==
  template<typename T>
  inline
  bool operator==(const ::std::tuple<T, T, T, T, T, T>& lhs,
                  const ::std::tuple<T, T, T, T, T, T>& rhs) {
    return (Spheral::compare(lhs, rhs) == 0);
  }

  // operator!=
  template<typename T>
  inline
  bool operator!=(const ::std::tuple<T, T, T, T, T, T>& lhs,
                  const ::std::tuple<T, T, T, T, T, T>& rhs) {
    return (Spheral::compare(lhs, rhs) != 0);
  }

  // operator<
  template<typename T>
  inline
  bool operator<(const ::std::tuple<T, T, T, T, T, T>& lhs,
                 const ::std::tuple<T, T, T, T, T, T>& rhs) {
    return (Spheral::compare(lhs, rhs) == -1);
  }

  // operator<<
  template<typename T>
  inline
  std::ostream&
  operator<<(std::ostream& os, const ::std::tuple<T, T, T, T, T, T>& x) {
    using namespace std;
    os << "(" 
       << get<0>(x) << " " << get<1>(x) << " " << get<2>(x) 
       << get<3>(x) << " " << get<4>(x) << " " << get<5>(x) 
       << ")";
    return os;
  }
}

#endif
