//---------------------------------Spheral++------------------------------------
// Compute the volume per point using an inverse convex hull.
// Optionally this volume can then be clipped to the Voronoi.
//------------------------------------------------------------------------------
#include "computeHullVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/range.hh"
#include "Geometry/PolyClipperUtilities.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
#include <utility>
#include <limits>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

namespace {  // anonymous
//------------------------------------------------------------------------------
// We put some helper methods for clipping to the Voronoi volume in the
// anonymous namespace.  Mostly lifted from what we do in computeVoronoiVolume.
//
// Trait class for local Dimension
//------------------------------------------------------------------------------
template<typename Dimension> struct ClippingType;

//..............................................................................
// 2D
template<>
struct
ClippingType<Dim<2>> {
  using Scalar = Dim<2>::Scalar;
  using Vector = Dim<2>::Vector;
  using SymTensor = Dim<2>::SymTensor;
  using FacetedVolume = Dim<2>::FacetedVolume;
  using Plane = PolyClipperPlane2d;
  using PolyVolume = PolyClipperPolygon;
  
  // Unit circle as template shape.
  static PolyVolume unitPolyVolume() {
    const auto nverts = 18;
    PolyVolume cell0;
    {
      const auto dtheta = 2.0*M_PI/nverts;
      vector<Vector> verts0(nverts);
      vector<vector<unsigned>> facets0(nverts, vector<unsigned>(2));
      for (auto j = 0; j != nverts; ++j) {
        const auto theta = j*dtheta;
        verts0[j].x(cos(theta));
        verts0[j].y(sin(theta));
        facets0[j][0] = j;
        facets0[j][1] = (j + 1) % nverts;
      }
      convertToPolyClipper(cell0, FacetedVolume(verts0, facets0));
    }
    return cell0;
  }

  // clipping operation
  static void clip(PolyVolume& cell, const std::vector<Plane>& planes) {
    PolyClipper::clipPolygon(cell, planes);
  }

  // moment computation
  static void moments(double& vol, Vector& cent, const PolyVolume& cell) {
    PolyClipper::moments(vol, cent, cell);
  }

  // collapse degenerate points
  static void collapseDegenerates(PolyVolume& cell, const double tol) {
    PolyClipper::collapseDegenerates(cell, tol);
  }

  // Generate the reduced void point stencil -- up to 4 for 2D
  static std::vector<Vector> createEtaVoidPoints(const Vector& etaVoidAvg,
                                                 const int nvoid,
                                                 const double rin,
                                                 const SymTensor& /*Hi*/,
                                                 const SymTensor& /*Hinvi*/,
                                                 const PolyVolume& /*celli*/) {
    std::vector<Vector> result;
    const auto nverts = 18;
    const auto thetaVoidAvg = atan2(etaVoidAvg.y(), etaVoidAvg.x());
    const auto nv = max(1U, min(4U, unsigned(4.0*double(nvoid)/double(nverts))));
    for (unsigned k = 0; k != nv; ++k) {
      const auto theta = thetaVoidAvg + (0.5*k - 0.25*(nv - 1))*M_PI;
      result.push_back(Vector(0.5*rin*cos(theta), 0.5*rin*sin(theta)));
    }
    ENSURE(result.size() == nv);
    return result;
  }
    
  // toString
  static std::string toString(const PolyVolume& celli) {
    return PolyClipper::polygon2string(celli);
  }

};

//..............................................................................
// 3D
template<>
struct
ClippingType<Dim<3>> {
  using Scalar = Dim<3>::Scalar;
  using Vector = Dim<3>::Vector;
  using SymTensor = Dim<3>::SymTensor;
  using FacetedVolume = Dim<3>::FacetedVolume;
  using Plane = PolyClipperPlane3d;
  using PolyVolume = PolyClipperPolyhedron;
  
  // Build an approximation of the starting kernel shape (in eta space) as an icosahedron with vertices
  static PolyVolume unitPolyVolume() {
    const auto t = (1.0 + sqrt(5.0)) / 2.0;
    const vector<Vector> vertsIco = {           // Array of vertex coordinates.
      Vector(-1,  t,  0)/sqrt(1 + t*t),
      Vector( 1,  t,  0)/sqrt(1 + t*t),
      Vector(-1, -t,  0)/sqrt(1 + t*t),
      Vector( 1, -t,  0)/sqrt(1 + t*t),
      Vector( 0, -1,  t)/sqrt(1 + t*t),
      Vector( 0,  1,  t)/sqrt(1 + t*t),
      Vector( 0, -1, -t)/sqrt(1 + t*t),
      Vector( 0,  1, -t)/sqrt(1 + t*t),
      Vector( t,  0, -1)/sqrt(1 + t*t),
      Vector( t,  0,  1)/sqrt(1 + t*t),
      Vector(-t,  0, -1)/sqrt(1 + t*t),
      Vector(-t,  0,  1)/sqrt(1 + t*t)
    };
    const vector<vector<unsigned>> facesIco = {
      // 5 faces around point 0
      {0, 11, 5},
      {0, 5, 1},
      {0, 1, 7},
      {0, 7, 10},
      {0, 10, 11},
      // 5 adjacent faces
      {1, 5, 9},
      {5, 11, 4},
      {11, 10, 2},
      {10, 7, 6},
      {7, 1, 8},
      // 5 faces around point 3
      {3, 9, 4},
      {3, 4, 2},
      {3, 2, 6},
      {3, 6, 8},
      {3, 8, 9},
      // 5 adjacent faces
      {4, 9, 5},
      {2, 4, 11},
      {6, 2, 10},
      {8, 6, 7},
      {9, 8, 1}
    };
    PolyVolume cell0;
    convertToPolyClipper(cell0, FacetedVolume(vertsIco, facesIco));
    ENSURE(cell0.size() == 12);
    return cell0;
  }

  // clipping operation
  static void clip(PolyVolume& cell, const std::vector<Plane>& planes) {
    PolyClipper::clipPolyhedron(cell, planes);
  }

  // moment computation
  static void moments(double& vol, Vector& cent, const PolyVolume& cell) {
    PolyClipper::moments(vol, cent, cell);
  }

  // collapse degenerate points
  static void collapseDegenerates(PolyVolume& cell, const double tol) {
    PolyClipper::collapseDegenerates(cell, tol);
  }

  // In 3D we simply use any unclipped original vertices as void generators
  static std::vector<Vector> createEtaVoidPoints(const Vector& /*etaVoidAvg*/,
                                                 const int /*nvoid*/,
                                                 const double rin,
                                                 const SymTensor& Hi,
                                                 const SymTensor& /*Hinvi*/,
                                                 const PolyVolume& celli) {
    std::vector<Vector> result;
    for (const auto& vert: celli) {
      const auto peta = Hi*vert.position;
      if (peta.magnitude2() > rin*rin) {
        result.push_back(0.5*rin*peta.unitVector());
      }
    }
    return result;
  }
    
  // toString
  static std::string toString(const PolyVolume& celli) {
    return PolyClipper::polyhedron2string(celli);
  }

};

//------------------------------------------------------------------------------
// Clip a Polygon/Polyhedron to the Voronoi about the origin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
clipHullToVoronoi(typename Dimension::FacetedVolume& poly) {
  using CT = ClippingType<Dimension>;
  using PolyVolume = typename CT::PolyVolume;
  using Plane = typename CT::Plane;

  // Convert Spheral::Poly -> PolyClipper::Poly
  PolyVolume PCpoly;
  convertToPolyClipper(PCpoly, poly);
  
  // Build the clipping planes
  vector<Plane> planes;
  for (const auto& v: PCpoly) {
    const auto d = v.position.magnitude();
    if (d > 1e-30) {  // skip if one of the vertices is the origin
      planes.push_back(Plane(0.5*d, -v.position.unitVector()));
    }
  }

  // Clip it (and collapse resulting degneracies)
  CT::clip(PCpoly, planes);
  CT::collapseDegenerates(PCpoly, 1.0e-10);

  // Convert back to Spheral poly
  convertFromPolyClipper(poly, PCpoly);
}

//..............................................................................
// 1D specialization
template<>
inline
void
clipHullToVoronoi<Dim<1>>(Dim<1>::FacetedVolume& poly) {
}

}

//------------------------------------------------------------------------------
// computeHullVolume
// The method we're actually providing
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeHullVolume(const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const ConnectivityMap<Dimension >& connectivityMap,
                  const bool clipToVoronoi,
                  FieldList<Dimension, int>& surfacePoint,
                  FieldList<Dimension, typename Dimension::Scalar>& vol,
                  FieldList<Dimension, typename Dimension::FacetedVolume>& cells) {

  TIME_FUNCTION;

  // Pre-conditions
  REQUIRE(vol.size() == position.size());

  using Vector = typename Dimension::Vector;
  using FacetedVolume = typename Dimension::FacetedVolume;
  // using PCVector = typename ClippingType<Dimension>::Vector;
  // using Plane = typename ClippingType<Dimension>::Plane;
  // using PolyVolume = typename ClippingType<Dimension>::PolyVolume;

  const auto numGens = position.numElements();
  const auto numNodeLists = position.size();
  const auto numGensGlobal = allReduce(numGens, SPHERAL_OP_SUM);
  const auto returnSurface = surfacePoint.size() == numNodeLists;
  const auto returnCells = cells.size() == numNodeLists;

  if (returnSurface) surfacePoint = 0;

  if (numGensGlobal > 0) {

    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = position[nodeListi]->numInternalElements();

    // Do each point independently
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        const auto& ri = position(nodeListi, i);

        // Build the set of inverse positions in eta space about this point (including itself as the origin)
        vector<Vector> invPositions = {Vector::zero()};
        const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, i);
        CHECK(connectivity.size() == numNodeLists);
        for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
          for (auto j: connectivity[nodeListj]) {
            const auto etaji = (position(nodeListj, j) - ri);
            invPositions.push_back(safeInv(etaji.magnitude2()) * etaji);
          }
        }

        // Compute the inverse convex hull (in 1/eta space)
        const FacetedVolume invHull(invPositions);

        // Now we can reconstruct the inner hull in proper coordinates
        auto surface = false;
        vector<Vector> verts;
        const auto& vertsInv = invHull.vertices();
        for (const auto& vinv: vertsInv) {
          const auto vimag2 = vinv.magnitude2();
          if (vimag2 < 1.0e-30) {
            verts.push_back(Vector::zero());
            surface = true;
          } else {
            verts.push_back(vinv/vimag2);
          }
        }

        // Construct the hull in normal coordinates
        FacetedVolume hull(verts);

        // If requested, clip to the Voronoi volume
        if (clipToVoronoi) clipHullToVoronoi<Dimension>(hull);

        // Put together our return values
        if (surface) {
          if (returnSurface) surfacePoint(nodeListi, i) = 1;
        } else {
          vol(nodeListi, i) = hull.volume();
        }
        if (returnCells) cells(nodeListi, i) = hull;
      }
    }
  }
}

}
