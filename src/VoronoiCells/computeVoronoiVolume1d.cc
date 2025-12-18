//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/FastMath.hh"

#include <limits>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

using namespace FastMath;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
template<>
void
computeVoronoiVolume(const FieldList<Dim<1>, Dim<1>::Vector>& position,
                     const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                     const ConnectivityMap<Dim<1> >& cm,
                     const FieldList<Dim<1>, Dim<1>::SymTensor>& damage,
                     const std::vector<Dim<1>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<1>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<1>>*>&,
                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                     FieldList<Dim<1>, int>& surfacePoint,
                     FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                     FieldList<Dim<1>, Dim<1>::Vector>& deltaMedian,
                     FieldList<Dim<1>, vector<Dim<1>::Vector>>& etaVoidPoints,
                     FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells,
                     FieldList<Dim<1>, std::vector<CellFaceFlag>>& cellFaceFlags) {

  // Pre-conditions
  REQUIRE(facetedBoundaries.size() == 0 or facetedBoundaries.size() == position.size());
  REQUIRE(vol.size() == position.size());
  REQUIRE(deltaMedian.size() == position.size());
  REQUIRE(holes.size() == facetedBoundaries.size());

  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::FacetedVolume FacetedVolume;

  const auto numGens = position.numElements();
  const auto numNodeLists = position.size();
  const auto haveFacetedBoundaries = facetedBoundaries.size() == numNodeLists;
  const auto haveWeights = weight.size() == numNodeLists;
  const auto haveDamage = damage.size() == numNodeLists;
  const auto returnSurface = surfacePoint.size() == numNodeLists;
  const auto returnCells = cells.size() == numNodeLists;
  const auto returnCellFaceFlags = cellFaceFlags.size() == numNodeLists;

  // Zero out return fields.
  deltaMedian = Vector::zero();
  if (returnSurface) {
    surfacePoint = 0;
    etaVoidPoints = vector<Vector>();
  }

  // Zero out the cell face flags.
  if (returnCellFaceFlags) {
    cellFaceFlags = vector<CellFaceFlag>();
  }

  // Copy the input positions to single list, and sort it.
  // Note our logic here relies on ghost nodes already being built, including parallel nodes.
  typedef pair<double, pair<unsigned, unsigned> > PointCoord;
  vector<PointCoord> coords;
  coords.reserve(numGens);
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = position[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      coords.push_back(make_pair(position(nodeListi, i).x(), make_pair(nodeListi, i)));
    }
  }
  sort(coords.begin(), coords.end(), ComparePairsByFirstElement<PointCoord>());

  // A local function to check if two points are neighbors
  std::set<size_t> pairHashes;
  if (returnCellFaceFlags) {
    const auto& pairs = cm.nodePairList();
    for (const auto& p: pairs) pairHashes.insert(p.hash());
  }
  auto areNeighbors = [&](const size_t il, const size_t i, const size_t jl, const size_t j) -> bool { return pairHashes.find(NodePairIdxType(i, il, j, jl).hash()) != pairHashes.end(); };

#pragma omp parallel
  {
  // Prepare some scratch variables.
  unsigned nodeListj1 = 0, nodeListj2 = 0, j1 = 0, j2 = 0;
  Scalar rin, Hi, H1, H2, x1, x2, xi, etamax, weighti, weightj, wij, xmin, xmax;
  Scalar xbound0 = -std::numeric_limits<Scalar>::max();
  Scalar xbound1 =  std::numeric_limits<Scalar>::max();

  // Now walk our sorted point and set the volumes and surface flags.
  const auto& nodeListPtrs = position.nodeListPtrs();
  const auto ntot = coords.size();
#pragma omp for
//  firstprivate(nodeListj1, nodeListj2, j1, j2,
//               rin, Hi, H1, H2, x1, x2, xi,
//               etamax, weighti, weightj, wij,
//               xbound0, xbound1)
  for (auto k = 0u; k < ntot; ++k) {
    const auto nodeListi = coords[k].second.first;
    const auto i = coords[k].second.second;
    // const bool barf = i == 0;
    if (i < nodeListPtrs[nodeListi]->firstGhostNode()) {

      // Is there a bounding volume for this NodeList?
      if (haveFacetedBoundaries) {
        xbound0 = facetedBoundaries[nodeListi].xmin().x();
        xbound1 = facetedBoundaries[nodeListi].xmax().x();
      }

      // Use the nperh to determine the cutoff distance.
      rin = 2.0/nodeListPtrs[nodeListi]->nodesPerSmoothingScale();

      // Grab our state.
      xi = position(nodeListi, i).x();
      Hi = H(nodeListi, i).xx();
      weighti = haveWeights ? weight(nodeListi, i) : 1.0;

      // phi = 1.0;
      if (k == 0) {
        x1 = xbound0 - xi;
        H1 = Hi;
        if (haveFacetedBoundaries) {
          xmin = xbound0;
        } else {
          xmin = xi - 0.5 * vol(nodeListi, i);
        }
        surfacePoint(nodeListi, i) |= 1;
        etaVoidPoints(nodeListi, i).push_back(Vector(-0.5*rin));
        if (returnCellFaceFlags) cellFaceFlags(nodeListi, i).push_back(CellFaceFlag(0, // cell face
                                                                                    -1, // void/bound
                                                                                    -1)); // void/bound
        // cerr << "Surface condition 1: " << nodeListi << " " << i << " " << surfacePoint(nodeListi, i) << endl;
      } else {
        nodeListj1 = coords[k-1].second.first;
        j1 = coords[k-1].second.second;
        H1 = H(nodeListj1, j1).xx();
        weightj = haveWeights ? weight(nodeListj1, j1) : 1.0;
        wij = weighti/(weighti + weightj);
        // the direction is here reversed from 2d, so it should weight on the i side
        x1 = wij*(position(nodeListj1, j1).x() - position(nodeListi, i).x());
        xmin = max(xbound0, x1 + xi);
        if (nodeListj1 != nodeListi) {
          surfacePoint(nodeListi, i) |= (1 << (nodeListj1 + 1));
          // cerr << "Surface condition 3: " << nodeListi << " " << i << " " << surfacePoint(nodeListi, i) << endl;
        }
        if (returnCellFaceFlags and areNeighbors(nodeListi, i, nodeListj1, j1)) cellFaceFlags(nodeListi, i).push_back(CellFaceFlag(0,           // cell face
                                                                                                                                   nodeListj1,  // other NodeList
                                                                                                                                   j1));        // other node index
      }

      if (k == ntot - 1) {
        x2 = xbound1 - xi;
        H2 = Hi;
        if (haveFacetedBoundaries) {
          xmax = xbound1;
        } else {
          xmax = xi + 0.5 * vol(nodeListi, i);
        }
        surfacePoint(nodeListi, i) |= 1;
        etaVoidPoints(nodeListi, i).push_back(Vector(0.5*rin));
        if (returnCellFaceFlags) cellFaceFlags(nodeListi, i).push_back(CellFaceFlag(1, // cell face
                                                                                    -1, // void/bound
                                                                                    -1)); // void/bound
        // cerr << "Surface condition 4: " << nodeListi << " " << i << " " << surfacePoint(nodeListi, i) << endl;
      } else {
        nodeListj2 = coords[k+1].second.first;
        j2 = coords[k+1].second.second;
        H2 = H(nodeListj2, j2).xx();
        x2 = 0.5*(position(nodeListj2, j2).x() - position(nodeListi, i).x());
        xmax = xi + x2;
        if (nodeListj2 != nodeListi) {
          surfacePoint(nodeListi, i) |= (1 << (nodeListj2 + 1));
          // cerr << "Surface condition 6: " << nodeListi << " " << i << " " << surfacePoint(nodeListi, i) << endl;
        }
        if (returnCellFaceFlags and areNeighbors(nodeListi, i, nodeListj2, j2)) cellFaceFlags(nodeListi, i).push_back(CellFaceFlag(1,           // cell face
                                                                                                                                   nodeListj2,  // other NodeList
                                                                                                                                   j2));        // other node index
      }

      CHECK(x1 <= 0.0 and x2 >= 0.0);
      // CHECK(phi >= 0.0 and phi <= 1.0);
      etamax = Hi*max(-x1, x2);
      if (etamax < rin) {
        if (surfacePoint(nodeListi, i) == 0) {
          vol(nodeListi, i) = x2 - x1;
          deltaMedian(nodeListi, i).x(0.5*(x2 - x1));
          // Check if the candidate motion is still in the boundary.  If not, project back.
          if (haveFacetedBoundaries) {
            const Vector ri = Vector(xi);
            if (not facetedBoundaries[nodeListi].contains(ri + deltaMedian(nodeListi, i))) {
              deltaMedian(nodeListi, i) = facetedBoundaries[nodeListi].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
            }
            for (unsigned ihole = 0; ihole != holes[nodeListi].size(); ++ihole) {
              if (holes[nodeListi][ihole].contains(ri + deltaMedian(nodeListi, i))) {
                deltaMedian(nodeListi, i) = holes[nodeListi][ihole].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
              }
            }
          }

          // {
          //   const Scalar x0 = -dx1,
          //     x1 = dx2,
          //     xm = deltaMedian(nodeListi, i).x(),
          //     m1 = -rhoi*x0 + 0.5*b*x0*x0 + rhoi*xm + 0.5*b*xm*xm,
          //     m2 = -rhoi*xm + 0.5*b*xm*xm + rhoi*x1 + 0.5*b*x1*x1;
          //   cerr << " --> " << i << " " << x0 << " " << x1 << " " << xm << " " << m1 << " " << m2 << endl;
          // }
        }
      } else {
        surfacePoint(nodeListi, i) |= 1;
        if (-Hi*x1 >= rin) etaVoidPoints(nodeListi, i).push_back(Vector(max(Hi*x1, -0.5*rin)));
        if ( Hi*x2 >= rin) etaVoidPoints(nodeListi, i).push_back(Vector(min(Hi*x2,  0.5*rin)));
        // cerr << "Surface condition 7: " << nodeListi << " " << i << " " << surfacePoint(nodeListi, i) << endl;
      }

      // If this point is fully damaged, we force creation of void points.
      if (haveDamage and damage(nodeListi, i).xx() > 1.0 - 1.0e-5) {
        etaVoidPoints(nodeListi, i).clear();
        etaVoidPoints(nodeListi, i).push_back(Vector(-0.5*rin));
        etaVoidPoints(nodeListi, i).push_back(Vector( 0.5*rin));
        surfacePoint(nodeListi, i) |= 1;
      }

      if (returnCells) {
        cells(nodeListi, i) = FacetedVolume({Vector(xmin), Vector(xmax)});
      }
    }
    CHECK2(((surfacePoint(nodeListi, i) & 1) == 1 and
            (etaVoidPoints(nodeListi, i).size() == 1 or etaVoidPoints(nodeListi, i).size() == 2)) or
           (surfacePoint(nodeListi, i) & 1) == 0,
           "(" << nodeListi << " " << i << ") " << xi << " " << surfacePoint(nodeListi, i) << " " << etaVoidPoints(nodeListi, i).size());

    // if (barf) cerr << "  " << i << " " << vol(nodeListi, i) << " " << surfacePoint(nodeListi, i) << " "
    //                << " ---- " << position(nodeListj1, j1).x() << " " << position(nodeListi, i) << " " << position(nodeListj2, j2).x() 
    //                << endl;

  }

  CONTRACT_VAR(H1);
  CONTRACT_VAR(H2);
  
  } // End of omp parallel region.
    
//   // Flag any points that overlap other NodeLists as multi-material.
//   if (returnSurface) {
//     for (auto nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
//       const unsigned n = position[nodeListi]->numInternalElements();
// #pragma omp parallel for
//       for (auto i = 0U; i < n; ++i) {
//         const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
//         for (auto nodeListj = 0U; nodeListj != numNodeLists; ++nodeListj) {
//           if (nodeListj != nodeListi and not fullConnectivity[nodeListj].empty()) surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));
//         }
//       }
//     }
//   }
}

}
