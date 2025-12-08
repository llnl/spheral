//------------------------------------------------------------------------------
// Compute the Hull mass density summation.
//------------------------------------------------------------------------------
#include "computeHullSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/comparisons.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/DamagedNodeCouplingWithFrags.hh"

#include "polytope/polytope.hh"
#include "polytope/convexHull_2d.hh"

namespace Spheral {

namespace {
#ifdef SPHERAL1D
//------------------------------------------------------------------------------
// Worker methods specialized by dimension to compute the mass density based
// on convex hulls.
//------------------------------------------------------------------------------
// 1D
inline
double hullMassDensity(const std::vector<Dim<1>::Vector>& posInv,
                       const std::vector<Dim<1>::Scalar>& mass) {
  REQUIRE(posInv.size() == mass.size());
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;

  if (posInv.size() < 2) return 0.0;

  // Remember the mass of point i which comes in as the first value.
  const Scalar mi = mass[0];

  // Copy the two vectors to a single vector<pair> for sorting.
  vector<pair<Vector, Scalar> > stuff;
  const unsigned n = posInv.size();
  for (unsigned i = 0; i != n; ++i) stuff.push_back(make_pair(posInv[i], mass[i]));
  sort(stuff.begin(), stuff.end(), ComparePairByFirstElement<pair<Vector, Scalar> >());

  // Add up the masses.
  Scalar msum = 0.5*(stuff[0].second + stuff[n-1].second);
  if (stuff[0].first != Vector::zero() and stuff[n-1].first != Vector::zero()) msum += mi;

  // Figure out the volume.
  const Scalar vol = safeInv(stuff[n-1].first.x()) - safeInv(stuff[0].first.x());
  CHECK(vol > 0.0);

  // We've got it.
  return msum*safeInv(vol);
}
#endif

//..............................................................................
// 2D
//..............................................................................
#ifdef SPHERAL2D
inline
double hullMassDensity(const std::vector<Dim<2>::Vector>& posInv,
                       const std::vector<Dim<2>::Scalar>& mass) {
  REQUIRE(posInv.size() == mass.size());
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Find the appropriate renormalization so that we can do the convex hull
  // in a unit box.
  Vector xmin, xmax;
  boundingBox(posInv, xmin, xmax);
  const double fscale = (xmax - xmin).maxElement();
  CHECK(fscale > 0.0);

  // Copy the point coordinates to a polytope point array.
  const unsigned n = posInv.size();
  vector<double> points_polytope;
  vector<Vector> pos;
  points_polytope.reserve(2*n);
  pos.reserve(n);
  for (const Vector& vec: posInv) {
    points_polytope.push_back((vec.x() - xmin.x())/fscale);
    points_polytope.push_back((vec.y() - xmin.y())/fscale);
    const Scalar mag2 = vec.magnitude2();
    if (mag2 == 0.0) {
      pos.push_back(Vector::zero());
    } else {
      pos.push_back(vec/mag2);
    }
  }
  CHECK(points_polytope.size() == 2*n);
  CHECK(pos.size() == n);

  // Call the polytope method for computing the convex hull.
  vector<double> low(2, 0.0);
  polytope::PLC<2, double> plc = polytope::convexHull_2d(points_polytope, &(*low.begin()), 1.0e-15);
  const unsigned numVertices = plc.facets.size();
  CHECK(numVertices >= 3);

  // Build a polygon in non-inverse coordinates.
  Scalar msum = 0.0;
  vector<unsigned> flags(n, 0);
  vector<Vector> verts(numVertices);
  vector<vector<unsigned> > facets(numVertices, vector<unsigned>(2));
  for (unsigned k = 0; k != numVertices; ++k) {
    const unsigned i = plc.facets[k][0],
                   j = plc.facets[k][1],
                   kk = (k + 1) % numVertices;
    CHECK(i < n);
    CHECK(j < n);
    CHECK(kk < numVertices and plc.facets[kk][0] == j);
    verts[k] = pos[i];
    facets[k][0] = k;
    facets[k][1] = kk;

    // Each point on the hull contributes a fraction of mass.
    const Vector v1 = pos[plc.facets[kk][1]] - pos[plc.facets[kk][0]],
                 v2 = pos[i] - pos[j],
                 v3 = pos[plc.facets[kk][1]] - pos[j];
    // const Scalar theta = atan2(v3.x(), v3.y());
    const Scalar theta = acos(max(-1.0, min(1.0, v1.dot(v2)*safeInv(v1.magnitude()*v2.magnitude()))));
    msum += theta/(2.0*M_PI) * mass[j];
    flags[j] = 1;
  }
  const FacetedVolume hull(verts, facets);

  // Check for any points contained in the hull and add their mass.
  // We exclude boundary points here since we took care of them in the
  // previous loop.
  for (unsigned i = 0; i != n; ++i) {
    if (flags[i] == 0) {
      if (pointOnPolygon(pos[i], verts)) {
        msum += 0.5*mass[i];
      } else if (hull.contains(pos[i], false)) {
        msum += mass[i];
      }
    }
  }

  // And finally we're there.
  const Scalar vol = hull.volume();
  CHECK2(msum > 0.0 and vol > 0.0, "Bad density estimate: " << msum << " " << vol);
  return msum*safeInv(vol);
}
#endif

//..............................................................................
// 3D
//..............................................................................
#ifdef SPHERAL3D
inline
double hullMassDensity(const std::vector<Dim<3>::Vector>& posInv,
                       const std::vector<Dim<3>::Scalar>& mass) {
  REQUIRE(posInv.size() == mass.size());
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  VERIFY(false);
}
#endif

}

//------------------------------------------------------------------------------
// The public method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeHullSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                          const TableKernel<Dimension>& W,
                          const FieldList<Dimension, typename Dimension::Vector>& position,
                          const FieldList<Dimension, typename Dimension::Scalar>& mass,
                          const FieldList<Dimension, typename Dimension::SymTensor>& H,
                          const NodeCoupling& nodeCoupling,
                          FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  const Scalar kernelExtent2 = FastMath::square(W.kernelExtent());

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar rhoMin = nodeList.rhoMin();
    const Scalar rhoMax = nodeList.rhoMax();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;
      if (i < firstGhostNodei) {

        // Get the state for node i.
        const Vector& ri = position(nodeListi, i);
        const Scalar mi = mass(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);

        // Get the neighbors for this node (in this NodeList).  We use the approximation here
        // that nodes from other NodeLists do not contribute to the density of this one.
        const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListi, i)[nodeListi];

        // Copy the neighbor positions & masses.
        vector<Vector> positionsInv(1, Vector::zero());
        vector<Scalar> masses(1, mi);
        positionsInv.reserve(connectivity.size() + 1);
        masses.reserve(connectivity.size() + 1);
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const unsigned j = *jItr;
          if (nodeCoupling(NodePairIdxType(i, nodeListi, j, nodeListi)) > 0.0) {
            const Vector& rj = position(nodeListi, j);
            const Vector rji = rj - ri;
            const Scalar etai2 = (Hi*rji).magnitude2();
            if (etai2 < kernelExtent2) {
              const Vector rjiHat = rji.unitVector();
              positionsInv.push_back(1.0/sqrt(rji.magnitude2() + 1.0e-30) * rjiHat);
              masses.push_back(mass(nodeListi, *jItr));
            }
          }
        }
        CHECK(masses.size() == positionsInv.size());

        // Delegate to specialized methods.
        massDensity(nodeListi, i) = max(rhoMin, min(rhoMax, hullMassDensity(positionsInv, masses)));
        CHECK(massDensity(nodeListi, i) > 0.0);
      }
    }
  }
}

}

