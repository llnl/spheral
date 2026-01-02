#include "Hydro/HydroFieldNames.hh"
#include "GaussLegendreValues.hh"
#include "SymmetricTriangularValues.hh"
#include "SymmetricTetrahedralValues.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return a value for the given coefficients at the nodal centers
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
KernelIntegrator<Dimension>::
coefficientsToValue(const FieldList<Dimension, DataType>& coeffs,
                    FieldList<Dimension, DataType>& value) const {
  VERIFY(mFlatConnectivity.indexingInitialized());
  VERIFY(mState);
  VERIFY(mState->fieldNameRegistered(HydroFieldNames::position) &&
         mState->fieldNameRegistered(HydroFieldNames::H) &&
         mState->fieldNameRegistered(HydroFieldNames::volume));
                                    
  CHECK(coeffs.size() == value.size());
  
  const auto numInternalNodes = mFlatConnectivity.numInternalNodes();
  const auto position = mState->fields(HydroFieldNames::position, Vector::zero());
  const auto H = mState->fields(HydroFieldNames::H, SymTensor::zero());
  const auto volume = mState->fields(HydroFieldNames::volume, 0.0);
  auto& kvals = mScratchData.values;
  auto& dkvals = mScratchData.dvalues;
  auto& indices = mScratchData.indices;
  auto& nodeIndices = mScratchData.nodeIndices;

  for (auto i = 0; i < numInternalNodes; ++i) {
    const auto pairi = mFlatConnectivity.localToNode(i);
    const auto nodeListi = pairi.first;
    const auto nodei = pairi.second;

    // Get the indices
    const auto numNeighbors = mFlatConnectivity.numNeighbors(i);
    mFlatConnectivity.neighborIndices(i, indices);
    CHECK(indices.size() == size_t(numNeighbors));
    nodeIndices.resize(numNeighbors);
    for (auto j = 0; j < numNeighbors; ++j) {
      nodeIndices[j] = mFlatConnectivity.localToNode(indices[j]);
    }

    // Evaluate the kernel
    const auto xi = position(nodeListi, nodei);
    kvals.resize(numNeighbors);
    dkvals.resize(numNeighbors);
    mKernel->evaluate(xi, nodeIndices, position, H, volume, mHmult,
                      kvals, dkvals);
    
    // Sum up the contributions
    auto val = DataTypeTraits<DataType>::zero();
    for (auto j = 0; j < numNeighbors; ++j) {
      const auto pairj = nodeIndices[j];
      const auto nodeListj = pairj.first;
      const auto nodej = pairj.second;
      const auto vj = volume(nodeListj, nodej);
        
      val += vj * kvals[j] * coeffs(nodeListj, nodej);
    }
    value(nodeListi, nodei) = val;
  }
}

//------------------------------------------------------------------------------
// Initialize the quadrature
//------------------------------------------------------------------------------
template<>
inline
void
KernelIntegrator<Dim<1>>::
initializeQuadrature() {
  mNumOrdinates = GaussLegendreValues::numOrdinatesForOrder(mIntegrationOrder);
  GaussLegendreValues::getQuadrature(mNumOrdinates,
                                     mBaseWeights,
                                     mBaseOrdinates);
  CHECK(mBaseWeights.size() == size_t(mNumOrdinates) &&
        mBaseOrdinates.size() == size_t(mNumOrdinates));
  mNumSurfaceOrdinates = 1;
  mBaseSurfaceWeights = {1.};
  mBaseSurfaceOrdinates = {Vector(0.)};
  CHECK(mBaseSurfaceWeights.size() == size_t(mNumSurfaceOrdinates) &&
        mBaseSurfaceOrdinates.size() == size_t(mNumSurfaceOrdinates));
}

template<>
inline
void
KernelIntegrator<Dim<2>>::
initializeQuadrature() {
  mNumOrdinates = SymmetricTriangularValues::numOrdinatesForOrder(mIntegrationOrder);
  SymmetricTriangularValues::getQuadrature(mNumOrdinates,
                                           mBaseWeights,
                                           mBaseOrdinates);
  CHECK(mBaseWeights.size() == size_t(mNumOrdinates) &&
        mBaseOrdinates.size() == size_t(mNumOrdinates));
  mNumSurfaceOrdinates = GaussLegendreValues::numOrdinatesForOrder(mIntegrationOrder);
  GaussLegendreValues::getQuadrature(mNumSurfaceOrdinates,
                                     mBaseSurfaceWeights,
                                     mBaseSurfaceOrdinates);
  CHECK(mBaseSurfaceWeights.size() == size_t(mNumSurfaceOrdinates) &&
        mBaseSurfaceOrdinates.size() == size_t(mNumSurfaceOrdinates));
}

template<>
inline
void
KernelIntegrator<Dim<3>>::
initializeQuadrature() {
  mNumOrdinates = SymmetricTetrahedralValues::numOrdinatesForOrder(mIntegrationOrder);
  SymmetricTetrahedralValues::getQuadrature(mNumOrdinates,
                                            mBaseWeights,
                                            mBaseOrdinates);
  CHECK(mBaseWeights.size() == size_t(mNumOrdinates) &&
        mBaseOrdinates.size() == size_t(mNumOrdinates));
  mNumSurfaceOrdinates = SymmetricTriangularValues::numOrdinatesForOrder(mIntegrationOrder);
  SymmetricTriangularValues::getQuadrature(mNumSurfaceOrdinates,
                                           mBaseSurfaceWeights,
                                           mBaseSurfaceOrdinates);
  CHECK(mBaseSurfaceWeights.size() == size_t(mNumSurfaceOrdinates) &&
        mBaseSurfaceOrdinates.size() == size_t(mNumSurfaceOrdinates));
}

// //------------------------------------------------------------------------------
// // Get the subfacets for the surface
// //------------------------------------------------------------------------------
// template<>
// inline
// void
// KernelIntegrator<Dim<1>>::
// getSubfacets(const Facet& facet,
//              std::vector<Subfacet>& subfacets) {
//   subfacets = {{facet.point()}};
// }

// template<>
// inline
// void
// KernelIntegrator<Dim<2>>::
// getSubfacets(const Facet& facet,
//              std::vector<Subfacet>& subfacets) {
//   subfacets = {{facet.point1(), facet.point2()}};
// }

// template<>
// inline
// void
// KernelIntegrator<Dim<3>>::
// getSubfacets(const Facet& facet,
//              std::vector<Subfacet>& subfacets) {
//   const auto numPoints = facet.ipoints().size();
//   switch (numPoints) {
//   case 3:
//     // Return the input triangle
//     subfacets = {{facet.point(0), facet.point(1), facet.point(2)}};
//     break;
//   case 4:
//     // Split the quadrilateral into two triangles
//     subfacets = {{facet.point(0), facet.point(1), facet.point(2)},
//                  {facet.point(2), facet.point(3), facet.point(0)}};
//     break;
//   default:
//     // Create centroid and make a new facet for each point
//     const auto centroid = facet.position();
//     subfacets.resize(numPoints);
//     for (auto i = 0; i < numPoints; ++i) {
//       subfacets[i] = {centroid,
//                       facet.point(i),
//                       facet.point((i+1) % numPoints)};
//     }
//     break;
//   }
// }

// //------------------------------------------------------------------------------
// // Get the subcells for a given faceted volume
// //------------------------------------------------------------------------------
// template<>
// inline
// void
// KernelIntegrator<Dim<1>>::
// getSubcells(const FacetedVolume& cell,
//             std::vector<FacetedVolume>& subcells) {
//   subcells = {cell};
// }

// template<>
// inline
// void
// KernelIntegrator<Dim<2>>::
// getSubcells(const FacetedVolume& cell,
//             std::vector<FacetedVolume>& subcells) {
//   // This version was having issues
//   // const auto& vertices = cell.vertices();
//   // const auto numVertices = vertices.size();
//   // switch (numVertices) {
//   // case 3:
//   //   // If this is a triangle, we don't need to split it
//   //   subcells = {cell};
//   //   break;
//   // case 4:
//   //   {
//   //     // If this is a quadrilaterial, split into two triangles
//   //     const std::vector<Vector> points1 = {vertices[0], vertices[1], vertices[2]};
//   //     const std::vector<Vector> points2 = {vertices[2], vertices[3], vertices[0]};
//   //     subcells = {FacetedVolume(points1), FacetedVolume(points2)};
//   //   }
//   //   break;
//   // default:
//   //   {
//   //     // Otherwise, create a centroid and make a new triangle for each facet
//   //     const auto& facets = cell.facets();
//   //     const auto numFacets = facets.size();
//   //     const auto centroid = cell.centroid();
//   //     subcells.resize(numFacets);
//   //     for (auto f = 0; f < numFacets; ++f) {
//   //       const auto& facet = facets[f];
//   //       const std::vector<Vector> points = {centroid, facet.point1(), facet.point2()};
//   //       subcells[f] = FacetedVolume(points);
//   //     }
//   //   }
//   //   break;
//   // }

//   // Use the normal Spheral version
//   const auto numFacets = cell.facets().size();
//   subcells.resize(numFacets);
//   for (auto f = 0; f < numFacets; ++f) {
//     subcells[f] = cell.facetSubVolume(f);
//   }

//   BEGIN_CONTRACT_SCOPE
//   {
//     const auto volume = cell.volume();
//     auto subvolume = 0.0;
//     for (auto& subcell : subcells) {
//       subvolume += subcell.volume();
//       CHECK(subcell.facets().size() == 3);
//     }
    
//     ENSURE(fuzzyEqual(volume, subvolume));
//   }
//   END_CONTRACT_SCOPE
// }

// template<>
// inline
// void
// KernelIntegrator<Dim<3>>::
// getSubcells(const FacetedVolume& cell,
//             std::vector<FacetedVolume>& subcells) {
//   const auto& vertices = cell.vertices();
//   const auto numVertices = vertices.size();
//   switch (numVertices) {
//   // case 4:
//   //   // We already have a tetrahedron
//   //   subcells = {cell};
//   //   break;
//   default:
//     {
//       // Split each facet into triangles and create the subcells using these triangles and the centroid
//       const auto centroid = cell.centroid();
//       const auto& facets = cell.facets();
//       const auto numFacets = facets.size();
//       subcells.clear();
//       subcells.reserve(numFacets);
//       for (auto f = 0; f < numFacets; ++f) {
//         const auto& facet = facets[f];
//         std::vector<Subfacet> subfacets;
//         getSubfacets(facet, subfacets);
        
//         for (auto& subfacet : subfacets) {
//           std::vector<Vector> points = {centroid, subfacet[0], subfacet[1], subfacet[2]};
//           subcells.push_back(FacetedVolume(points));
//         }
//       }
//     }
//     break;
//   }
  
//   // // The normal Spheral version isn't working, since it produces pyramids
//   // const auto numFacets = cell.facets().size();
//   // subcells.resize(numFacets);
//   // for (auto f = 0; f < numFacets; ++f) {
//   //   subcells[f] = cell.facetSubVolume(f);
//   // }

//   BEGIN_CONTRACT_SCOPE
//   {
//     const auto volume = cell.volume();
//     auto subvolume = 0.0;
//     for (auto& subcell : subcells) {
//       subvolume += subcell.volume();
//       CHECK(subcell.facets().size() == 4);
//     }
    
//     ENSURE(fuzzyEqual(volume, subvolume));
//   }
//   END_CONTRACT_SCOPE
// }

//------------------------------------------------------------------------------
// Get a quadrature for the given geometry
// This must be a line segment, triangle, or tetrahedron
// Base quadrature in 1d is a line segment (-1, 1)
// Base quadrature in 2d is a triangle (0,0), (1,0), (0,1)
// Base quadrature in 3d is a tetrahedron (0,0,0), (1,0,0), (0,1,0), (0,0,1)
//------------------------------------------------------------------------------
template<>
inline
void
KernelIntegrator<Dim<1>>::
getQuadrature(const FacetedVolume& region,
              std::vector<Scalar>& weights,
              std::vector<Vector>& ordinates) {
  CHECK(weights.size() == size_t(mNumOrdinates) &&
        ordinates.size() == size_t(mNumOrdinates) &&
        mBaseWeights.size() == size_t(mNumOrdinates) &&
        mBaseOrdinates.size() == size_t(mNumOrdinates));
  CHECK(region.vertices().size() == 2);
  const auto p0 = region.xmin();
  const auto p1 = region.xmax();
  const auto J = Tensor(0.5 * (p1.x() - p0.x()));
  const auto Jdet = std::abs(J.Determinant());
  for (auto i = 0; i < mNumOrdinates; ++i) {
    weights[i] = mBaseWeights[i] * Jdet;
    ordinates[i] = p0 + J * (mBaseOrdinates[i] + Vector::one());
  }
}

template<>
inline
void
KernelIntegrator<Dim<2>>::
getQuadrature(const FacetedVolume& region,
              std::vector<Scalar>& weights,
              std::vector<Vector>& ordinates) {
  CHECK(weights.size() == size_t(mNumOrdinates) &&
        ordinates.size() == size_t(mNumOrdinates) &&
        mBaseWeights.size() == size_t(mNumOrdinates) &&
        mBaseOrdinates.size() == size_t(mNumOrdinates));
  const auto& vertices = region.vertices();
  CHECK(vertices.size() == 3);
  const auto& p0 = vertices[0];
  const auto& p1 = vertices[1];
  const auto& p2 = vertices[2];
  const auto J = Tensor(p1.x() - p0.x(), p2.x() - p0.x(),
                        p1.y() - p0.y(), p2.y() - p0.y());
  const auto Jdet = std::abs(J.Determinant());
  for (auto i = 0; i < mNumOrdinates; ++i) {
    weights[i] = mBaseWeights[i] * Jdet;
    ordinates[i] = p0 + J * mBaseOrdinates[i];
  }
}

template<>
inline
void
KernelIntegrator<Dim<3>>::
getQuadrature(const FacetedVolume& region,
              std::vector<Scalar>& weights,
              std::vector<Vector>& ordinates) {
  CHECK(weights.size() == size_t(mNumOrdinates) &&
        ordinates.size() == size_t(mNumOrdinates) &&
        mBaseWeights.size() == size_t(mNumOrdinates) &&
        mBaseOrdinates.size() == size_t(mNumOrdinates));
  const auto& vertices = region.vertices();
  CHECK(vertices.size() == 4);
  const auto& p0 = vertices[0];
  const auto& p1 = vertices[1];
  const auto& p2 = vertices[2];
  const auto& p3 = vertices[3];
  const auto J = Tensor(p1.x() - p0.x(), p2.x() - p0.x(), p3.x() - p0.x(),
                        p1.y() - p0.y(), p2.y() - p0.y(), p3.y() - p0.y(),
                        p1.z() - p0.z(), p2.z() - p0.z(), p3.z() - p0.z());
  const auto Jdet = std::abs(J.Determinant());
  for (auto i = 0; i < mNumOrdinates; ++i) {
    weights[i] = mBaseWeights[i] * Jdet;
    ordinates[i] = p0 + J * mBaseOrdinates[i];
  }
}

//------------------------------------------------------------------------------
// Get a surface quadrature for the given geometry
// The input is a set of points defining the surface with size Dimension::nDim
// See the getQuadrature method for the base geometries
//------------------------------------------------------------------------------
template<>
inline
void
KernelIntegrator<Dim<1>>::
getSurfaceQuadrature(const Subfacet& subfacet,
                     std::vector<Scalar>& weights,
                     std::vector<Vector>& ordinates) {
  weights = {1.0};
  ordinates = {subfacet[0]};

  ENSURE(weights.size() == size_t(mNumSurfaceOrdinates));
  ENSURE(ordinates.size() == size_t(mNumSurfaceOrdinates));
}

// \bm{x}\left(\xi\right)=\bm{p}_{0}+\left(\bm{p}_{1}-\bm{p}_{0}\right)\frac{\left(\xi+1\right)}{2}
// \bm{n}=\frac{\partial\bm{x}}{\partial\xi}
template<>
inline
void
KernelIntegrator<Dim<2>>::
getSurfaceQuadrature(const Subfacet& subfacet,
                     std::vector<Scalar>& weights,
                     std::vector<Vector>& ordinates) {
  CHECK(weights.size() == size_t(mNumSurfaceOrdinates));
  CHECK(ordinates.size() == size_t(mNumSurfaceOrdinates));
  const auto p0 = subfacet[0];
  const auto p1 = subfacet[1];
  const auto p10h = 0.5 * (p1 - p0);
  const auto normalMag = p10h.magnitude();
  for (auto i = 0; i < mNumSurfaceOrdinates; ++i) {
    weights[i] = mBaseSurfaceWeights[i] * normalMag;
    ordinates[i] = p0 + p10h * (1. + mBaseSurfaceOrdinates[i][0]);
  }
}

// \bm{x}\left(\bm{\xi}\right)=\bm{p}_{0}+\xi_{1}\left(\bm{p}_{1}-\bm{p}_{0}\right)+\xi_{2}\left(\bm{p}_{2}-\bm{p}_{0}\right)
// \bm{n}=\frac{\partial\bm{x}}{\partial\xi_{1}}\times\frac{\partial\bm{x}}{\partial\xi_{2}}
template<>
inline
void
KernelIntegrator<Dim<3>>::
getSurfaceQuadrature(const Subfacet& subfacet,
                     std::vector<Scalar>& weights,
                     std::vector<Vector>& ordinates) {
  CHECK(weights.size() == size_t(mNumSurfaceOrdinates) &&
        ordinates.size() == size_t(mNumSurfaceOrdinates) &&
        mBaseSurfaceWeights.size() == size_t(mNumSurfaceOrdinates) &&
        mBaseSurfaceOrdinates.size() == size_t(mNumSurfaceOrdinates));
  const auto p0 = subfacet[0];
  const auto p1 = subfacet[1];
  const auto p2 = subfacet[2];
  const auto p10 = p1 - p0;
  const auto p20 = p2 - p0;
  const auto normal = p10.cross(p20);
  const auto normalMag = normal.magnitude();
  for (auto i = 0; i < mNumSurfaceOrdinates; ++i) {
    weights[i] = mBaseSurfaceWeights[i] * normalMag;
    ordinates[i] = p0 + p10 * mBaseSurfaceOrdinates[i][0] + p20 * mBaseSurfaceOrdinates[i][1];
  }
}

} // end namespace Spheral
