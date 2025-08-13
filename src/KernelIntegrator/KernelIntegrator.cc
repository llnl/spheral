//---------------------------------Spheral++----------------------------------//
// KernelIntegrator
//
// Performs integrals of kernels and coefficients
//----------------------------------------------------------------------------//
#include "KernelIntegrator.hh"

#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
KernelIntegrator<Dimension>::
KernelIntegrator(const int integrationOrder,
                 const std::shared_ptr<IntegrationKernel<Dimension>> kernel,
                 const DataBase<Dimension>& dataBase,
                 const FlatConnectivity<Dimension>& flatConnectivity):
  mIntegrationOrder(integrationOrder),
  mKernel(kernel),
  mDataBase(dataBase),
  mFlatConnectivity(flatConnectivity),
#if REPLACEOVERLAP
  mHmult(2.0),
#else
  mHmult(1.0),
#endif
  mNumOrdinates(),
  mBaseWeights(),
  mBaseOrdinates(),
  mNumSurfaceOrdinates(),
  mBaseSurfaceWeights(),
  mBaseSurfaceOrdinates(),
  mState(),
  mIntegrals(),
  mScratchData(),
  mTotalNumSubcells(),
  mTotalNumSubfacets() {
  initializeQuadrature();
}

//------------------------------------------------------------------------------
// Set the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
KernelIntegrator<Dimension>::
setState(const double time,
         const State<Dimension>& state) {
  mTime = time;
  mState = std::make_unique<State<Dimension>>(state);
  VERIFY(mState->fieldNameRegistered(HydroFieldNames::position) &&
         mState->fieldNameRegistered(HydroFieldNames::H) &&
         mState->fieldNameRegistered(HydroFieldNames::volume) &&
         mState->fieldNameRegistered(HydroFieldNames::cells) &&
         mState->fieldNameRegistered(HydroFieldNames::cellFaceFlags));
}

//------------------------------------------------------------------------------
// Add integrals 
//------------------------------------------------------------------------------
template<typename Dimension>
void
KernelIntegrator<Dimension>::
addIntegral(std::shared_ptr<KernelIntegralBase<Dimension>> integral) {
  for (auto val : mIntegrals) {
    if (integral == val) {
      std::cerr << "KernelIntegrator: tried to add integral that already exists" << std::endl;
      return;
    }
  }
  mIntegrals.push_back(integral);
}

//------------------------------------------------------------------------------
// Perform integration
//------------------------------------------------------------------------------
template<typename Dimension>
void
KernelIntegrator<Dimension>::
performIntegration() {
  VERIFY(mFlatConnectivity.indexingInitialized());
  VERIFY2(mIntegrals.size() > 0, "there are no integrals to do");
#pragma omp parallel
  VERIFY2(omp_get_num_threads() == 1, "integration fails for > 1 OpenMP thread");
  
  // Get some data out of database and state
  VERIFY(mState);
  const auto numNodeLists = mDataBase.numFluidNodeLists();
  CONTRACT_VAR(numNodeLists);
  const auto position = mState->fields(HydroFieldNames::position, Vector::zero);
  const auto H = mState->fields(HydroFieldNames::H, SymTensor::zero);
  const auto volume = mState->fields(HydroFieldNames::volume, 0.0);
  const auto cells = mState->fields(HydroFieldNames::cells, FacetedVolume());

  // Since we may be receiving state from user, we need to make sure FieldLists aren't empty
  VERIFY(position.size() >= numNodeLists &&
         H.size() >= numNodeLists &&
         cells.size() == numNodeLists);
  
  // Initialize the integration quadrature for a single triangle
  std::vector<Scalar> weights(mNumOrdinates);
  std::vector<Vector> ordinates(mNumOrdinates);
  std::vector<Scalar> surfaceWeights(mNumSurfaceOrdinates);
  std::vector<Vector> surfaceOrdinates(mNumSurfaceOrdinates);
  ArrayDim normalArray;

  // Initialize the subcells, which are needed so we don't have to specialize in dimension
  std::vector<FacetedVolume> subcells;
  std::vector<Facet> facets;
  std::vector<Subfacet> subfacets;
  
  // Initialize the integration data
  KernelIntegrationData<Dimension> data;
  data.time = mTime;

  // Zero out the integrals
  for (auto integral : mIntegrals) {
    integral->initialize(mFlatConnectivity);
  }

  // Get information about the integrals
  auto haveBilinearIntegral = false;
  auto haveVolumeIntegral = false;
  auto haveSurfaceIntegral = false;
  for (auto integral : mIntegrals) {
    if (integral->bilinear()) haveBilinearIntegral = true;
    if (integral->volume())   haveVolumeIntegral = true;
    if (integral->surface())  haveSurfaceIntegral = true;
  }

#if !REPLACEOVERLAP
  if (haveBilinearIntegral) {
    VERIFY(mFlatConnectivity.overlapIndexingInitialized());
  }
#endif
  if (haveSurfaceIntegral) {
    VERIFY(mFlatConnectivity.surfaceIndexingInitialized());
  }
  
  // Zero out the integral count
  mTotalNumSubcells = 0;
  mTotalNumSubfacets = 0;
  
  // Top level loop is over Voronoi integration regions
#if REPLACEOVERLAP
  std::vector<int> indicesTemp;
  // const auto extent = mKernel->extent(mHmult);
#endif
  const auto numNodes = mFlatConnectivity.numNodes();
  for (auto i = 0; i < numNodes; ++i) {
    const auto pairi = mFlatConnectivity.localToNode(i);
    const auto nodeListi = pairi.first;
    const auto nodei = pairi.second;
// #if REPLACEOVERLAP
//     const auto xi = position(nodeListi, nodei);
//     const auto volumei = volume(nodeListi, nodei);
//     const auto deltai = std::pow(volumei, 1. / Dimension::nDim); // An approximation of cell length
// #endif
    const auto& cell = cells(nodeListi, nodei);
    
    // Get the number of neighbors for the Voronoi integration region
    const auto numNeighborsIncSelf = mFlatConnectivity.numNeighbors(i);
    
    // Fill the indices: these are the standard connectivity, not the overlap
    data.index0 = i;
    data.nodeIndex0 = pairi;
// #if REPLACEOVERLAP
//     mFlatConnectivity.neighborIndices(i, indicesTemp);
//     CHECK(indicesTemp.size() == numNeighborsIncSelf);
//     data.indices.clear();
//     data.nodeIndices.clear();
//     for (auto j = 0; j < numNeighborsIncSelf; ++j) {
//       const auto pairj = mFlatConnectivity.localToNode(indicesTemp[j]);
//       const auto nodeListj = pairj.first;
//       const auto nodej = pairj.second;
//       const auto Hj = H(nodeListj, nodej);
//       const auto xj = position(nodeListj, nodej);

//       // Add the cell half-length to the extent to find overlap w/ cell instead of intersection w/ point
//       const auto deltaj = deltai * Hj.eigenValues().maxElement();
//       const auto extentj = extent + deltaj;

//       // Check whether support j intersects with cell i
//       // if (mHmult * (Hj * (xi - xj)).magnitude() <= extentj) {
//         data.indices.push_back(indicesTemp[j]);
//         data.nodeIndices.push_back(pairj);
//       // }
//     }

//     // Update the number of neighbors to be the support-cell intersections
//     const auto numNeighborsRefined = data.indices.size();
//     CHECK(data.nodeIndices.size() == numNeighborsRefined);
//     data.values.resize(numNeighborsRefined);
//     data.dvalues.resize(numNeighborsRefined);
// #else
    mFlatConnectivity.neighborIndices(i, data.indices);
    CHECK(data.indices.size() == size_t(numNeighborsIncSelf));
    data.nodeIndices.resize(numNeighborsIncSelf);
    data.volume.resize(numNeighborsIncSelf);
    for (auto j = 0; j < numNeighborsIncSelf; ++j) {
      const auto pairj = mFlatConnectivity.localToNode(data.indices[j]);
      const auto nodeListj = pairj.first;
      const auto nodej = pairj.second;
      data.nodeIndices[j] = pairj;
      data.volume[j] = volume(nodeListj, nodej);
    }
    data.values.resize(numNeighborsIncSelf);
    data.dvalues.resize(numNeighborsIncSelf);
// #endif
    
    // Fill the bilinear mapping for integrals
    // This is a lot of storage, but otherwise we end up doing it a lot of times
    if (haveBilinearIntegral) {
#if REPLACEOVERLAP
      // data.localIndex.resize(numNeighborsRefined * numNeighborsRefined);
      // for (auto j = 0; j < numNeighborsRefined; ++j) {
      //   for (auto k = 0; k < numNeighborsRefined; ++k) {
      //     data.localIndex[k + numNeighborsRefined*j] = mFlatConnectivity.localToFlat(data.indices[j],
      //                                                                                data.indices[k]);
      //   }
      // }
      data.localIndex.resize(numNeighborsIncSelf * numNeighborsIncSelf);
      for (auto j = 0; j < numNeighborsIncSelf; ++j) {
        for (auto k = 0; k < numNeighborsIncSelf; ++k) {
          data.localIndex[k + numNeighborsIncSelf*j] = mFlatConnectivity.localToFlat(data.indices[j],
                                                                                     data.indices[k]);
        }
      }
#else
      data.localIndex.resize(numNeighborsIncSelf * numNeighborsIncSelf);
      for (auto j = 0; j < numNeighborsIncSelf; ++j) {
        for (auto k = 0; k < numNeighborsIncSelf; ++k) {
          data.localIndex[k + numNeighborsIncSelf*j] = mFlatConnectivity.localToFlatOverlap(data.indices[j],
                                                                                            data.indices[k]);
        }
      }
#endif
    }

    // if (Process::getRank() == 0) {
    //   std::cout << i << "\t";
    //   for (int j = 0; j < numNeighborsIncSelf; ++j) {
    //     std::cout << data.indices[j] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    
    // Decompose the Voronoi integration region into subcells
    if (haveVolumeIntegral) {
      cell.decompose(subcells);
      // getSubcells(cell, subcells);
      const auto numSubcells = subcells.size();
      mTotalNumSubcells += numSubcells;
      for (size_t c = 0; c < numSubcells; ++c) {
        // For each subcell, create a quadrature and loop over integration points
        const auto& subcell = subcells[c];
        getQuadrature(subcell, weights, ordinates);
        for (auto q = 0; q < mNumOrdinates; ++q) {
          data.weight = weights[q];
          data.ordinate = ordinates[q];
          mKernel->evaluate(data.ordinate, data.nodeIndices, position, H, volume, mHmult, data.values, data.dvalues);
          
          for (auto integral : mIntegrals) {
            integral->addToIntegral(data);
          }
        } // for q
      } // for c
    } // if haveVolumeIntegral
      
      // Do the surface integration
    if (haveSurfaceIntegral) {
      // Loop over the void surfaces for this Voronoi cell, if any
      const auto numSurfacesForCell = mFlatConnectivity.numSurfacesForCell(i);
      if (numSurfacesForCell > 0) {
// #if REPLACEOVERLAP
//         data.surfaceIndex.resize(numNeighborsRefined);
//         data.numSurfaces.resize(numNeighborsRefined);
// #else
        data.surfaceIndex.resize(numNeighborsIncSelf);
        data.numSurfaces.resize(numNeighborsIncSelf);
// #endif        
        const auto& facets = cell.facets();
        for (auto flats = 0; flats < numSurfacesForCell; ++flats) {
          const auto f = mFlatConnectivity.surfaceIndexForCell(i, flats);
          CHECK(size_t(f) < facets.size());
          const auto& facet = facets[f];
          const auto& normalArea = facet.normal();
          data.normal = normalArea.unitVector();
          mFlatConnectivity.normalToArray(data.normal,
                                          normalArray);
          facet.decompose(subfacets);
          
          // Get the surface indices for the integrals
          getSurfaceIndices(normalArray,
                            data.indices,
                            data.surfaceIndex,
                            data.numSurfaces);
          CHECK(data.indices[0] == i);
          data.surfaceIndex0 = data.surfaceIndex[0];
          const auto numSubfacets = subfacets.size();
          mTotalNumSubfacets += numSubfacets;
          for (auto s = 0u; s < numSubfacets; ++s) {
            const auto& subfacet = subfacets[s];
            getSurfaceQuadrature(subfacet, surfaceWeights, surfaceOrdinates);
            for (auto q = 0; q < mNumSurfaceOrdinates; ++q) {
              // Evaluate functions at quadrature point
              data.weight = surfaceWeights[q];
              data.ordinate = surfaceOrdinates[q];
              mKernel->evaluate(data.ordinate, data.nodeIndices, position, H, volume, mHmult, data.values, data.dvalues);
              
              // Perform integration for all surface integrals
              for (auto integral : mIntegrals) {
                integral->addToSurfaceIntegral(data);
              }
            } // for q
          } // for s
        } // for c
      } // if flags.size() > 0
    } // if haveSurfaceIntegral
  } // for locali

  // Perform any post-integration work
  for (auto integral : mIntegrals) {
    integral->finalize(mFlatConnectivity);
  }
}

//------------------------------------------------------------------------------
// Get surface indices for each point
//------------------------------------------------------------------------------
template<typename Dimension>
void
KernelIntegrator<Dimension>::
getSurfaceIndices(const ArrayDim& normal,
                  const std::vector<int>& indices,
                  std::vector<int>& surfaceIndex,
                  std::vector<int>& numSurfaces) {
  const auto numIndices = indices.size();
  for (auto i = 0u; i < numIndices; ++i) {
    surfaceIndex[i] = mFlatConnectivity.surfaceIndex(indices[i],
                                                     normal);
    numSurfaces[i] = mFlatConnectivity.numSurfaces(indices[i]);
  }
}

} // end namespace Spheral
