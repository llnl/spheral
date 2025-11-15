//------------------------------------------------------------------------------
// Implement Lloyd's algorithm for centroidal relaxation of fluid points.
//------------------------------------------------------------------------------
#include "centroidalRelaxNodesImpl.hh"
#include "VoronoiCells/computeVoronoiVolume.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/gradientRK.hh"

#include <ctime>
using std::vector;
using std::tuple;
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

template<typename Dimension>
unsigned
centroidalRelaxNodesImpl(DataBase<Dimension>& db,
                         const std::vector<typename Dimension::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<typename Dimension::FacetedVolume> >& holes,
                         const TableKernel<Dimension>& W,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, typename Dimension::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<Boundary<Dimension>*>& boundaries,
                         const unsigned maxIterations,
                         const double maxFracTol,
                         const double avgFracTol,
                         const RKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dimension, double>& vol,
                         FieldList<Dimension, int>& surfacePoint,
                         FieldList<Dimension, typename Dimension::FacetedVolume>& cells) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Read some state.
  const auto numNodeLists = db.numNodeLists();
  const bool useBounds = (volumeBoundaries.size() == numNodeLists);
  const bool useHoles = (holes.size() == numNodeLists);
  auto pos = db.fluidPosition();
  auto H = db.fluidHfield();
  auto mass = db.fluidMass();
  auto rhof = db.fluidMassDensity();
  auto D = db.solidDamage();

  // Prepare the storage for the point-wise fields.
  auto gradRhof = db.newFluidFieldList(Vector::zero(), "mass density gradient");
  auto deltaCentroid = db.newFluidFieldList(Vector::zero(), "delta centroid");
  auto corrections0 = db.newFluidFieldList(RKCoefficients<Dimension>(), "corrections0");
  auto corrections = db.newFluidFieldList(RKCoefficients<Dimension>(), "corrections");

  // Temporary until we decide to propagate void info to this method.
  auto etaVoidPoints = db.newFluidFieldList(vector<Vector>(), "eta void points");

  // Make a dummy set of cells so we don't ask computeVoronoiVolume to compute the return FacetedVolumes every step.
  FieldList<Dimension, FacetedVolume> dummyCells;
  FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;

  // Make sure the density starts out consistently, and kick-start the volume using m/rho.
  for (auto nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
    const auto n = rhof[nodeListi]->numInternalElements();
    for (auto i = 0U; i != n; ++i) {
      rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
      CHECK(rhof(nodeListi, i) > 0.0);
      vol(nodeListi, i) = mass(nodeListi, i)/rhof(nodeListi, i);
    }
  }

  // // Update the H tensors a bit.
  // iterateIdealH(db, boundaries, W, ASPHSmoothingScale<Dimension>(), 5);

  // Iterate until we converge or max out.
  auto iter = 0u;
  auto avgdelta = 2.0*avgFracTol;
  auto maxdelta = 2.0*maxFracTol;
  while (iter < maxIterations and (avgdelta > avgFracTol or maxdelta > maxFracTol)) {
    std::clock_t t0 = std::clock();
    iter += 1;

    // Remove any old ghost nodes info, and update the mass density
    {
      unsigned nodeListi = 0U;
      for (auto& nodes: db.fluidNodeListPtrs()) {
        nodes->numGhostNodes(0);
        nodes->neighbor().updateNodes();
        if (not rhoConst) for (unsigned i = 0; i != nodes->numInternalNodes(); ++i) rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
      }
      ++nodeListi;
    }

    // Create the new ghost nodes.
    for (auto& bc: boundaries) bc->setAllGhostNodes(db);
    for (auto& bc: boundaries) bc->finalizeGhostBoundary();
    for (auto& nodes: db.fluidNodeListPtrs()) nodes->neighbor().updateNodes();

    // Compute the new connectivity.
    std::clock_t tcm = std::clock();
    db.updateConnectivityMap(false, false, false);
    const auto& cm = db.connectivityMap();
    tcm = std::clock() - tcm;
    { // BLAGO
      double avgneighbors = 0.0;
      unsigned ntot = 0;
      for (auto nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
        const auto n = rhof[nodeListi]->numInternalElements();
        ntot += n;
        for (auto i = 0U; i != n; ++i) avgneighbors += cm.numNeighborsForNode(nodeListi, i);
      }
      ntot = allReduce(ntot, SPHERAL_OP_SUM);
      avgneighbors = allReduce(avgneighbors, SPHERAL_OP_SUM)/ntot;
      if (Process::getRank() == 0) cerr << "Avergage number of neighbors per node: " << avgneighbors << " " << ntot << endl;
    } // BLAGO

    // Compute the new volumes and centroids (note this uses the old rho gradient, not quite right,
    // but expedient/efficient).
    std::clock_t tvoro = std::clock();
    computeVoronoiVolume(pos, H, cm, D, volumeBoundaries, holes, boundaries,
                         FieldList<Dimension, typename Dimension::Scalar>(),  // no weights
                         surfacePoint, vol, deltaCentroid, etaVoidPoints, dummyCells, cellFaceFlags);
    tvoro = std::clock() - tvoro;
     
    // Apply boundary conditions.
    for (auto& bc: boundaries) {
      bc->applyFieldListGhostBoundary(vol);
      bc->applyFieldListGhostBoundary(rhof);
    }
    for (auto& bc: boundaries) bc->finalizeGhostBoundary();

    // If the density is constant we can entirely skip finding the gradient.
    if (not rhoConst) {

      // If the user provided a gradrho method, we can use it.  Otherwise we need to numerically evaluate
      // the density gradient.
      if (useGradRhoFunc) {
        for (unsigned nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
          const auto n = rhof[nodeListi]->numInternalElements();
          for (auto i = 0U; i != n; ++i) {
            gradRhof(nodeListi, i) = gradrhofunc(pos(nodeListi, i));
          }
        }

      } else {
        // Use RK to numerically compute the new mass density gradient.
        ReproducingKernel<Dimension> WR(W, correctionOrder);
        WR.computeCorrections(cm, vol, pos, H, false, corrections0, corrections);
        gradRhof.assignFields(gradientRK(rhof, pos, vol, H, cm, WR, corrections));
      }
    }
     
    // Displace the points and update point masses.
    avgdelta = 0.0;
    maxdelta = 0.0;
    for (unsigned nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = rhof[nodeListi]->numInternalElements();
      for (auto i = 0U; i != n; ++i) {
        auto delta = centroidFrac * deltaCentroid(nodeListi, i);
        if (useBounds) {
          while (not volumeBoundaries[nodeListi].contains(pos(nodeListi, i) + delta, false)) delta *= 0.9;
          if (useHoles) {
            for (const auto& hole: holes[nodeListi]) {
              while (hole.contains(pos(nodeListi, i) + delta, false)) delta *= 0.9;
            }
          }
        }
        const auto deltai = delta.magnitude()/(Dimension::nDim/H(nodeListi, i).Trace());
        avgdelta += deltai;
        maxdelta = std::max(maxdelta, deltai);
        // if (vol(nodeListi, i) > 0.0) H(nodeListi, i) = SymTensor::one() / std::min(hmax, 2.0*Dimension::rootnu(vol(nodeListi, i)));  
// Not correct, but hopefully good enough for our iterative Voronoi purposes.
        pos(nodeListi, i) += delta;
        rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
        if (vol(nodeListi, i) > 0.0) mass(nodeListi, i) = rhof(nodeListi,i)*vol(nodeListi,i);
      }
    }
    avgdelta = (allReduce(avgdelta, SPHERAL_OP_SUM)/
                allReduce(db.numInternalNodes(), SPHERAL_OP_SUM));
    maxdelta = allReduce(maxdelta, SPHERAL_OP_MAX);
    if (Process::getRank() == 0) cerr << "centroidalRelaxNodes iteration " << iter
                                      << ", avg delta frac " << avgdelta 
                                      << ", max delta frac " << maxdelta 
                                      << ", required " << ((std::clock() - t0)/CLOCKS_PER_SEC) 
                                      << " seconds (" << (tvoro/CLOCKS_PER_SEC) << " in computeVoronoiVolume, "
                                      << (tcm/CLOCKS_PER_SEC) << " in ConnectivityMap)." << endl;
        
    // // Update the H tensors a bit.
    // iterateIdealH(db, boundaries, W, ASPHSmoothingScale<Dimension>(), 2);
  }

  // If requested to return the FacetedVolumes, make one last call to fill 'em in.
  if (cells.size() > 0) {
    const auto& cm = db.connectivityMap();
    computeVoronoiVolume(pos, H, cm, D, volumeBoundaries, holes, boundaries,
                         FieldList<Dimension, typename Dimension::Scalar>(),  // no weights
                         surfacePoint, vol, deltaCentroid, etaVoidPoints, cells, cellFaceFlags);
  }

  // Return how many iterations we actually took.
  return iter;
}

}
