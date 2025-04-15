//---------------------------------Spheral++----------------------------------//
// SPHSmoothingScale
//
// Implements the standard SPH scalar smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#include "SmoothingScale/SPHSmoothingScale.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/Timer.hh"

#include <cmath>
#include <vector>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::vector;

namespace {

//------------------------------------------------------------------------------
// Convert a given number of neighbors to the equivalent 1D "radius" in nodes.
//------------------------------------------------------------------------------
template<typename Dimension> inline double equivalentRadius(const double n);

// 1D
template<>
inline double
equivalentRadius<Dim<1> >(const double n) {
  return 0.5*n;
}

// 2D
template<>
inline double
equivalentRadius<Dim<2> >(const double n) {
  return std::sqrt(n/M_PI);
}

// 3D
template<>
inline double
equivalentRadius<Dim<3> >(const double n) {
  return Dim<3>::rootnu(3.0*n/(4.0*M_PI));
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHSmoothingScale<Dimension>::
SPHSmoothingScale(const HEvolutionType HUpdate,
                  const TableKernel<Dimension>& W):
  SmoothingScaleBase<Dimension>(HUpdate),
  mWT(W),
  mZerothMoment(FieldStorageType::CopyFields),
  mFirstMoment(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHSmoothingScale<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Make sure our FieldLists are correctly sized.
  SmoothingScaleBase<Dimension>::initializeProblemStartup(dataBase);
  dataBase.resizeFluidFieldList(mZerothMoment, 0.0, HydroFieldNames::massZerothMoment, false);
  dataBase.resizeFluidFieldList(mFirstMoment, Vector::zero, HydroFieldNames::massFirstMoment, false);
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHSmoothingScale<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  SmoothingScaleBase<Dimension>::registerDerivatives(dataBase, derivs);
  derivs.enroll(mZerothMoment);
  derivs.enroll(mFirstMoment);
}

//------------------------------------------------------------------------------
// Time derivative of the smoothing scale.
// We depend on a previous package evaluating the velocity gradient (DvDx)
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHSmoothingScale<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("SPHSmoothingScaleDerivs");

  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  etaMax = mWT.kernelExtent();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  CHECK(position.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DHDt = derivs.fields(IncrementBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  massZerothMoment = derivs.fields(HydroFieldNames::massZerothMoment, 0.0);
  auto  massFirstMoment = derivs.fields(HydroFieldNames::massFirstMoment, Vector::zero);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(massZerothMoment.size() == numNodeLists);
  CHECK(massFirstMoment.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar mi, mj, ri, rj, mRZi, mRZj, rhoi, rhoj, WSPHi, WSPHj, etaMagi, etaMagj, fweightij, fispherical, fjspherical;
    Vector xij, etai, etaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto massZerothMoment_thread = massZerothMoment.threadCopy(threadStack);
    auto massFirstMoment_thread = massFirstMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      mi = mass(nodeListi, i);
      rhoi = massDensity(nodeListi, i);
      const auto& xi = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);

      auto& massZerothMomenti = massZerothMoment_thread(nodeListi, i);
      auto& massFirstMomenti = massFirstMoment_thread(nodeListi, i);

      // Get the state for node j
      mj = mass(nodeListj, j);
      rhoj = massDensity(nodeListj, j);
      const auto& xj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);

      auto& massZerothMomentj = massZerothMoment_thread(nodeListj, j);
      auto& massFirstMomentj = massFirstMoment_thread(nodeListj, j);

      // Node displacement.
      xij = xi - xj;
      etai = Hi*xij;
      etaj = Hj*xij;
      etaMagi = etai.magnitude();
      etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Compute the node-node weighting
      fweightij = 1.0;
      fispherical = 1.0;
      fjspherical = 1.0;
      if (nodeListi != nodeListj) {
        if (GeometryRegistrar::coords() == CoordinateType::RZ) {
          ri = abs(xi.y());
          rj = abs(xj.y());
          mRZi = mi/(2.0*M_PI*ri);
          mRZj = mj/(2.0*M_PI*rj);
          fweightij = mRZj*rhoi/(mRZi*rhoj);
        } else {
          fweightij = mj*rhoi/(mi*rhoj);
        }
      } else if (GeometryRegistrar::coords() == CoordinateType::Spherical) {
        const auto eii = Hi.xx()*xi.x();
        const auto eji = Hi.xx()*xj.x();
        const auto ejj = Hj.xx()*xj.x();
        const auto eij = Hj.xx()*xi.x();
        fispherical = (eii > etaMax ? 1.0 :
                       eii < eji ? 2.0 :
                       0.0);
        fjspherical = (ejj > etaMax ? 1.0 :
                       ejj < eij ? 2.0 :
                       0.0);
      }

      // Symmetrized kernel weight
      WSPHi = mWT.kernelValueSPH(etaMagi);
      WSPHj = mWT.kernelValueSPH(etaMagj);

      // Moments of the node distribution -- used for the ideal H calculation.
      massZerothMomenti +=     fweightij*WSPHi * fispherical;
      massZerothMomentj += 1.0/fweightij*WSPHj * fjspherical;
      massFirstMomenti -=     fweightij*WSPHi*etai;
      massFirstMomentj += 1.0/fweightij*WSPHj*etaj;
    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region

  // Finish up the derivatives now that we've walked all pairs
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& Hi = H(nodeListi, i);
      const auto& DvDxi = DvDx(nodeListi, i);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      auto& massZerothMomenti = massZerothMoment(nodeListi, i);
      massZerothMomenti = Dimension::rootnu(max(0.0, massZerothMomenti));

      // Time derivative of H
      DHDt(nodeListi, i) = -Hi/(Dimension::nDim)*DvDxi.Trace();

      // Determine the current effective number of nodes per smoothing scale.
      const auto currentNodesPerSmoothingScale = (fuzzyEqual(massZerothMomenti, 0.0) ?  // Is this node isolated (no neighbors)?
                                                  0.5*nPerh :
                                                  mWT.equivalentNodesPerSmoothingScale(massZerothMomenti));
      CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

      // The ratio of the desired to current nodes per smoothing scale.
      const auto s = std::min(4.0, std::max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
      CHECK(s > 0.0);

      // Now determine how to scale the current H to the desired value.
      const auto a = (s < 1.0 ? 
                      0.4*(1.0 + s*s) :
                      0.4*(1.0 + 1.0/(s*s*s)));
      CHECK(1.0 - a + a*s > 0.0);
      const auto hi0 = 1.0/Hi.xx();
      const auto hi1 = std::min(hmax, std::max(hmin, hi0*(1.0 - a + a*s)));
      CHECK(hi1 > 0.0);
      Hideal(nodeListi, i) = 1.0/hi1 * SymTensor::one;
    }
  }
  TIME_END("SPHSmoothingScaleDerivs");
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHSmoothingScale<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  SmoothingScaleBase<Dimension>::dumpState(file, pathName);
  file.write(mZerothMoment, pathName + "/zerothMoment");
  file.write(mFirstMoment, pathName + "/firstMoment");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHSmoothingScale<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  SmoothingScaleBase<Dimension>::restoreState(file, pathName);
  file.read(mZerothMoment, pathName + "/zerothMoment");
  file.read(mFirstMoment, pathName + "/firstMoment");
}

}
