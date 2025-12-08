//---------------------------------Spheral++----------------------------------//
// Hydro -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "RK/computeVoronoiVolume.hh"
#include "RK/computeHullVolumes.hh"
#include "RK/computeRKSumVolume.hh"
#include "RK/computeHVolumes.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/correctSPHSumMassDensity.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Strength/SolidFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

#include "FVCRKHydroBase.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <tuple>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::tuple;
using std::make_tuple;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
FVCRKHydroBase<Dimension>::
FVCRKHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
               DataBase<Dimension>& dataBase,
               ArtificialViscosity<Dimension>& Q,
               const RKOrder order,
               const double filter,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool compatibleEnergyEvolution,
               const bool evolveTotalEnergy,
               const bool XSPH,
               const MassDensityType densityUpdate,
               const HEvolutionType HUpdate,
               const double epsTensile,
               const double nTensile):
  CRKSPHHydroBase<Dimension>(smoothingScaleMethod,
                             dataBase,
                             Q,
                             order,
                             filter,
                             cfl,
                             useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution,
                             evolveTotalEnergy,
                             XSPH,
                             densityUpdate,
                             HUpdate,
                             epsTensile,
                             nTensile) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
FVCRKHydroBase<Dimension>::
~FVCRKHydroBase() {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVCRKHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

//   // Get the ArtificialViscosity.
//   auto& Q = this->artificialViscosity();

//   // The kernels and such.
//   const auto& W = this->kernel();
//   const auto& WQ = this->PiKernel();

//   // A few useful constants we'll use in the following loop.
//   const double tiny = 1.0e-30;
//   const auto  compatibleEnergy = this->compatibleEnergyEvolution();
//   const auto  evolveTotalEnergy = this->evolveTotalEnergy();
//   const auto  XSPH = this->XSPH();
//   const auto& smoothingScaleMethod = this->smoothingScaleMethod();

//   // The connectivity.
//   const auto& connectivityMap = dataBase.connectivityMap();
//   const auto& nodeLists = connectivityMap.nodeLists();
//   const auto  numNodeLists = nodeLists.size();
//   const auto  order = this->correctionOrder();

//   // Get the state and derivative FieldLists.
//   // State FieldLists.
//   const auto mass = state.fields(HydroFieldNames::mass, 0.0);
//   const auto volume = state.fields(HydroFieldNames::volume, 0.0);
//   const auto position = state.fields(HydroFieldNames::position, Vector::zero());
//   const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
//   const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
//   const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
//   const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
//   const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
//   const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
//   const auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
//   const auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero());
//   const auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero());
//   const auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero());
//   const auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero());
//   const auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero());
//   const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
//   const auto cells = state.fields(HydroFieldNames::cells, FacetedVolume());
//   CHECK(mass.size() == numNodeLists);
//   CHECK(position.size() == numNodeLists);
//   CHECK(velocity.size() == numNodeLists);
//   CHECK(massDensity.size() == numNodeLists);
//   CHECK(specificThermalEnergy.size() == numNodeLists);
//   CHECK(H.size() == numNodeLists);
//   CHECK(pressure.size() == numNodeLists);
//   CHECK(soundSpeed.size() == numNodeLists);
//   CHECK(A.size() == numNodeLists);
//   CHECK(B.size() == numNodeLists or order == RKOrder::ZerothOrder);
//   CHECK(C.size() == numNodeLists or order != RKOrder::QuadraticOrder);
//   CHECK(gradA.size() == numNodeLists);
//   CHECK(gradB.size() == numNodeLists or order == RKOrder::ZerothOrder);
//   CHECK(gradC.size() == numNodeLists or order != RKOrder::QuadraticOrder);
//   CHECK(surfacePoint.size() == numNodeLists);
//   CHECK(cells.size() == numNodeLists);

//   // Derivative FieldLists.
//   auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero());
//   auto DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
//   auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
//   auto DepsDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
//   auto DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero());
//   auto localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
//   auto DHDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero());
//   auto Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero());
//   auto maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
//   auto effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
//   auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
//   auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
//   auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
//   auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
//   auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero());
//   auto gradRho = derivatives.fields(HydroFieldNames::massDensityGradient, Vector::zero());
//   CHECK(DxDt.size() == numNodeLists);
//   CHECK(DrhoDt.size() == numNodeLists);
//   CHECK(DvDt.size() == numNodeLists);
//   CHECK(DepsDt.size() == numNodeLists);
//   CHECK(DvDx.size() == numNodeLists);
//   CHECK(localDvDx.size() == numNodeLists);
//   CHECK(DHDt.size() == numNodeLists);
//   CHECK(Hideal.size() == numNodeLists);
//   CHECK(maxViscousPressure.size() == numNodeLists);
//   CHECK(effViscousPressure.size() == numNodeLists);
//   CHECK(viscousWork.size() == numNodeLists);
//   CHECK(pairAccelerations.size() == numNodeLists);
//   CHECK(XSPHDeltaV.size() == numNodeLists);
//   CHECK(weightedNeighborSum.size() == numNodeLists);
//   CHECK(massSecondMoment.size() == numNodeLists);
//   CHECK(gradRho.size() == numNodeLists);

//   // Size up the pair-wise accelerations before we start.
//   if (compatibleEnergy) {
//     auto nodeListi = 0;
//     for (auto itr = dataBase.fluidNodeListBegin();
//          itr != dataBase.fluidNodeListEnd();
//          ++itr, ++nodeListi) {
//       for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
//         const size_t n = connectivityMap.numNeighborsForNode(*itr, i);
//         pairAccelerations(nodeListi, i).reserve(n);
//       }
//     }
//   }

//   // Some scratch variables.
//   Scalar Ai, Aj;
//   Vector gradAi, gradAj, forceij, forceji;
//   Vector Bi = Vector::zero(), Bj = Vector::zero();
//   Tensor Ci = Tensor::zero(), Cj = Tensor::zero();
//   Tensor gradBi = Tensor::zero(), gradBj = Tensor::zero();
//   ThirdRankTensor gradCi = ThirdRankTensor::zero(), gradCj = ThirdRankTensor::zero();
//   Scalar gWi, gWj, Wi, Wj, gW0i, gW0j, W0i, W0j;
//   Vector gradWi, gradWj, gradW0i, gradW0j, gradWij;
//   Vector deltagrad;

//   // Prepare the node coupling.
//   const NodeCoupling couple;

//   // Start our big loop over all FluidNodeLists.
//   size_t nodeListi = 0;
//   for (auto itr = dataBase.fluidNodeListBegin();
//        itr != dataBase.fluidNodeListEnd();
//        ++itr, ++nodeListi) {
//     const auto& nodeList = **itr;
//     const auto firstGhostNodei = nodeList.firstGhostNode();
//     const auto hmin = nodeList.hmin();
//     const auto hmax = nodeList.hmax();
//     const auto hminratio = nodeList.hminratio();
//     const auto maxNumNeighbors = nodeList.maxNumNeighbors();
//     const auto nPerh = nodeList.nodesPerSmoothingScale();

//     // Get the work field for this NodeList.
//     auto& workFieldi = nodeList.work();

//     // Iterate over the internal nodes in this NodeList.
//     const auto ni = connectivityMap.numNodes(nodeListi);
// #pragma omp parallel for                                  \
//   firstprivate(Ai, Aj,                                    \
//                gradAi, gradAj, forceij, forceji,          \
//                Bi, Bj,                                    \
//                Ci, Cj,                                    \
//                gradBi, gradBj,                            \
//                gradCi, gradCj,                            \
//                gWi, gWj, Wi, Wj, gW0i, gW0j, W0i, W0j,    \
//                gradWi, gradWj, gradW0i, gradW0j, gradWij, \
//                deltagrad,                                 \
//                couple)
//     for (auto iItr = 0; iItr < ni; ++iItr) {
//       const auto i = connectivityMap.ithNode(nodeListi, iItr);

//       // const bool barf = false; // ((nodeListi == 0 and i >= 98) or (nodeListi == 1 and i <= 1));
//       // if (barf) printf("  --> (%d, %d) :\n", nodeListi, i);

//       // Prepare to accumulate the time.
//       const auto start = Timing::currentTime();

//       // Get the state for node i.
//       const auto& ri = position(nodeListi, i);
//       const auto mi = mass(nodeListi, i);
//       const auto& vi = velocity(nodeListi, i);
//       const auto rhoi = massDensity(nodeListi, i);
//       const auto epsi = specificThermalEnergy(nodeListi, i);
//       const auto Pi = pressure(nodeListi, i);
//       const auto& Hi = H(nodeListi, i);
//       const auto ci = soundSpeed(nodeListi, i);
//       const auto  surfi = surfacePoint(nodeListi, i);
//       Ai = A(nodeListi, i);
//       gradAi = gradA(nodeListi, i);
//       if (order != RKOrder::ZerothOrder) {
//         Bi = B(nodeListi, i);
//         gradBi = gradB(nodeListi, i);
//       }
//       if (order == RKOrder::QuadraticOrder) {
//         Ci = C(nodeListi, i);
//         gradCi = gradC(nodeListi, i);
//       }
//       const auto Hdeti = Hi.Determinant();
//       const auto weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
//       CHECK2(mi > 0.0, i << " " << mi);
//       CHECK2(rhoi > 0.0, i << " " << rhoi);
//       // CHECK2(Ai > 0.0, i << " " << Ai);
//       CHECK2(Hdeti > 0.0, i << " " << Hdeti);
//       CHECK2(weighti > 0.0, i << " " << weighti);

//       auto& DxDti = DxDt(nodeListi, i);
//       auto& DrhoDti = DrhoDt(nodeListi, i);
//       auto& DvDti = DvDt(nodeListi, i);
//       auto& DepsDti = DepsDt(nodeListi, i);
//       auto& DvDxi = DvDx(nodeListi, i);
//       auto& localDvDxi = localDvDx(nodeListi, i);
//       auto& DHDti = DHDt(nodeListi, i);
//       auto& Hideali = Hideal(nodeListi, i);
//       auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
//       auto& effViscousPressurei = effViscousPressure(nodeListi, i);
//       auto& viscousWorki = viscousWork(nodeListi, i);
//       auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
//       auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
//       auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
//       auto& massSecondMomenti = massSecondMoment(nodeListi, i);
//       auto& gradRhoi = gradRho(nodeListi, i);
//       auto& worki = workFieldi(i);

//       // Get the connectivity info for this node.
//       const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

//       // Walk all NodeLists
//       for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {

//         // Connectivity of this node with this NodeList.  We only need to proceed if
//         // there are some nodes in this list.
//         const auto& connectivity = fullConnectivity[nodeListj];
//         if (connectivity.size() > 0) {

//           // Loop over the neighbors.
// #if defined __INTEL_COMPILER
// #pragma vector always
// #endif
//           for (auto jItr = connectivity.begin();
//                jItr != connectivity.end();
//                ++jItr) {
//             const int j = *jItr;

//             // Get the state for node j
//             const auto& rj = position(nodeListj, j);
//             const auto  mj = mass(nodeListj, j);
//             const auto& vj = velocity(nodeListj, j);
//             const auto  rhoj = massDensity(nodeListj, j);
//             const auto  epsj = specificThermalEnergy(nodeListj, j);
//             const auto  Pj = pressure(nodeListj, j);
//             const auto& Hj = H(nodeListj, j);
//             const auto  cj = soundSpeed(nodeListj, j);
//             const auto  surfj = surfacePoint(nodeListj, j);
//             Aj = A(nodeListj, j);
//             gradAj = gradA(nodeListj, j);
//             if (order != RKOrder::ZerothOrder) {
//               Bj = B(nodeListj, j);
//               gradBj = gradB(nodeListj, j);
//             }
//             if (order == RKOrder::QuadraticOrder) {
//               Cj = C(nodeListj, j);
//               gradCj = gradC(nodeListj, j);
//             }
//             const auto Hdetj = Hj.Determinant();
//             const auto weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
//             CHECK(mj > 0.0);
//             CHECK(rhoj > 0.0);
//             CHECK(Hdetj > 0.0);
//             CHECK(weightj > 0.0);

//             // // Find the effective weights of i->j and j->i.
//             // // const auto wi = 2.0*weighti*weightj/(weighti + weightj);
//             // const auto wij = 0.5*(weighti + weightj);

//             // // Are both (i,j) surface points?
//             // // Note we are supposed to have trimmed the topology before this point, so either j
//             // // is an internal point of nodeListi or both (i,j) are surface points.
//             // const bool surfTestij = ((surfi && (1 << (nodeListj + 1)) > 0) and 
//             //                          (surfj && (1 << (nodeListi + 1)) > 0));
//             // CHECK2((nodeListj == nodeListi) or surfTestij or
//             //        (surfi && (1 << (nodeListj + 1)) > 0) or
//             //        (surfj && (1 << (nodeListi + 1)) > 0), "(" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") : " << surfi << " " << surfj << " " << surfTestij);

//             // Node displacement.
//             const auto rij = ri - rj;
//             const auto etai = Hi*rij;
//             const auto etaj = Hj*rij;
//             const auto etaMagi = etai.magnitude();
//             const auto etaMagj = etaj.magnitude();
//             CHECK(etaMagi >= 0.0);
//             CHECK(etaMagj >= 0.0);
//             const auto vij = vi - vj;

//             // Symmetrized kernel weight and gradient.
//             CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, order,  rij,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi);
//             CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, order, -rij, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj);
//             deltagrad = gradWj - gradWi;
//             const auto gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);

//             // Zero'th and second moment of the node distribution -- used for the
//             // ideal H calculation.
//             const auto fweightij = nodeListi == nodeListj ? 1.0 : mj*rhoi/(mi*rhoj);
//             const auto rij2 = rij.magnitude2();
//             const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
//             weightedNeighborSumi +=     fweightij*std::abs(gWi);
//             massSecondMomenti +=        fweightij*gradWSPHi.magnitude2()*thpt;

//             // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
//             const auto QPiij = Q.Piij(nodeListi, i, nodeListj, j,
//                                       ri, etai, vi, rhoi, ci, Hi,
//                                       rj, etaj, vj, rhoj, cj, Hj);
//             const auto Qaccij = (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second).dot(deltagrad);
//             const auto workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);
//             const auto Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
//             maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                        // We need tighter timestep controls on the Q with CRK
//             effViscousPressurei += weightj * Qi * Wj;
//             viscousWorki += 0.5*weighti*weightj/mi*workQi;

//             // Velocity gradient.
//             DvDxi -= weightj*vij.dyad(gradWj);
//             if (nodeListi == nodeListj) localDvDxi -= weightj*vij.dyad(gradWj);

//             // Mass density gradient.
//             gradRhoi += weightj*(rhoj - rhoi)*gradWj;

//             // The force between the points depends on the surface test
//             // Momentum
//             forceij = (false ? //(surfi != 0 and surfj != 0) ?
//                        weighti*weightj*((Pi + Pj)*gradWj + rhoj*rhoj*QPiij.second.dot(gradWj)) :                   // surface
//                        0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij));                // Type III CRK interpoint force.
//             DvDti -= forceij/mi;

//             // if (barf) {
//             //   printf(" (%d, %d): ", nodeListj, j);
//             //   cout << "  " << DvDti << " " << -0.5*wij*wij*((Pi + Pj)*deltagrad + Qaccij)/mi;
//             //   if (surfTestij) {
//             //     cout << "  <--- surface\n";
//             //   } else {
//             //     cout << "\n";
//             //   }
//             // }

//             if (compatibleEnergy) pairAccelerationsi.push_back(-forceij/mi);

//             // Energy
//             DepsDti += (false ? //(surfi != 0 and surfj != 0) ?
//                         weighti*weightj*(Pj*vij.dot(gradWj) - rhoj*rhoj*QPiij.second.dot(vij).dot(gradWj))/mi :         // surface
//                         0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + workQi)/mi);         // CRK

//             // Estimate of delta v (for XSPH).
//             if (XSPH and (nodeListi == nodeListj)) XSPHDeltaVi -= weightj*Wj*vij;
//           }
//         }
//       }

//       const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
//       CHECK(not compatibleEnergy or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
//             (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
//             (pairAccelerationsi.size() == 2*numNeighborsi));

//       // if (barf) printf("\n");

//       // // For a surface point, add the RK thermal energy evolution.
//       // // DepsDti -= Pi/rhoi*DvDxi.Trace();
//       // if (surfacePoint(nodeListi, i) > 1) DepsDti -= Pi/rhoi*DvDxi.Trace();

//       // Get the time for pairwise interactions.
//       const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/max(1, numNeighborsi);

//       // Time evolution of the mass density.
//       DrhoDti = -rhoi*DvDxi.Trace();

//       // If needed finish the total energy derivative.
//       if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

//       // Complete the moments of the node distribution for use in the ideal H calculation.
//       weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
//       massSecondMomenti /= Hdeti*Hdeti;

//       // Determine the position evolution, based on whether we're doing XSPH or not.
//       if (XSPH) {
//         DxDti = vi + XSPHDeltaVi;
//       } else {
//         DxDti = vi;
//       }

//       // The H tensor evolution.
//       DHDti = smoothingScaleMethod.smoothingScaleDerivative(Hi,
//                                                             ri,
//                                                             DvDxi,
//                                                             hmin,
//                                                             hmax,
//                                                             hminratio,
//                                                             nPerh);
//       Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
//                                                        ri,
//                                                        weightedNeighborSumi,
//                                                        massSecondMomenti,
//                                                        W,
//                                                        hmin,
//                                                        hmax,
//                                                        hminratio,
//                                                        nPerh,
//                                                        connectivityMap,
//                                                        nodeListi,
//                                                        i);

//       // Increment the work for i.
//       worki += Timing::difference(start, Timing::currentTime());
//     }
//   }
}

}
