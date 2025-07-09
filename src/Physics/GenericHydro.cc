//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "GenericHydro.hh"

#include "Physics.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Strength/SolidFieldNames.hh"

#include <limits>
#include <algorithm>
#include <sstream>

using std::vector;
using std::string;
using std::to_string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::make_pair;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Helper method to compute the magnitude of the shear term in the given
// dvdx tensor.
//------------------------------------------------------------------------------
inline
double
computeShearMagnitude(const Dim<1>::Tensor& /*dvdx*/) {
  return 0.0;
}

inline
double
computeShearMagnitude(const Dim<2>::Tensor& dvdx) {
  return std::abs(dvdx(1,0) - dvdx(0,1));
}

inline
double
computeShearMagnitude(const Dim<3>::Tensor& dvdx) {
  return sqrt(FastMath::square(dvdx(2,1) - dvdx(1,2)) +
              FastMath::square(dvdx(2,0) - dvdx(0,2)) +
              FastMath::square(dvdx(1,0) - dvdx(0,1)));
}

// Provide to_string for Spheral Vector
template<typename Vector>
std::string
vec_to_string(const Vector& vec) {
  std::ostringstream oss;
  oss << vec << std::endl;
  return oss.str();
}

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
GenericHydro<Dimension>::
GenericHydro(ArtificialViscosityHandle<Dimension>& Q,
             const double cfl,
             const bool useVelocityMagnitudeForDt):
  Physics<Dimension>(),
  mArtificialViscosity(Q),
  mCfl(cfl),
  mUseVelocityMagnitudeForDt(useVelocityMagnitudeForDt),
  mMinMasterNeighbor(INT_MAX),
  mMaxMasterNeighbor(0),
  mSumMasterNeighbor(0),
  mMinCoarseNeighbor(INT_MAX),
  mMaxCoarseNeighbor(0),
  mSumCoarseNeighbor(0),
  mMinRefineNeighbor(INT_MAX),
  mMaxRefineNeighbor(0),
  mSumRefineNeighbor(0),
  mMinActualNeighbor(INT_MAX),
  mMaxActualNeighbor(0),
  mSumActualNeighbor(0),
  mNormMasterNeighbor(0),
  mNormCoarseNeighbor(0),
  mNormRefineNeighbor(0),
  mNormActualNeighbor(0),
  mDTrank(0u),
  mDTNodeList(0u),
  mDTnode(0u),
  mDTreason() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GenericHydro<Dimension>::~GenericHydro() {
}

//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename GenericHydro<Dimension>::TimeStepType
GenericHydro<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   typename Dimension::Scalar /*currentTime*/) const {

  const double tiny = std::numeric_limits<double>::epsilon();

  // Get some useful fluid variables from the DataBase.
  const auto  mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  const auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
                                                         this->requireOverlapConnectivity(),
                                                         this->requireIntersectionConnectivity());
  const auto  numNodeLists = connectivityMap.nodeLists().size();

  // Stuff from the artificial viscosity
  //const auto& Q = this->artificialViscosity();
  //const auto  Cl = Q.Cl();
  //const auto  Cq = Q.Cq();
  //const auto  Qeps2 = Q.epsilon2();

  // Check for deviatoric stress.
  const auto haveDS = state.fieldNameRegistered(SolidFieldNames::deviatoricStress);
  FieldList<Dimension, SymTensor> S;
  if (haveDS) S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);

  // Check if the longitudinal sound speed is registered.
#ifdef CXXONLY
  const auto haveLongCs = false; // ignore longitudinal sound speed for CInterface
#else
  const auto haveLongCs = state.fieldNameRegistered(SolidFieldNames::longitudinalSoundSpeed);
#endif
  FieldList<Dimension, Scalar> csl;
  if (haveLongCs) csl = state.fields(SolidFieldNames::longitudinalSoundSpeed, 0.0);

  // Initialize the return value to some impossibly high value.
  auto minDt = make_pair(std::numeric_limits<double>::max(), string());

  // Define a function for computing the velocity divergence (different for curvilinear coordinates)
  auto Fdiv = +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi.Trace(); };
  if (GeometryRegistrar::coords() == CoordinateType::Spherical) {
    Fdiv = +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi[0] + 2.0*veli[0]*safeInv(posi[0]); };
  } else if (GeometryRegistrar::coords() == CoordinateType::RZ) {
    Fdiv = +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi.Trace() + veli[1]*safeInv(posi[1]); };
  }

  // Loop over every fluid node.
  // #pragma omp declare reduction (MINPAIR : pair<double,string> : omp_out = (omp_out.first < omp_in.first ? omp_out : omp_in)) initializer(omp_priv = pair<double,string>(std::numeric_limits<double>::max(), string("null")))
  // #pragma omp parallel for reduction(MINPAIR:minDt) collapse(2)
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& fluidNodeList = **(dataBase.fluidNodeListBegin() + nodeListi);
    const auto nPerh = fluidNodeList.nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);

    // Check if we have a longitudinal sound speed for this material.
    const bool useCsl = haveLongCs and csl.haveNodeList(fluidNodeList);
    const Field<Dimension, Scalar>* cslptr = nullptr;
    if (useCsl) cslptr = *csl.fieldForNodeList(fluidNodeList);

    // Check if we have a deviatoric stress for this material.
    const bool useS = haveDS and S.haveNodeList(fluidNodeList);
    const Field<Dimension, SymTensor>* Sptr = nullptr;
    if (useS) Sptr = *S.fieldForNodeList(fluidNodeList);

    // Walk all the nodes in this FluidNodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
    const auto rank = Process::getRank();
    // #pragma omp parallel for reduction(MINPAIR:minDt)
#pragma omp parallel
    {
      auto minDt_local = minDt;
      auto DTrank_local = mDTrank;
      auto DTNodeList_local = mDTNodeList;
      auto DTnode_local = mDTnode;
      auto DTreason_local = mDTreason;
#pragma omp for
      for (auto k = 0; k < ni; ++k) {
        const auto i = connectivityMap.ithNode(nodeListi, k);

        // If this node is masked, don't worry about it.
        if (mask(nodeListi, i) == 1) {

          // Get this nodes minimum characteristic smoothing scale.
          CHECK2(H(nodeListi, i).Determinant() >  0.0,
                 "Bad H tensor : " << H(nodeListi, i) << " : " << fluidNodeList.name() << " " << i << " " << fluidNodeList.firstGhostNode());
          const auto& Hi = H(nodeListi, i);
          const Scalar nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;
          // const Scalar nodeScalei = 1.0/Dimension::rootnu(Hi.Determinant());
          //     const Scalar nodeScalei = nodeExtent(nodeListi, i).minElement()/kernelExtent;

          // Sound speed limit.
          const auto csi = cs(nodeListi, i);
          const auto csDt = nodeScalei/(csi + tiny);
          if (csDt < minDt_local.first) {
            minDt_local = make_pair(csDt, ("Sound speed limit: dt = " + to_string(csDt) + "\n" +
                                           "                   cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                           "            nodeScale = " + to_string(nodeScalei) + "\n" +
                                           "             material = " + fluidNodeList.name() + "\n" +
                                           "(nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                           "           @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "sound speed";
          }

          // Longitudinal sound speed limit.
          if (useCsl) {
            const auto csDt = nodeScalei/((*cslptr)(i) + tiny);
            if (csDt < minDt_local.first) {
              minDt_local = make_pair(csDt, ("Longitudinal sound speed limit: dt = " + to_string(csDt) + "\n" + 
                                             "                                cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                             "                               csl = " + to_string((*cslptr)(i)) + "\n" +
                                             "                         nodeScale = " + to_string(nodeScalei) + "\n" +
                                             "                          material = " + fluidNodeList.name() + "\n" +
                                             "             (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                             "                        @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "longitudinal sound speed";
            }
          }

          // Deviatoric stress limit.
          const auto rhoi = rho(nodeListi, i);
          CHECK(rhoi > 0.0);
          if (useS) {
            const auto csS = sqrt((*Sptr)(i).eigenValues().maxAbsElement()/rhoi);
            const auto SDt = nodeScalei/(csS + tiny);
            if (SDt < minDt_local.first) {
              minDt_local = make_pair(SDt, ("Deviatoric stress effective sound speed limit: dt = " + to_string(SDt) + "\n" +
                                            "                                               cs = " + to_string(cs(nodeListi, i)) + "\n" + 
                                            "                                              csS = " + to_string(csS) + "\n" +
                                            "                                              rho = " + to_string(rhoi) + "\n" +
                                            "                                        nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                                         material = " + fluidNodeList.name() + "\n" +
                                            "                            (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                                       @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "deviatoric stress effective sound speed";
            }
          }

          // Artificial viscosity effective sound speed.
          const auto csq = sqrt(maxViscousPressure(nodeListi, i)/rhoi);
          const auto csqDt = nodeScalei/(csq + tiny);
          if (csqDt < minDt_local.first) {
            minDt_local = make_pair(csqDt, ("Artificial viscosity sound speed limit: dt = " + to_string(csqDt) + "\n" + 
                                            "                                        cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                            "                                       csQ = " + to_string(csq) + "\n" +
                                            "                                       rho = " + to_string(rhoi) + "\n" +
                                            "                                 nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                                  material = " + fluidNodeList.name() + "\n" +
                                            "                     (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                                @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "artificial viscosity";
          }

          // Velocity divergence limit.
          const auto divVelocity = Fdiv(DvDx(nodeListi, i), position(nodeListi, i), velocity(nodeListi, i));
          const auto divvDt = 1.0/(std::abs(divVelocity) + tiny);
          if (divvDt < minDt_local.first) {
            minDt_local = make_pair(divvDt, ("Velocity divergence limit: dt = " + to_string(divvDt) + "\n" +
                                             "                 div velocity = " + to_string(divVelocity) + "\n" +
                                             "                     material = " + fluidNodeList.name() + "\n" +
                                             "        (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                             "                   @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "velocity divergence";
          }

          // Maximum velocity difference limit.
          //const auto& xi = position(nodeListi, i);
          const auto& vi = velocity(nodeListi, i);
          const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
          for (auto nodeListj = 0u; nodeListj != numNodeLists; ++nodeListj) {
            const auto& connectivity = fullConnectivity[nodeListj];
            for (auto jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) {
              const auto  j = *jItr;
              //const auto& xj = position(nodeListj, j);
              const auto& vj = velocity(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              //const auto  rhoj = rho(nodeListj, j);
              const Scalar nodeScalej = 1.0/Hj.eigenValues().maxElement()/nPerh;
              // const Scalar vij = std::abs((vj - vi).dot((xi - xj).unitVector()));
              const auto  vij = vi - vj;
              const auto  dtVelDiff = std::min(nodeScalei, nodeScalej)*safeInvVar(vij.magnitude(), 1e-30);
              if (dtVelDiff < minDt_local.first) {
                minDt_local = make_pair(dtVelDiff, ("Pairwise velocity difference limit: dt = " + to_string(dtVelDiff) + "\n" + 
                                                    "                              material = " + fluidNodeList.name() + "\n" +
                                                    "                  (nodeListi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                                    "                  (nodeListj, i, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(rank) + ")\n" +
                                                    "                   @ pos(nodeListi, i) = " + vec_to_string(position(nodeListi, i)) + "\n" +
                                                    "                   @ pos(nodeListj, j) = " + vec_to_string(position(nodeListj, j)) + "\n" +
                                                    "                                   vij = " + to_string(vij.magnitude()) + "\n" +
                                                    "                            nodeScalei = " + to_string(nodeScalei) + "\n" +
                                                    "                            nodeScalej = " + to_string(nodeScalej)));
                DTrank_local = rank;
                DTNodeList_local = nodeListi;
                DTnode_local = i;
                DTreason_local = "pairwise velocity difference";
              }

              // // We also use a pairwise condition modeled on the Monaghan-Gingold viscosity formulation.
              // const auto csj = cs(nodeListj, j);
              // const auto xij = (xi - xj);
              // const auto etai = Hi*xij;
              // const auto etaj = Hj*xij;
              // const auto mui = std::abs(vij.dot(etai)/(etai.magnitude2() + Qeps2));
              // const auto muj = std::abs(vij.dot(etaj)/(etaj.magnitude2() + Qeps2));
              // const auto ei = Cl*csi*mui + Cq*mui*mui;
              // const auto ej = Cl*csj*muj + Cq*muj*muj;
              // const auto csQpair = 2.0*std::min(rhoi*ei, rhoj*ej);
              // const auto dtQpair = std::min(nodeScalei, nodeScalej)/(csQpair + tiny);
              // if (dtQpair < minDt_local.first) {
              //   minDt_local = make_pair(dtQpair, ("                  Pairwise Q limit: dt = " + to_string(dtVelDiff) + "\n" + 
              //                                     "                              material = " + fluidNodeList.name() + "\n" +
              //                                     "                  (nodeListi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
              //                                     "                  (nodeListj, i, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(Process::getRank()) + ")\n" +
              //                                     "                   @ pos(nodeListi, i) = " + vec_to_string(position(nodeListi, i)) + "\n" +
              //                                     "                   @ pos(nodeListj, j) = " + vec_to_string(position(nodeListj, j)) + "\n" +
              //                                     "                               csQpair = " + to_string(csQpair) + "\n"
              //                                     "                                   mui = " + to_string(mui) + "\n"
              //                                     "                                   muj = " + to_string(muj) + "\n"
              //                                     "                                   vij = " + vec_to_string(vij) + "\n" +
              //                                     "                            nodeScalei = " + to_string(nodeScalei) + "\n" +
              //                                     "                            nodeScalej = " + to_string(nodeScalej)));
              // }
            }
          }

          //     // Eigenvalues of the stress-strain tensor.
          //     Vector eigenValues = DvDx(nodeListi, i).Symmetric().eigenValues();
          //     Scalar maxValue = -1.0;
          //     for (int i = 0; i < Dimension::nDim; ++i) {
          //       maxValue = max(maxValue, fabs(eigenValues(i)));
          //     }
          //     const double strainDt = 1.0/(maxValue + tiny);
          //     if (strainDt < minDt) {
          //       minDt = strainDt;
          //       reason = "strain limit";
          //     }

          //     // Limit by the velocity shear.
          //     const Scalar shearVelocity = computeShearMagnitude(DvDx(nodeListi, i));
          //     const double shearDt = 1.0/(shearVelocity + tiny);
          //     if (shearDt < minDt) {
          //       minDt = shearDt;
          //       reason = "velocity shear limit";
          //     }

          // Total acceleration limit.
          const auto vmagi = vi.magnitude();
          const auto dtAcc = 0.1*std::max(nodeScalei/(vmagi + tiny), vmagi/(DvDt(nodeListi, i).magnitude() + tiny));
          if (dtAcc < minDt_local.first) {
            minDt_local = make_pair(dtAcc, ("Total acceleration limit: dt = " + to_string(dtAcc) + "\n" + 
                                            "              |acceleration| = " + to_string(DvDt(nodeListi, i).magnitude()) + "\n" +
                                            "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                    material = " + fluidNodeList.name() + "\n" +
                                            "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                  @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "acceleration";
          }

          // If requested, limit against the absolute velocity.
          if (useVelocityMagnitudeForDt()) {
            const auto velDt = nodeScalei/(velocity(nodeListi, i).magnitude() + 1.0e-10);
            if (velDt < minDt_local.first) {
              minDt_local = make_pair(velDt, ("Velocity magnitude limit: dt = " + to_string(velDt) + "\n" +
                                              "                        |vi| = " + to_string(velocity(nodeListi, i).magnitude()) + "\n" +
                                              "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                              "                    material = " + fluidNodeList.name() + "\n" +
                                              "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                              "                  @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "velocity magnitude";
            }
          }
        }
      }

#pragma omp critical
      {
        if (minDt_local.first < minDt.first) {
          minDt = minDt_local;
          mDTrank = DTrank_local;
          mDTNodeList = DTNodeList_local;
          mDTnode = DTnode_local;
          mDTreason = DTreason_local;
        }
      }
    }
  }

  // Scale by the cfl safety factor.
  minDt.first *= cfl();
  return minDt;
}

}
