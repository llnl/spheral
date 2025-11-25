
//---------------------------------Spheral++----------------------------------//
// SubPointPressureHourglassControl
//
// Impose additional forces on each point using subdivisions of the Voronoi
// control volume to constrain the unphysical degrees of freedom in our hydro
// discretization and avoid spurious so-called houglass modes.
//----------------------------------------------------------------------------//
#include "VoronoiCells/SubPointPressureHourglassControl.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/CellFaceFlag.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/Timer.hh"
#include "Utilities/range.hh"

#include <limits>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <sstream>

namespace Spheral {

using FastMath::pow2;

using std::vector;
using std::unordered_map;
using std::tuple;
using std::string;
using std::to_string;

namespace {  // anonymous

template<typename Vector>
inline
double
limiter(const Vector& xi,
        const Vector& xn,
        const double  rhoi,
        const Vector& gradRhoi) {
  const auto delta = gradRhoi.dot(xn - xi);
  return (rhoi + delta < 0.0 ?
          abs(rhoi/delta) :
          1.0);
}

//------------------------------------------------------------------------------
// Center of mass (1D)
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
centerOfMass(const Dim<1>::FacetedVolume& celli,
             const Dim<1>::Vector& xi,
             const Dim<1>::Scalar rhoi,
             Dim<1>::Vector gradRhoi) {
  using Vector = Dim<1>::Vector;
  const auto& x1 = celli.xmin();
  const auto& x2 = celli.xmax();
  gradRhoi *= limiter(xi, x1, rhoi, gradRhoi);
  gradRhoi *= limiter(xi, x2, rhoi, gradRhoi);
  CHECK(rhoi + gradRhoi.dot(x1 - xi) >= 0.0);
  CHECK(rhoi + gradRhoi.dot(x2 - xi) >= 0.0);
  const auto c = xi.x();
  const auto ab = x1.x() + x2.x() - 2.0*c;
  const auto m = gradRhoi.x();
  const Vector result(c + ab*(2.0*m*ab + 3.0*rhoi)/(3.0*(m*ab + 2.0*rhoi)));
  ENSURE2(result.x() >= x1.x() and result.x() <= x2.x(), result << " not in [" << x1.x() << " " << x2.x() << "] : " << gradRhoi << " " << (rhoi + gradRhoi.dot(x1 - xi)) << " " << (rhoi + gradRhoi.dot(x2 - xi)));
  return result;
}

//------------------------------------------------------------------------------
// Center of mass (2D)
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
centerOfMass(const Dim<2>::FacetedVolume& celli,
             const Dim<2>::Vector& xi,
             const Dim<2>::Scalar rhoi,
             Dim<2>::Vector gradRhoi) {
  using Vector = Dim<2>::Vector;
  return Vector();
}

//------------------------------------------------------------------------------
// Center of mass (3D)
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
centerOfMass(const Dim<3>::FacetedVolume& celli,
             const Dim<3>::Vector& xi,
             const Dim<3>::Scalar rhoi,
             Dim<3>::Vector gradRhoi) {
  using Vector = Dim<3>::Vector;
  return Vector();
}

//------------------------------------------------------------------------------
// Compute the internal acceleration due to a single facet (1D)
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
subCellAcceleration(const Dim<1>::FacetedVolume& celli,
                    const int cellFace,
                    const Dim<1>::Vector& comi,
                    const Dim<1>::Vector& xi,
                    const Dim<1>::Scalar  Pi) {
  using Vector = Dim<1>::Vector;
  REQUIRE(cellFace == 0 or cellFace == 1);

  const auto vertx = cellFace == 0 ? celli.xmin().x() : celli.xmax().x();
  const auto fA = cellFace == 0 ? Vector(1.0) : Vector(-1.0);   // Inward pointing normal since we want -\grad P
  // const auto dA = (comi - vert).unitVector();   // Inward pointing normal since we want -\grad P
  const auto dA0 = vertx - comi.x();
  const auto dA1 = vertx - xi.x();
  const auto Psub = abs(Pi) * (dA0*safeInvVar(dA1) - 1.0);
  // const auto Psub = abs(Pi) * max(-1.0, min(1.0, 1.0 - dA1*safeInv(dA0)));
  // const auto Psub = abs(Pi * (vert.x() - comi.x())/(vert.x() - xi.x()));
  return Psub*fA;

  // // Define a function to increment the acceleration for each subcell
  // auto asub = [&](const Vector& vert) -> Vector {
  //               const auto dA = (comi - vert).unitVector();
  //               const auto Psub = abs(Pi * (vert.x() - comi.x())/(vert.x() - xi.x()));
  //               return Psub*dA;
  //             };

  // // Now we can sum up finite volume contribution to the acceleration for each subvolume
  // const auto Vi = celli.volume();
  // CHECK(Vi > 0.0);
  // return (asub(celli.xmin()) + asub(celli.xmax()))/(rhoi*Vi);
}

//------------------------------------------------------------------------------
// Compute the internal acceleration (2D)
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
subCellAcceleration(const Dim<2>::FacetedVolume& celli,
                    const int cellFace,
                    const Dim<2>::Vector& comi,
                    const Dim<2>::Vector& xi,
                    const Dim<2>::Scalar  Pi) {
  const auto& facets = celli.facets();
  REQUIRE(size_t(cellFace) < facets.size());
  // auto  Atot = 0.0;
  // for (const auto& f: facets) Atot += f.area();
  const auto& f = facets[cellFace];
  const auto  nA = -f.normal();  // Inward pointing area normal (has magnitude of facet area)
  // const auto Aref = Atot/6u;
  // const auto  nAref = Aref * nA.unitVector();
  const auto& p0 = f.point1();
  const auto  dV0 = (comi - p0).dot(nA);    // 2* actually
  const auto  dV1 = (xi   - p0).dot(nA);    // 2* actually
  const auto  Psub = abs(Pi) * FastMath::pow4(dV0*safeInvVar(dV1) - 1.0);
  // const auto  dV0 = (comi - p0).dot(nA);    // 2* actually
  // const auto  dV1 = (xi   - p0).dot(nA);    // 2* actually
  // const auto  Psub = abs(Pi) * (dV0*safeInvVar(dV1) - 1.0);
  return Psub*nA;

  // const auto comi = celli.centroid();

  // // Define a function to increment the acceleration for each subcell
  // auto asub = [&](const Vector& v1, const Vector& v2) -> Vector {
  //               const auto v12 = v2 - v1;
  //               const Vector dA(-v12.y(), v12.x());
  //               const auto Psub = abs(Pi * ((v1 - comi).cross(v2 - comi)).z()*safeInv(((v1 - xi).cross(v2 - xi)).z()));
  //               return Psub*dA;
  //             };

  // // Now we can sum up finite volume contribution to the acceleration for each subvolume.
  // Vector result;
  // const auto& facets = celli.facets();
  // for (auto& f: facets) result += asub(f.point1(), f.point2());
  // const auto Vi = celli.volume();
  // CHECK(Vi > 0.0);
  // result /= rhoi*Vi;
  // return result;
}

//------------------------------------------------------------------------------
// Compute the internal acceleration (3D)
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
subCellAcceleration(const Dim<3>::FacetedVolume& celli,
                    const int cellFace,
                    const Dim<3>::Vector& comi,
                    const Dim<3>::Vector& xi,
                    const Dim<3>::Scalar  Pi) {
  const auto& facets = celli.facets();
  REQUIRE(size_t(cellFace) < facets.size());
  const auto& f = facets[cellFace];
  const auto  fA = f.area();
  const auto  nhat = -f.normal();  // Inward pointing
  const auto  p0 = f.point(0u);
  const auto  dV0 = (comi - p0).dot(nhat) * fA;    // 3* actually
  const auto  dV1 = (xi   - p0).dot(nhat) * fA;    // 3* actually
  const auto  Psub = abs(Pi) * (dV0*safeInvVar(dV1) - 1.0);
  return Psub*fA * nhat;
}

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SubPointPressureHourglassControl<Dimension>::
SubPointPressureHourglassControl(const Scalar fHG,
                                 const Scalar xfilter):
  mfHG(fHG),
  mxfilter(xfilter),
  mDvDt(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SubPointPressureHourglassControl<Dimension>::
~SubPointPressureHourglassControl() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  mDvDt = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::ahgAcceleration);
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
}

//------------------------------------------------------------------------------
// Register the derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  // derivs.enroll(mDvDt);
}

//------------------------------------------------------------------------------
// Vote on the time step
//------------------------------------------------------------------------------
template<typename Dimension>
typename SubPointPressureHourglassControl<Dimension>::TimeStepType
SubPointPressureHourglassControl<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  auto dtMin = std::numeric_limits<Scalar>::max();
  size_t nodeListMin = 0u;
  size_t iMin = 0u;
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto numNodeLists = mDvDt.size();
  CHECK(H.size() == numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = H[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto ahat = mDvDt(k,i).unitVector();
      const auto hi = 1.0/(H(k,i).dot(ahat)).magnitude();
      const auto dti = hi*safeInvVar(mDvDt(k,i).magnitude(), 1.0e-10);
      if (dti < dtMin) {
        dtMin = dti;
        nodeListMin = k;
        iMin = i;
      }
    }
  }
  dtMin = 0.5*sqrt(dtMin)*safeInvVar(mfHG);
  return std::make_pair(dtMin, std::string("SubPointPressureHourglassControl on point (" + to_string(nodeListMin) + " " + to_string(iMin) + ")"));
}

//------------------------------------------------------------------------------
// Add our terms to the hydro derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("SubPointHGevalDerivs");

  // The connectivity.
  const auto& cm = dataBase.connectivityMap();
  const auto& nodeLists = cm.nodeLists();
  const auto& pairs = cm.nodePairList();
  const auto  numNodeLists = nodeLists.size();
  const auto  npairs = pairs.size();
  // const auto  nint = dataBase.numInternalNodes();
  CONTRACT_VAR(npairs);

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto P = state.fields(HydroFieldNames::pressure, 0.0);
  const auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  const auto cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, vector<CellFaceFlag>());
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  const auto gradRho = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero());
  CHECK(mass.size() == numNodeLists);
  CHECK(pos.size() == numNodeLists);
  CHECK(rho.size() == numNodeLists);
  CHECK(P.size() == numNodeLists);
  CHECK(cells.size() == numNodeLists);
  CHECK(cellFaceFlags.size() == numNodeLists);
  CHECK(gradRho.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto& pairAccelerations = derivs.get(HydroFieldNames::pairAccelerations, vector<Vector>());
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  const auto compatibleEnergy = pairAccelerations.size() > 0u;
  CHECK((not compatibleEnergy) or pairAccelerations.size() == npairs);

  // // BLAGO
  // {
  //   int nodeListi, i, nodeListj, j, cellFace;
  //   for (nodeListi = 0; nodeListi < int(numNodeLists); ++nodeListi) {
  //     const int n = cells[nodeListi]->numInternalElements();
  //     for (i = 0; i < n; ++i) {
  //       // const bool barf = Process::getRank() == 0 and i == 0;
  //       const auto& celli = cells(nodeListi, i);
  //       const auto& xi = pos(nodeListi, i);
  //       const auto  mi = mass(nodeListi, i);
  //       const auto  Pi = P(nodeListi, i);
  //       const auto  rhoi = P(nodeListi, i);
  //       // const auto& gradRhoi = gradRho(nodeListi, i);
  //       const auto  comi = celli.centroid(); // centerOfMass(celli, xi, rhoi, gradRhoi);
  //       // if (barf) cerr << i << " " << cellFaceFlags(nodeListi, i).size() << endl;
  //       for (const auto& flags: cellFaceFlags(nodeListi,i)) {
  //         cellFace = flags.cellFace;
  //         nodeListj = flags.nodeListj;
  //         j = flags.j;
  //         CHECK(nodeListj != -1 or (nodeListj == -1 and j == -1));
  //         // if (barf) cerr << cellFace << " " << nodeListj << " " << j << " : ";
  //         const auto deltaDvDtij = mfHG * subCellAcceleration(celli, cellFace, comi, xi, Pi)/mi;
  //         DvDt(nodeListi, i) += deltaDvDtij;
  //         DepsDt(nodeListi, i) -= vel(nodeListi, i).dot(deltaDvDtij);
  //         // if (barf) cerr << "[" << i << " " << j << "] : " << deltaDvDtij << " " << deltaDvDtij.dot(comi - xi) << " : " << DvDt(nodeListi, i) << " " << DvDt(nodeListj, j);
  //         // if (barf) cerr << endl;
  //       }
  //     }
  //   }
  // }

  // Walk the cell face flags, looking for pair interactions
  mDvDt = Vector::zero();
  {
    int nodeListi, i, nodeListj, j, cellFace;
    for (nodeListi = 0; nodeListi < int(numNodeLists); ++nodeListi) {
      const int n = cellFaceFlags[nodeListi]->numInternalElements();
      for (i = 0; i < n; ++i) {
        if (surfacePoint(nodeListi, i) == 0) {
          // const bool barf = Process::getRank() == 0 and i == 100;
          const auto& celli = cells(nodeListi, i);
          const auto& xi = pos(nodeListi, i);
          const auto  Pi = P(nodeListi, i);
          const auto  mi = mass(nodeListi, i);
          const auto  rhoi = rho(nodeListi, i);
          // const auto& gradRhoi = gradRho(nodeListi, i);
          const auto  comi = celli.centroid(); // centerOfMass(celli, xi, rhoi, gradRhoi);
          // if (barf) cerr << i << " " << xi << " " << cellFaceFlags(nodeListi, i).size() << endl;
          for (const auto& flags: cellFaceFlags(nodeListi,i)) {
            cellFace = flags.cellFace;
            nodeListj = flags.nodeListj;
            j = flags.j;
            CHECK(nodeListj != -1 or (nodeListj == -1 and j == -1));
            // const bool barf = (Process::getRank() == 0 and (i == 8 or j == 8));
            // if (barf) cerr << cellFace << " " << nodeListj << " " << j << " : ";
            if (nodeListj != -1 and     // Avoid external faces (with void)
                cm.calculatePairInteraction(nodeListi, i, nodeListj, j, nodeLists[nodeListj]->firstGhostNode())) {  // make sure we hit each pair only once
              const auto& cellj = cells(nodeListj, j);
              const auto& xj = pos(nodeListj, j);
              const auto  Pj = P(nodeListj, j);
              const auto  mj = mass(nodeListj, j);
              const auto  rhoj = rho(nodeListj, j);
              const auto  comj = cellj.centroid();
              const auto  aij =  mfHG * (subCellAcceleration(celli, cellFace, comi, xi, Pi)/rhoi +
                                         subCellAcceleration(celli, cellFace, comj, xj, Pj)*mj/(mi*rhoj));
              const auto  aji = -aij*mi/mj;
              CHECK2(fuzzyGreaterThanOrEqual(subCellAcceleration(celli, cellFace, comi, xi, Pi).dot(comi - xi), 0.0, 1.0e-5), 
                     subCellAcceleration(celli, cellFace, comi, xi, Pi)/rhoi << " " << xi << " " << comi << " : " << cellFace << " " << celli << " " << Pi << " : " << subCellAcceleration(celli, cellFace, comi, xi, Pi).dot(comi - xi)/rhoi);
              // const auto aij = mfHG * (subCellAcceleration(celli, cellFace, comi, xi, Pi)/mi +
              //                          subCellAcceleration(celli, cellFace, comj, xj, Pj)/mj);
              // const auto aji = -aij * mi/mj;
              // const bool barf = (Process::getRank() == 0 and (i == 8 or j == 8));
              // if (barf) {
              //   cerr << " --> " << i << " " << j << " : " << xi << " " << xj << " : " << comi << " " << comj << " : "
              //        // << celli << " " << cellj << " : "
              //        << aij << " " << aji << " : "
              //        << subCellAcceleration(celli, cellFace, comi, xi, Pi) << " " << subCellAcceleration(celli, cellFace, comj, xj, Pj) << " : "
              //        << aij.dot(comi - xi) << " " << aji.dot(comj - xj) << "\n";
              // }
              DvDt(nodeListi, i) += aij;
              DvDt(nodeListj, j) += aji;
              mDvDt(nodeListi, i) += aij;
              mDvDt(nodeListj, j) += aji;
              DepsDt(nodeListi, i) -= vel(nodeListi, i).dot(aij);
              DepsDt(nodeListj, j) -= vel(nodeListj, j).dot(aji);
              if (compatibleEnergy) {
                const auto kk = pairs.index(NodePairIdxType(i, nodeListi, j, nodeListj));
                CHECK((nodeListi == int(pairs[kk].i_list) and i == int(pairs[kk].i_node)) or
                      (nodeListi == int(pairs[kk].j_list) and i == int(pairs[kk].j_node)));
                const bool flip = (nodeListi == int(pairs[kk].j_list) and i == int(pairs[kk].j_node));
                pairAccelerations[kk] += (flip ? aji : aij);
              }
              // if (barf) cerr << "[" << i << " " << j << "] : " << aij << " " << aij.dot(comi - xi) << " : " << DvDt(nodeListi, i) << " " << DvDt(nodeListj, j);
            }
            // if (barf) cerr << endl;
          }

          // Optionally add direct filtering to the position update
          const auto vhat = vel(nodeListi, i).unitVector();
          const auto vcent = (celli.centroid() - xi).dot(vhat)*safeInv(dt) * vhat;
          const auto fcent = min(1.0, mxfilter*vel(nodeListi, i).magnitude()*safeInv(vcent.magnitude()));
          DxDt(nodeListi, i) += fcent * vcent;
        }
      }
    }
  }

  // // Scan the accelerations to build our timestep vote
  // mDtMin = std::numeric_limits<double>::max();
  // const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  // CHECK(H.size() == numNodeLists);
  // for (auto k = 0u; k < numNodeLists; ++k) {
  //   const auto n = H[k]->numInternalElements();
  //   for (auto i = 0u; i < n; ++i) {
  //     const auto ahat = DvDt(k,i).unitVector();
  //     const auto hi = 1.0/(H(k,i).dot(ahat)).magnitude();
  //     const auto dti = hi*safeInvVar(DvDt(k,i).magnitude(), 1.0e-10);
  //     if (dti < mDtMin) {
  //       mDtMin = dti;
  //       mNodeListMin = k;
  //       mImin = i;
  //     }
  //   }
  // }
  // mDtMin = 0.01*sqrt(mDtMin);

  TIME_END("SubPointHGevalDerivs");
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mDvDt, pathName + "/DvDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mDvDt, pathName + "/DvDt");
}

} // end namespace Spheral
