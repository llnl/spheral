//---------------------------------Spheral++----------------------------------//
// DamageModel -- Base class for the damage physics models.
// This class just provides the basic interface for damage models, and does 
// not fill out the complete physics package interface.
//
// Created by JMO, Thu Sep 29 13:31:57 PDT 2005
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DamageModel.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/NodeCoupling.hh"
#include "Damage/PairMaxDamageNodeCoupling.hh"
#include "Damage/DamageGradientNodeCoupling.hh"
#include "Damage/ThreePointDamagedNodeCoupling.hh"
#include "Utilities/Timer.hh"

#include <string>
#include <vector>
#include <sstream>

using std::vector;
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

//------------------------------------------------------------------------------
// Provide to_string for tensors
//------------------------------------------------------------------------------
template<int nDim>
std::string to_string(const GeomSymmetricTensor<nDim>& x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamageModel<Dimension>::
DamageModel(SolidNodeList<Dimension>& nodeList,
            const TableKernel<Dimension>& W,
            const double crackGrowthMultiplier,
            const DamageCouplingAlgorithm damageCouplingAlgorithm):
  Physics<Dimension>(),
  mNodeList(nodeList),
  mW(W),
  mCrackGrowthMultiplier(crackGrowthMultiplier),
  mDamageCouplingAlgorithm(damageCouplingAlgorithm),
  mExcludeNode("Nodes excluded from damage", nodeList, 0),
  mNodeCouplingPtr(new NodeCoupling()),
  mComputeIntersectConnectivity(false),
  mFreezeDamage(false),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Compute the generic Grady-Kipp (ala Benz-Asphaug) scalar damage time 
// derivative.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
computeScalarDDDt(const DataBase<Dimension>& /*dataBase*/,
                  const State<Dimension>& state,
                  const Scalar /*time*/,
                  const Scalar /*dt*/,
                  Field<Dimension, Scalar>& DDDt) const {

  // Pre-conditions.
  REQUIRE(DDDt.nodeListPtr() == &mNodeList);

  if (mFreezeDamage) {

    DDDt = 0.0;

  } else {

    // Get the state fields.
    const auto  clKey = State<Dimension>::buildFieldKey(SolidFieldNames::longitudinalSoundSpeed, mNodeList.name());
    const auto  HKey = State<Dimension>::buildFieldKey(HydroFieldNames::H, mNodeList.name());
    const auto& cl = state.field(clKey, 0.0);
    const auto& H = state.field(HKey, SymTensor::zero());

    // Constant multiplicative parameter for the crack growth.
    const auto A = mCrackGrowthMultiplier / mW.kernelExtent();

    // Iterate over the internal nodes.
    const auto ni = mNodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      if (mExcludeNode(i) == 1) {

        DDDt(i) = 0.0;

      } else {

        const double hrInverse = Dimension::rootnu(H(i).Determinant());
        DDDt(i) = A * cl(i) * hrInverse;

      }
      CHECK(DDDt(i) >= 0.0);
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  for (int i = 0; i != (int)mNodeList.numInternalNodes(); ++i) {
    if (mExcludeNode(i) == 1) {
      ENSURE(DDDt(i) == 0.0);
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// initialize (before evaluateDerivatives)
// This is where we update the damage coupling, which is stored in the
// connectivity pair list in f_couple.
// We can't do this during registerState because some coupling algorithms
// require ghost state to be updated first.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DamageModel<Dimension>::
initialize(const Scalar /*time*/,
           const Scalar /*dt*/,
           const DataBase<Dimension>& /*dataBase*/,
           State<Dimension>& state,
           StateDerivatives<Dimension>& /*derivs*/) {

  auto& connectivity = state.connectivityMap();
  auto& pairs = const_cast<NodePairList&>(connectivity.nodePairList());

  switch(mDamageCouplingAlgorithm) {
  case DamageCouplingAlgorithm::DirectDamage:
    break;

  case DamageCouplingAlgorithm::PairMaxDamage:
    mNodeCouplingPtr = std::make_shared<PairMaxDamageNodeCoupling<Dimension>>(state, pairs);
    break;

  case DamageCouplingAlgorithm::DamageGradient:
    mNodeCouplingPtr = std::make_shared<DamageGradientNodeCoupling<Dimension>>(state, mW, this->boundaryBegin(), this->boundaryEnd(), pairs);
    break;

  case DamageCouplingAlgorithm::ThreePointDamage:
    mNodeCouplingPtr = std::make_shared<ThreePointDamagedNodeCoupling<Dimension>>(state, mW, pairs);
    break;

  default:
    VERIFY2(false, "DamageModel ERROR: unhandled damage coupling algorithm case");
  }
  connectivity.coupling(mNodeCouplingPtr);
  return false;
}

//------------------------------------------------------------------------------
// finalize (end of physics step)
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
finalize(const Scalar /*time*/, 
         const Scalar /*dt*/,
         DataBase<Dimension>& dataBase, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& /*derivs*/) {

  TIME_BEGIN("DamageModel_finalize");

  // For 3pt damage, check if we should switch to using full intersection data
  // from the ConnectivityMap
  if (mDamageCouplingAlgorithm == DamageCouplingAlgorithm::ThreePointDamage) {
    const auto D = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());
    auto nD = 0u;
    const auto numNodeLists = D.numFields();
#pragma omp parallel
    {
      auto nD_thread = 0u;
      for (auto il = 0u; il < numNodeLists; ++il) {
        const auto n = D[il]->numInternalElements();
#pragma omp for
        for (auto i = 0u; i < n; ++i) {
          if (D(il,i).Trace() > 1.0e-3) ++nD_thread;
        }
      }
#pragma omp critical
      {
        nD += nD_thread;
      }
    }
    nD = allReduce(nD, SPHERAL_OP_SUM);
    const auto ntot = std::max<size_t>(1u, dataBase.globalNumInternalNodes());
    const auto dfrac = double(nD)/double(ntot);
    mComputeIntersectConnectivity = (dfrac > 0.2);  // Should tune this number...
    // if (Process::getRank() == 0) std::cout << "DamageModel dfrac = " << nD << "/" << ntot << " = " << dfrac << " : " << mComputeIntersectConnectivity << std::endl;
  }
  TIME_END("DamageModel_finalize");
}

//------------------------------------------------------------------------------
// Return the maximum state change we care about for checking for convergence
// in the implicit integration methods.
//------------------------------------------------------------------------------
template<typename Dimension>
typename DamageModel<Dimension>::ResidualType
DamageModel<Dimension>::
maxResidual(const DataBase<Dimension>& dataBase, 
            const State<Dimension>& state1,
            const State<Dimension>& state0,
            const Scalar tol) const {
  REQUIRE(tol > 0.0);

  // Define some functions to compute residuals
  auto fresT = [](const SymTensor& x1, const SymTensor& x2, const Scalar tol) { auto dx = std::abs((x2 - x1).Trace()); return dx/std::max(std::abs(x1.Trace()) + std::abs(x2.Trace()), Dimension::nDim*tol); };

  // Initialize the return value to some impossibly high value.
  auto result = ResidualType(-1.0, "You should not see me!");

  // Grab the state we're comparing
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, mNodeList.name()); };
  const auto& D0 = state0.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());
  const auto& D1 = state1.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());
  
  // Walk the nodes
  const auto n = mNodeList.numInternalNodes();
  const auto rank = Process::getRank();
#pragma omp parallel
  {
    auto maxResidual_local = result;
#pragma omp for
    for (auto i = 0u; i < n; ++i) {

      // We limit by the change in damage
      const auto Dres = fresT(D0(i), D1(i), tol);
      if (Dres > maxResidual_local.first) {
        maxResidual_local = ResidualType(Dres, ("Damage change: residual = " + std::to_string(Dres) + "\n" +
                                                "                     D0 = " + to_string(D0(i)) + 
                                                "                     D1 = " + to_string(D1(i)) + 
                                                "      (nodeList, i, rank) = (" + mNodeList.name() + " " + std::to_string(i) + " " + std::to_string(rank) + ")\n"));
      }
    }

#pragma omp critical
    {
      if (maxResidual_local.first > result.first) {
        result = maxResidual_local;
      }
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Access the set of nodes to be excluded from damage.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
DamageModel<Dimension>::
excludeNodes() const {
  vector<int> result;
  for (auto i = 0u; i != mNodeList.numInternalNodes(); ++i) {
    if (mExcludeNode(i) == 1.0) result.push_back(i);
  }
  return result;
}

template<typename Dimension>
void
DamageModel<Dimension>::
excludeNodes(vector<int> ids) {
  mExcludeNode = 0;
  for (vector<int>::const_iterator itr = ids.begin();
       itr != ids.end();
       ++itr) {
    REQUIRE(*itr >= 0 && *itr < (int)mNodeList.numInternalNodes());
    mExcludeNode(*itr) = 1;
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mCrackGrowthMultiplier, pathName + "/crackGrowthMultiplier");
  file.write(mExcludeNode, pathName + "/excludeNode");
  file.write(mComputeIntersectConnectivity, pathName + "/computeIntersectConnectivity");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mCrackGrowthMultiplier, pathName + "/crackGrowthMultiplier");
  file.read(mExcludeNode, pathName + "/excludeNode");
  file.read(mComputeIntersectConnectivity, pathName + "/computeIntersectConnectivity");
}

}

