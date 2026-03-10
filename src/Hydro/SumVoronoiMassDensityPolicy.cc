//---------------------------------Spheral++----------------------------------//
// SumVoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase 
// specialized for the updating the mass density according to the specific 
// volume from the Voronoi tesselation.  This version uses a weighted sum of the
// local mass and volume to do the deed.
//
// Created by JMO, Mon Aug  1 10:48:03 PDT 2011
//----------------------------------------------------------------------------//
#include "SumVoronoiMassDensityPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SumVoronoiMassDensityPolicy<Dimension>::
SumVoronoiMassDensityPolicy(const TableKernel<Dimension>& W, 
                            const Physics<Dimension>& package,
                            const double rhoMin, const double rhoMax):
  UpdatePolicyBase<Dimension>({HydroFieldNames::mass,
                               HydroFieldNames::volume}),
  mW(W),
  mPackage(package),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax) {
  REQUIRE(rhoMin <= rhoMax);
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SumVoronoiMassDensityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto massDensity = state.fields(fieldKey, Scalar());
  const auto numFields = massDensity.numFields();
  const auto W0 = mW.kernelValue(0.0, 1.0);

  // Intitialize the mass density.
  massDensity = 0.0;

  // Grab the state we need to do our job.
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto& cm = state.connectivityMap();
  const auto& pairs = cm.nodePairList();
  const auto  npairs = pairs.size();

  // // The volume still needs to have boundary conditions enforced.
  // typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  // for (ConstBoundaryIterator boundaryItr = mPackage.boundaryBegin();
  //      boundaryItr != mPackage.boundaryEnd();
  //      ++boundaryItr) {
  //   (*boundaryItr)->applyGhostBoundary(mass);
  //   (*boundaryItr)->applyGhostBoundary(volume);
  // }
  // for (ConstBoundaryIterator boundaryItr = mPackage.boundaryBegin();
  //      boundaryItr != mPackage.boundaryEnd();
  //      ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // Prepare the effective volume storage.
  FieldList<Dimension, Scalar> volEff(mass);
  volEff.copyFields();
  volEff = 0.0;

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Hdeti, Hdetj, etai, etaj, Wi, Wj;
    Vector xij;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto volEff_thread = volEff.threadCopy(threadStack);

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      const auto& Vi = volume(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& xi = pos(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      Hdeti = Hi.Determinant();
      auto& rhoi = massDensity(nodeListi, i);
      auto& Veffi = volEff_thread(nodeListi, i);

      const auto& Vj = volume(nodeListi, j);
      const auto& mj = mass(nodeListi, j);
      const auto& xj = pos(nodeListi, j);
      const auto& Hj = H(nodeListi, j);
      Hdetj = Hj.Determinant();
      auto& rhoj = massDensity(nodeListi, j);
      auto& Veffj = volEff_thread(nodeListj, j);

      xij = xi - xj;
      etai = (Hi*xij).magnitude();
      etaj = (Hj*xij).magnitude();
      Wi = mW.kernelValue(etai, Hdeti);
      Wj = mW.kernelValue(etaj, Hdetj);

      Veffi += Vj*Wi;
      rhoi +=  mj*Wi;

      Veffj += Vi*Wj;
      rhoj +=  mi*Wj;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }

  // Finalize the density for node i.
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto ni = mass[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto& Vi = volume(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto Hdeti = Hi.Determinant();
      auto& rhoi = massDensity(nodeListi, i);
      rhoi = max(mRhoMin, min(mRhoMax, (rhoi + mi*Hdeti*W0) * safeInv(volEff(nodeListi, i) + Vi*Hdeti*W0)));
    }
  }
}

//------------------------------------------------------------------------------
// Update the density as an increment.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SumVoronoiMassDensityPolicy<Dimension>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double /*t*/,
                  const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the matching derivative field from the StateDerivatives.
  auto incrementKey = IncrementState<Dimension, Scalar>::prefix() + fieldKey;
  auto       f = state.fields(fieldKey, 0.0);
  const auto df = derivs.fields(incrementKey, 0.0);

  const auto numFields = f.numFields();
  for (unsigned i = 0u; i < numFields; ++i) {
    const auto ni = f[i]->numInternalElements();
#pragma omp parallel for
    for (auto j = 0u; j < ni; ++j) {
      f(i,j) = max(mRhoMin, min(mRhoMax, f(i,j) + multiplier*(df(i,j))));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SumVoronoiMassDensityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const SumVoronoiMassDensityPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

