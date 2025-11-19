//---------------------------------Spheral++----------------------------------//
// CRKSPHRZ -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Thu May 12 15:25:24 PDT 2016
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/RKFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "RK/ContinuityVolumePolicyRZ.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/range.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/GeometryRegistrar.hh"

#include "CRKSPHRZ.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
CRKSPHRZ::
CRKSPHRZ(DataBase<Dimension>& dataBase,
                  ArtificialViscosityHandle<Dimension>& Q,
                  const RKOrder order,
                  const double cfl,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool XSPH,
                  const MassDensityType densityUpdate,
                  const double epsTensile,
                  const double nTensile):
  CRKSPHBase<Dim<2>>(dataBase,
                     Q,
                     order,
                     cfl,
                     useVelocityMagnitudeForDt,
                     compatibleEnergyEvolution,
                     evolveTotalEnergy,
                     XSPH,
                     densityUpdate,
                     epsTensile,
                     nTensile),
  mPairAccelerationsPtr() {
}

//------------------------------------------------------------------------------
// On problem start up some dependent state needs to be calculated
//------------------------------------------------------------------------------
void
CRKSPHRZ::
initializeProblemStartupDependencies(DataBase<Dim<2>>& dataBase,
                                     State<Dim<2>>& state,
                                     StateDerivatives<Dim<2>>& derivs) {
  
  // Correct the mass to mass/r.
  auto mass = dataBase.fluidMass();
  const auto pos = dataBase.fluidPosition();
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Do general initializations.
  CRKSPHBase<Dim<2>>::initializeProblemStartupDependencies(dataBase, state, derivs);

  // Convert back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
CRKSPHRZ::
registerState(DataBase<Dim<2>>& dataBase,
              State<Dim<2>>& state) {

  // The base class does most of it.
  CRKSPHBase<Dim<2>>::registerState(dataBase, state);

  // Reregister the volume update
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  state.enroll(vol, make_policy<ContinuityVolumePolicyRZ>());

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy.
  // Note in RZ we require the specific thermal energy go before the position so we can use the r position
  // during update.  This is why we make position update dependent on the thermal energy in SPHBase.
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    state.enroll(specificThermalEnergy, make_policy<RZNonSymmetricSpecificThermalEnergyPolicy>(dataBase));

  } else if (evolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
void
CRKSPHRZ::
registerDerivatives(DataBase<Dim<2>>& dataBase,
                    StateDerivatives<Dim<2>>& derivs) {
  CRKSPHBase<Dimension>::registerDerivatives(dataBase, derivs);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
  }
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
void
CRKSPHRZ::
preStepInitialize(const DataBase<Dim<2>>& dataBase, 
                  State<Dim<2>>& state,
                  StateDerivatives<Dim<2>>& derivs) {

  // Convert the mass to mass per unit length first.
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Base class does most of the work.
  CRKSPHBase<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Now convert back to true masses.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto& xi = pos(nodeListi, i);
      const auto circi = 2.0*M_PI*abs(xi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
CRKSPHRZ::
evaluateDerivatives(const Dimension::Scalar time,
                    const Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Depending on the type of the ArtificialViscosity, dispatch the call to
  // the secondDerivativesLoop
  auto& Qhandle = this->artificialViscosity();
  if (Qhandle.QPiTypeIndex() == std::type_index(typeid(Scalar))) {
      const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Scalar>&>(Qhandle);
      this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  } else {
    CHECK(Qhandle.QPiTypeIndex() == std::type_index(typeid(Tensor)));
    const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Tensor>&>(Qhandle);
    this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  }
}
  
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename QType>
void
CRKSPHRZ::
evaluateDerivativesImpl(const Dim<2>::Scalar /*time*/,
                        const Dim<2>::Scalar /*dt*/,
                        const DataBase<Dim<2>>& dataBase,
                        const State<Dim<2>>& state,
                        StateDerivatives<Dim<2>>& derivs,
                        const QType& Q) const {

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  //const auto  order = this->correctionOrder();
  const auto& WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(mOrder));

  // A few useful constants we'll use in the following loop.
  //const auto tiny = 1.0e-30;
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto XSPH = this->XSPH();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto corrections = state.fields(RKFieldNames::rkCorrections(mOrder), RKCoefficients<Dimension>());
  const auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  const auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero(), true);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(corrections.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK((compatibleEnergy and pairAccelerationsPtr->size() == npairs) or not compatibleEnergy);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, Wj, Qi, Qj;
    QPiType QPiij, QPiji;
    Vector gradWi, gradWj;
    Vector deltagrad, forceij, forceji;
    Vector xij, vij, etai, etaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto  ri = abs(posi.y());
      const auto  circi = 2.0*M_PI*ri;
      const auto  mi = mass(nodeListi, i);
      const auto  mRZi = mi/circi;
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      const auto  weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      const auto  zetai = abs((Hi*posi).y());
      //const auto  hri = ri*safeInv(zetai);
      CHECK2(ri > 0.0, i << " " << ri);
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);
      CHECK2(weighti > 0.0, i << " " << weighti);

      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);

      // Get the state for node j
      const auto& posj = position(nodeListj, j);
      const auto  rj = abs(posj.y());
      const auto  circj = 2.0*M_PI*rj;
      const auto  mj = mass(nodeListj, j);
      const auto  mRZj = mj/circj;
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto& correctionsj = corrections(nodeListj, j);
      const auto  weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
      const auto  zetaj = abs((Hj*posj).y());
      CHECK2(rj > 0.0, j << " " << rj);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(weightj > 0.0);

      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);

      // Node displacement.
      xij = posi - posj;
      vij = vi - vj;
      etai = Hi*xij;
      etaj = Hj*xij;

      // Symmetrized kernel weight and gradient.
      std::tie(Wj, gradWj) = WR.evaluateKernelAndGradient( xij, Hj, correctionsi);  // Hj because we compute RK using scatter formalism
      std::tie(Wi, gradWi) = WR.evaluateKernelAndGradient(-xij, Hi, correctionsj);
      deltagrad = gradWj - gradWi;

      // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
      Q.QPiij(QPiij, QPiji, Qi, Qj,
              nodeListi, i, nodeListj, j,
              posi, Hi, etai, vi, rhoi, ci,  
              posj, Hj, etaj, vj, rhoj, cj,
              fClQ, fCqQ, DvDxQ); 
      const auto Qaccij = (rhoi*rhoi*QPiij + rhoj*rhoj*QPiji)*deltagrad;
      const auto workQi = (rhoj*rhoj*QPiji*vij).dot(deltagrad);                // CRK
      const auto workQj = (rhoi*rhoi*QPiij*vij).dot(deltagrad);                // CRK
      maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                     // We need tighter timestep controls on the Q with CRK
      maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
      effViscousPressurei += weightj * Qi * Wj;
      effViscousPressurej += weighti * Qj * Wi;

      // Velocity gradient.
      DvDxi -= weightj*vij.dyad(gradWj);
      DvDxj += weighti*vij.dyad(gradWi);
      if (nodeListi == nodeListj) {
        localDvDxi -= weightj*vij.dyad(gradWj);
        localDvDxj += weighti*vij.dyad(gradWi);
      }

      // Acceleration (CRKSPH form).
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto forceij  = 0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij); // <- Type III, with CRKSPH Q forces
      DvDti -= forceij/mRZi; //CRK Acceleration
      DvDtj += forceij/mRZj; //CRK Acceleration
      if (compatibleEnergy) {
        (*pairAccelerationsPtr)[kk][0] = -forceij/mRZi;
        (*pairAccelerationsPtr)[kk][1] =  forceij/mRZj;
      }

      DepsDti += 0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + workQi)/mRZi;    // CRK Q
      DepsDtj += 0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + workQj)/mRZj;    // CRK Q

      // Estimate of delta v (for XSPH).
      if ((XSPH and (nodeListi == nodeListj)) or min(zetai, zetaj) < 1.0) {
        XSPHDeltaVi -= weightj*Wj*vij;
        XSPHDeltaVj += weighti*Wi*vij;
      }
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dim<2>>(threadStack);

  }   // OMP parallel

  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto  ri = abs(posi.y());
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  zetai = abs((Hi*posi).y());
      const auto  hri = ri*safeInv(zetai);
      const auto  riInv = safeInv(ri, 0.25*hri);
      CHECK2(ri > 0.0, i << " " << ri);
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);

      // Time evolution of the mass density.
      const auto vri = vi.y(); // + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(DvDxi.Trace() + vri*riInv);

      // Finish the specific thermal energy evolution.
      DepsDti -= Pi/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
CRKSPHRZ::
applyGhostBoundaries(State<Dim<2>>& state,
                     StateDerivatives<Dim<2>>& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero());
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary CRKSPH BCs.
  CRKSPHBase<Dim<2>>::applyGhostBoundaries(state, derivs);
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Scale back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
CRKSPHRZ::
enforceBoundaries(State<Dim<2>>& state,
                  StateDerivatives<Dim<2>>& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero());
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary CRKSPH BCs.
  CRKSPHBase<Dim<2>>::enforceBoundaries(state, derivs);

  // Scale back to mass.
  // We also ensure no point approaches the z-axis too closely.
  FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    //const Scalar nPerh = mass[nodeListi]->nodeList().nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      Vector& posi = pos(nodeListi, i);
      const Scalar circi = 2.0*M_PI*abs(posi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

}
