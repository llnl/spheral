//---------------------------------Spheral++----------------------------------//
// SPHHydroBaseGSRZ -- The SPH/ASPH hydrodynamic package for Spheral++,
//                   specialized for 2D RZ (cylindrical) geometry.
//
// Based on the methodology described in the dissertation
// de Catalunya Departament de Física i Enginyeria Nuclear, U. P., García Senz, D., (2012).
// AxisSPH:devising and validating an axisymmetric smoothed particle hydrodynamics code.
// TDX (Tesis Doctorals en Xarxa). Universitat Politècnica de Catalunya.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Tue Apr 26 16:06:21 PDT 2016
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosityView.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

#include "SPHHydroBaseGSRZ.hh"

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
SPHHydroBaseGSRZ::
SPHHydroBaseGSRZ(const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
                 DataBase<Dimension>& dataBase,
                 ArtificialViscosityView<Dim<2> >& Q,
                 const TableKernel<Dim<2> >& W,
                 const TableKernel<Dim<2> >& WPi,
                 const double filter,
                 const double cfl,
                 const bool useVelocityMagnitudeForDt,
                 const bool compatibleEnergyEvolution,
                 const bool evolveTotalEnergy,
                 const bool gradhCorrection,
                 const bool XSPH,
                 const bool correctVelocityGradient,
                 const bool sumMassDensityOverAllNodeLists,
                 const MassDensityType densityUpdate,
                 const HEvolutionType HUpdate,
                 const double epsTensile,
                 const double nTensile,
                 const Vector& xmin,
                 const Vector& xmax):
  SPHHydroBase<Dim<2> >(smoothingScaleMethod,
                        dataBase,
                        Q,
                        W,
                        WPi,
                        filter,
                        cfl,
                        useVelocityMagnitudeForDt,
                        compatibleEnergyEvolution,
                        evolveTotalEnergy,
                        gradhCorrection,
                        XSPH,
                        correctVelocityGradient,
                        sumMassDensityOverAllNodeLists,
                        densityUpdate,
                        HUpdate,
                        epsTensile,
                        nTensile,
                        xmin,
                        xmax) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SPHHydroBaseGSRZ::
~SPHHydroBaseGSRZ() {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SPHHydroBaseGSRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // The base class does most of it.
  SPHHydroBase<Dim<2> >::registerState(dataBase, state);

  // Are we using the compatible energy evolution scheme?
  // If so we need to override the ordinary energy registration with a specialized version.
  if (mCompatibleEnergyEvolution) {
    FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    PolicyPointer thermalEnergyPolicy(new RZNonSymmetricSpecificThermalEnergyPolicy(dataBase));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SPHHydroBaseGSRZ::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // If we're going to do the SPH summation density, we need to convert the mass
  // to mass per unit length first.
  if (densityUpdate() == MassDensityType::RigorousSumDensity or
      densityUpdate() == MassDensityType::CorrectedSumDensity) {
    FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
    const unsigned numNodeLists = mass.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = mass[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
        mass(nodeListi, i) *= safeInvVar(circi);
      }
    }
  }

  // Base class finalization does most of the work.
  SPHHydroBase<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Now convert back to true masses and mass densities.  We also apply the RZ
  // correction factor to the mass density.
  if (densityUpdate() == MassDensityType::RigorousSumDensity or
      densityUpdate() == MassDensityType::CorrectedSumDensity) {
    //const TableKernel<Dimension>& W = this->kernel();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const unsigned numNodeLists = massDensity.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = massDensity[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& xi = position(nodeListi, i);
        //const SymTensor& Hi = H(nodeListi, i);
        //const Scalar zetai = abs((Hi*xi).y());
        //const Scalar fi = W.f1(zetai);
        const Scalar circi = 2.0*M_PI*abs(xi.y());
        mass(nodeListi, i) *= circi;
        // massDensity(nodeListi, i) *= fi;
        // massDensity(nodeListi, i) *= fi*safeInvVar(2.0*M_PI*abs(xi.y()));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SPHHydroBaseGSRZ::
evaluateDerivatives(const Dim<2>::Scalar /*time*/,
                    const Dim<2>::Scalar /*dt*/,
                    const DataBase<Dim<2> >& dataBase,
                    const State<Dim<2> >& state,
                    StateDerivatives<Dim<2> >& derivatives) const {

  // Get the ArtificialViscosityView.
  ArtificialViscosityView<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const double rhoTiny = 1.0e-10;
  const double W0 = W(0.0, 1.0);

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  FieldList<Dimension, Tensor> localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, Scalar> effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  FieldList<Dimension, Scalar> viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Scalar> XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (auto i = 0u; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodeList = **itr;
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();

    // Iterate over the internal nodes in this NodeList.
    for (ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // RZ correction factors for node i.
      const Vector& posi = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar ri = abs(posi.y());
      const Scalar zetai = abs((Hi*posi).y());
      const Scalar hrInvi = zetai*safeInvVar(ri);
      Scalar f1i, f2i, gradf1i, gradf2i;
      W.f1Andf2(zetai, f1i, f2i, gradf1i, gradf2i);
      // f1i = 1.0; f2i = 1.0;
      // gradf1i = 0.0; gradf2i = 0.0;
      gradf1i *= hrInvi;
      gradf2i *= hrInvi;
      const Scalar circi = 2.0*M_PI*ri;
      const Scalar circInvi = safeInvVar(circi);
      const Scalar rhoi = f1i*massDensity(nodeListi, i);
      const Scalar rhoRZi = rhoi*circi;

      // Get the state for node i.
      const Scalar mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar vri = vi.y();
      const Scalar vzi = vi.x();
      const Scalar Pi = f1i*pressure(nodeListi, i);
      const Scalar ci = f1i*soundSpeed(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      Scalar& rhoSumi = rhoSum(nodeListi, i);
      Scalar& normi = normalization(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      Tensor& Mi = M(nodeListi, i);
      Tensor& localMi = localM(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      Scalar& effViscousPressurei = effViscousPressure(nodeListi, i);
      Scalar& viscousWorki = viscousWork(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Scalar& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#if defined __INTEL_COMPILER
#pragma vector always
#endif
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // RZ correction factors for node j.
              const Vector& posj = position(nodeListj, j);
              const Scalar rj = abs(posj.y());
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar zetaj = abs((Hj*posj).y());
              const Scalar hrInvj = zetaj*safeInvVar(rj);
              Scalar f1j, f2j, gradf1j, gradf2j;
              W.f1Andf2(zetaj, f1j, f2j, gradf1j, gradf2j);
              // f1j = 1.0; f2j = 1.0;
              // gradf1j = 0.0; gradf2j = 0.0;
              gradf1j *= hrInvj;
              gradf2j *= hrInvj;
              const Scalar circj = 2.0*M_PI*rj;
              const Scalar rhoj = f1j*massDensity(nodeListj, j);
              const Scalar rhoRZj = rhoj*circj;

              // Get the state for node j
              const Scalar mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar vrj = vj.y();
              const Scalar vzj = vj.x();
              const Scalar Pj = f1j*pressure(nodeListj, j);
              const Scalar cj = f1j*soundSpeed(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Scalar& normj = normalization(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Tensor& Mj = M(nodeListj, j);
              Tensor& localMj = localM(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              Scalar& effViscousPressurej = effViscousPressure(nodeListj, j);
              Scalar& viscousWorkj = viscousWork(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Scalar& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
              Vector& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Node displacement.
              const Vector xij = posi - posj;
              const Vector etai = Hi*xij;
              const Vector etaj = Hj*xij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // Symmetrized kernel weight and gradient.
              const Vector Hetai = Hi*etai.unitVector();
              Scalar Wi, gWi;
              W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
              const Vector gradWi = gWi*Hetai;
              const Vector gradWQi = WQ.gradValue(etaMagi, Hdeti) * Hetai;

              const Vector Hetaj = Hj*etaj.unitVector();
              Scalar Wj, gWj;
              W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
              const Vector gradWj = gWj*Hetaj;
              const Vector gradWQj = WQ.gradValue(etaMagj, Hdetj) * Hetaj;

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double xij2 = xij.magnitude2();
              const SymTensor thpt = xij.selfdyad()/max(tiny, xij2*FastMath::square(Dimension::pownu12(xij2)));
              weightedNeighborSumi += fweightij*std::abs(gWi);
              weightedNeighborSumj += fweightij*std::abs(gWj);
              massSecondMomenti += fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density.
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
                normi += mi/rhoRZi*Wi;
                normj += mj/rhoRZj*Wj;
              }

              // Compute the pair-wise artificial viscosity.
              const Vector vij = vi - vj;
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        posi, etai, vi, rhoRZi, ci, Hi,
                                                        posj, etaj, vj, rhoRZj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWQi) * 2.0*M_PI*ri;
              const Vector Qaccj = 0.5*(QPiij.second*gradWQj) * 2.0*M_PI*rj;
              // const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWQi) * 2.0*M_PI*f1i*ri;
              // const Scalar workQj = 0.5*(QPiij.second*vij).dot(gradWQj) * 2.0*M_PI*f1j*rj;
              const Scalar workQi = vij.dot(Qacci);
              const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoRZi*rhoRZi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoRZj*rhoRZj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += mj/rhoRZj * Qi * Wi;
              effViscousPressurej += mi/rhoRZi * Qj * Wj;
              viscousWorki += mj*workQi;
              viscousWorkj += mi*workQj;

              // Acceleration.
              CHECK(rhoRZi > 0.0);
              CHECK(rhoRZj > 0.0);
              const double Prhoi = f1i*Pi*ri*safeInv(rhoRZi*rhoRZi, rhoTiny);
              const double Prhoj = f1j*Pj*rj*safeInv(rhoRZj*rhoRZj, rhoTiny);
              const Vector deltaDvDt = 2.0*M_PI*(Prhoi*gradWi + Prhoj*gradWj) + Qacci + Qaccj;
              DvDti -= mj*deltaDvDt;
              DvDtj += mi*deltaDvDt;

              // Specific thermal energy evolution.
              DepsDti += mj*(2.0*M_PI*Pi*ri*safeInv(rhoRZi*rhoRZi, rhoTiny)*
                             ( (f1i*vri - f2i*vrj)*gradWi.y() + f1i*(vzi - vzj)*gradWi.x() + (gradf1i*vri - gradf2i*vrj)*Wi) + workQi);
              DepsDtj += mi*(2.0*M_PI*Pj*rj*safeInv(rhoRZj*rhoRZj, rhoTiny)*
                             (-(f1j*vrj - f2j*vri)*gradWj.y() - f1j*(vzj - vzi)*gradWi.x() + (gradf1j*vrj - gradf2j*vri)*Wj) + workQj);
              if (mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(-mj*deltaDvDt);
                pairAccelerationsj.push_back( mi*deltaDvDt);
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = mj*vij.dyad(gradWi);
              const Tensor deltaDvDxj = mi*vij.dyad(gradWj);
              DvDxi -= deltaDvDxi;
              DvDxj -= deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= deltaDvDxi;
                localDvDxj -= deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (nodeListi == nodeListj) {
                const double fXSPH = max(0.0, min(1.0, abs(vij.dot(xij)*safeInv(vij.magnitude()*xij.magnitude()))));
                CHECK(fXSPH >= 0.0 and fXSPH <= 1.0);
                XSPHWeightSumi += fXSPH*mj*safeInv(rhoRZj, rhoTiny)*Wi;
                XSPHWeightSumj += fXSPH*mi*safeInv(rhoRZi, rhoTiny)*Wj;
                XSPHDeltaVi -= fXSPH*mj*safeInv(rhoRZj, rhoTiny)*Wi*vij;
                XSPHDeltaVj += fXSPH*mi*safeInv(rhoRZi, rhoTiny)*Wj*vij;
              }

              // Linear gradient correction term.
              Mi -= mj*xij.dyad(gradWi);
              Mj -= mi*xij.dyad(gradWj);
              if (nodeListi == nodeListj) {
                localMi -= mj*xij.dyad(gradWi);
                localMj -= mi*xij.dyad(gradWj);
              }
            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CONTRACT_VAR(firstGhostNodei);
      CHECK(not mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Add the self-contribution to density sum.
      rhoSumi = (rhoSumi + mi*W0*Hdeti)*f1i*circInvi;
      normi += mi*safeInv(rhoRZi, rhoTiny)*W0*Hdeti;

      // Finish the acceleration, adding the hoop terms.
      // DvDti.y(DvDti.y() + 1.01525*Pi/rhoRZi - Pi*ri/(rhoRZi*f1i)*gradf1i);   // Fiddled with a magic number on the hoop stress term.
      const Scalar deltaDvDtSelf = 2.0*M_PI*(Pi*safeInv(rhoRZi, rhoTiny) - Pi*ri*safeInv(rhoRZi*f1i, rhoTiny)*gradf1i);
      DvDti.y(DvDti.y() + deltaDvDtSelf);
      // DvDti.y(f1i*f1i*DvDti.y() - 0.5*(1.0 - f1i*f1i)*vri*safeInv(dt));
      if (mCompatibleEnergyEvolution) pairAccelerationsi.push_back(Vector(0.0, deltaDvDtSelf));

      // Finish the specific thermal energy derivative.
      DepsDti -= 2.0*M_PI*Pi*vri*safeInv(rhoRZi, rhoTiny);

      // If needed, convert to the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Finish the gradient of the velocity.
      CHECK(rhoRZi > 0.0);
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoRZi;
      }
      if (this->mCorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoRZi;
      }

      // Evaluate the continuity equation.
      DrhoDti = -rhoi*(DvDxi.Trace() + vri/ri);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      XSPHWeightSumi += Hdeti*mi/rhoRZi*W0;
      CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      if (mXSPH) {
        DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      } else {
        DxDti = f1i*f1i*vi + (1.0 - f1i*f1i)*(vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi));
      }
      // if (mXSPH) {
      //   XSPHWeightSumi += Hdeti*mi/rhoRZi*W0;
      //   CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      //   DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      // } else {
      //   DxDti = vi;
      // }

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             posi,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        posi,
                                                        weightedNeighborSumi,
                                                        massSecondMomenti,
                                                        W,
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        connectivityMap,
                                                        nodeListi,
                                                        i);
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SPHHydroBaseGSRZ::
applyGhostBoundaries(State<Dim<2> >& state,
                     StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary SPH BCs.
  SPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

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
SPHHydroBaseGSRZ::
enforceBoundaries(State<Dim<2> >& state,
                  StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary SPH BCs.
  SPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

  // Scale back to mass.
  // We also ensure no point approaches the z-axis too closely.
  FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    const Scalar nPerh = mass[nodeListi]->nodeList().nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      Vector& posi = pos(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar zetai = (Hi*posi).y();
      const Scalar ri = posi.y();
      const Scalar hrInvi = zetai*safeInvVar(ri);
      const Scalar rmin = 0.5/(nPerh*hrInvi);
      if (ri < rmin) posi.y(2.0*rmin - ri);
      const Scalar circi = 2.0*M_PI*abs(posi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

}
