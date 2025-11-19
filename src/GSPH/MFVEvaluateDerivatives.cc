namespace Spheral {


//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
secondDerivativesLoop(const typename Dimension::Scalar time,
                      const typename Dimension::Scalar dt,
                      const DataBase<Dimension>& dataBase,
                      const State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs) const {

  const auto& riemannSolver = this->riemannSolver();
  
  // A few useful constants we'll use in the following loop.
  const auto tiny = std::numeric_limits<Scalar>::epsilon();
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto gradType = this->gradientType();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  const auto  numNodeLists = nodeLists.size();
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();

  // kernel
  const auto& W = this->kernel();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto riemannDpDx = state.fields(GSPHFieldNames::RiemannPressureGradient,Vector::zero);
  const auto riemannDvDx = state.fields(GSPHFieldNames::RiemannVelocityGradient,Tensor::zero);
  
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(riemannDpDx.size() == numNodeLists);
  CHECK(riemannDvDx.size() == numNodeLists);

  // Derivative FieldLists.
  const auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  DrhoDx = derivs.fields(GSPHFieldNames::densityGradient, Vector::zero);
  const auto  DepsDx = derivs.fields(GSPHFieldNames::specificThermalEnergyGradient, Vector::zero);
  auto  normalization = derivs.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DvolDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::volume, 0.0);
  auto  DmDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::mass, 0.0);
  auto  DEDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy, 0.0);
  auto  DpDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + GSPHFieldNames::momentum, Vector::zero);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  newRiemannDpDx = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannPressureGradient,Vector::zero);
  auto  newRiemannDvDx = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannVelocityGradient,Tensor::zero);
  auto* pairAccelerationsPtr = (compatibleEnergy ?
                                &derivs.template get<PairAccelerationsType>(HydroFieldNames::pairAccelerations) :
                                nullptr);
  auto* pairDepsDtPtr = (compatibleEnergy ?
                         &derivs.template get<PairWorkType>(HydroFieldNames::pairWork) :
                         nullptr);
  auto* pairMassFluxPtr = (compatibleEnergy ?
                           &derivs.template get<PairMassFluxType>(GSPHFieldNames::pairMassFlux) :
                           nullptr);
  
  CHECK(DepsDx.size() == numNodeLists);
  CHECK(DrhoDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvolDt.size() == numNodeLists);
  CHECK(DEDt.size() == numNodeLists);
  CHECK(DpDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(newRiemannDpDx.size() == numNodeLists);
  CHECK(newRiemannDvDx.size() == numNodeLists);
  CHECK(not compatibleEnergy or pairAccelerationsPtr->size() == npairs);
  CHECK(not compatibleEnergy or pairDepsDtPtr->size() == npairs);
  CHECK(not compatibleEnergy or pairMassFluxPtr->size() == npairs);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, Wj, gWi, gWj, Pstar, rhostari, rhostarj, epsstari, epsstarj;
    Vector vstar;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvolDt_thread = DvolDt.threadCopy(threadStack);
    auto DmDt_thread = DmDt.threadCopy(threadStack);
    auto DEDt_thread = DEDt.threadCopy(threadStack);
    auto DpDt_thread = DpDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto newRiemannDpDx_thread = newRiemannDpDx.threadCopy(threadStack);
    auto newRiemannDvDx_thread = newRiemannDvDx.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      const auto& mi = mass(nodeListi, i);
      const auto& mj = mass(nodeListj, j);

      if( mi >tiny or mj > tiny){

      // Get the state for node i.
      const auto& riemannDpDxi = riemannDpDx(nodeListi, i);
      const auto& riemannDvDxi = riemannDvDx(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& voli = volume(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(voli > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DvolDti = DvolDt_thread(nodeListi,i);
      auto& DmDti = DmDt_thread(nodeListi, i);
      auto& DEDti = DEDt_thread(nodeListi, i);
      auto& DpDti = DpDt_thread(nodeListi, i);
      const auto& DxDti = DxDt(nodeListi,i);
      auto& newRiemannDpDxi = newRiemannDpDx_thread(nodeListi, i);
      auto& newRiemannDvDxi = newRiemannDvDx_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      const auto& gradRhoi = DrhoDx(nodeListi, i);
      const auto& gradEpsi = DepsDx(nodeListi, i);
      const auto& Mi = M(nodeListi,i);

      // Get the state for node j
      const auto& riemannDpDxj = riemannDpDx(nodeListj, j);
      const auto& riemannDvDxj = riemannDvDx(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& volj = volume(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(rhoj > 0.0);
      CHECK(volj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DvolDtj = DvolDt_thread(nodeListj,j);
      auto& DmDtj = DmDt_thread(nodeListj, j);
      auto& DEDtj = DEDt_thread(nodeListj, j);
      auto& DpDtj = DpDt_thread(nodeListj, j);
      const auto& DxDtj = DxDt(nodeListj,j);
      auto& newRiemannDpDxj = newRiemannDpDx_thread(nodeListj,j);
      auto& newRiemannDvDxj = newRiemannDvDx_thread(nodeListj,j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      const auto& gradRhoj = DrhoDx(nodeListj, j);
      const auto& gradEpsj = DepsDx(nodeListj, j);
      const auto& Mj = M(nodeListj,j);

      // Node displacement.
      const auto rij = ri - rj;
      const auto rhatij =rij.unitVector();
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);


      // Symmetrized kernel weight and gradient.
      //------------------------------------------------------
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      const auto gradPsii = voli * Mi.Transpose()*gradWi;
      const auto gradPsij = volj * Mj.Transpose()*gradWj;

      const auto Astar = voli*gradPsii + volj*gradPsij;

      // Determine an effective pressure including a term to fight the tensile instability.
      //const auto fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
      const auto fij = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
      const auto Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
      const auto Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
      const auto Peffi = Pi + Ri;
      const auto Peffj = Pj + Rj;

      // Reimann Solver and Fluxes
      //------------------------------------------------------
      // we'll clean this up when we have a gradient 
      // implementation we're in love with
      auto gradPi = riemannDpDxi;
      auto gradPj = riemannDpDxj;
      auto gradVi = riemannDvDxi;
      auto gradVj = riemannDvDxj;
      if (gradType==GradientType::SPHSameTimeGradient or
          gradType==GradientType::SPHUncorrectedGradient){
        gradPi = newRiemannDpDx(nodeListi,i);
        gradPj = newRiemannDpDx(nodeListj,j);
        gradVi = newRiemannDvDx(nodeListi,i);
        gradVj = newRiemannDvDx(nodeListj,j);
      }
      // need grad rho and grad eps
      riemannSolver.interfaceState(ri,           rj, 
                                   Hi,           Hj, 
                                   rhoi,         rhoj,
                                   epsi,         epsj, 
                                   ci,           cj, 
                                   Peffi,        Peffj, 
                                   vi,           vj, 
                                   gradRhoi,     gradRhoj,
                                   gradEpsi,     gradEpsj, // dumby for gradeps
                                   gradPi,       gradPj, 
                                   gradVi,       gradVj, 
                                   Pstar,     //output
                                   vstar,     //output
                                   rhostari,  //output
                                   rhostarj,  //output
                                   epsstari,  //output
                                   epsstarj); //output

      const auto fluxSwitch = (nodeListi==nodeListj ? 1.0 : 0.0);
      const auto vframe = (DxDti+DxDtj)*0.5; // EVENTUALLY ACCOUNT FOR Hi,Hj difference
      const auto vflux = vstar-vframe;
      const auto fluxTowardsNodei = vflux.dot(rhatij) > 0;
      const auto rhostar = (fluxTowardsNodei ? rhostarj : rhostari); // we'll need to fix these later
      const auto epsstar = (fluxTowardsNodei ? epsstarj : epsstari); // we'll need to fix these later

      const auto massFlux = fluxSwitch * rhostar * vflux.dot(Astar);
      const auto momentumFlux = massFlux * vstar;
      const auto energyFlux = massFlux * epsstar;

      // mass
      //------------------------------------------------------
      DmDti -= massFlux;
      DmDtj += massFlux;

      // momentum
      //------------------------------------------------------
      const auto deltaDvDt = Pstar*Astar + momentumFlux; 
      DpDti -= deltaDvDt;
      DpDtj += deltaDvDt;

      // energy
      //------------------------------------------------------
      const auto deltaDepsDti = Pstar*Astar.dot(vi-vstar) - energyFlux;
      const auto deltaDepsDtj = Pstar*Astar.dot(vstar-vj) + energyFlux;

      DEDti += deltaDepsDti;
      DEDtj += deltaDepsDtj;
     
      if(compatibleEnergy){
        (*pairMassFluxPtr)[kk] = massFlux;
        (*pairAccelerationsPtr)[kk] = deltaDvDt;
        (*pairDepsDtPtr)[kk][0] = deltaDepsDti;
        (*pairDepsDtPtr)[kk][1] = deltaDepsDtj;
      }

      // volume change based on nodal velocity
      //-----------------------------------------------------
      DvolDti -= (DxDti-DxDtj).dot(gradPsii);
      DvolDtj -= (DxDti-DxDtj).dot(gradPsij);

      // gradients
      //------------------------------------------------------
      const auto deltaDvDxi = 2.0*(vi-vstar).dyad(gradPsii);
      const auto deltaDvDxj = 2.0*(vstar-vj).dyad(gradPsij);

      // based on riemann soln
      DvDxi -= deltaDvDxi;
      DvDxj -= deltaDvDxj;

      // while we figure out what we want ...
      switch(gradType){ 
        case GradientType::RiemannGradient: // default grad based on riemann soln
          newRiemannDvDxi -= deltaDvDxi;
          newRiemannDvDxj -= deltaDvDxj;
          newRiemannDpDxi -= 2.0*(Pi-Pstar)*gradPsii;
          newRiemannDpDxj -= 2.0*(Pstar-Pj)*gradPsij;
          break;
        case GradientType::HydroAccelerationGradient: // based on hydro accel for DpDx
          newRiemannDvDxi -= deltaDvDxi;
          newRiemannDvDxj -= deltaDvDxj;
          newRiemannDpDxi += Pstar*Astar/voli;
          newRiemannDpDxj -= Pstar*Astar/volj;
          break;
        case GradientType::SPHGradient: // raw gradients
          newRiemannDvDxi -= (vi-vj).dyad(gradPsii);
          newRiemannDvDxj -= (vi-vj).dyad(gradPsij);
          newRiemannDpDxi -= (Pi-Pj)*gradPsii;
          newRiemannDpDxj -= (Pi-Pj)*gradPsij;
          break;
        case GradientType::MixedMethodGradient: // raw gradient for P riemann gradient for v
          newRiemannDvDxi -= deltaDvDxi;
          newRiemannDvDxj -= deltaDvDxj;
          newRiemannDpDxi -= (Pi-Pj)*gradPsii;
          newRiemannDpDxj -= (Pi-Pj)*gradPsij;
          break;       
        default:
          break;
        }

      } //if statement
    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  } // OpenMP parallel region


  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  ni = nodeList.numInternalNodes();

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& voli = volume(nodeListi,i);
      CHECK(voli > 0.0);

      auto& DvolDti = DvolDt(nodeListi, i);
      DvolDti *= voli;

      // If needed finish the total energy derivative.
      //if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);
    } // nodes loop
  }   // nodeLists loop
}     // eval derivs method 

//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFV<Dimension>::
firstDerivativesLoop(const typename Dimension::Scalar /*time*/,
                     const typename Dimension::Scalar dt,
                     const DataBase<Dimension>& dataBase,
                     const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const {

  const auto tiny = std::numeric_limits<Scalar>::epsilon();
  const auto nodeMotionCoeff = this->nodeMotionCoefficient();
  const auto cfl = this->cfl();

  const auto calcSpatialGradients =  (this->gradientType() == GradientType::SPHSameTimeGradient 
                                  or  this->gradientType() == GradientType::SPHUncorrectedGradient);
  const auto correctSpatialGradients = (this->gradientType() == GradientType::SPHSameTimeGradient);

  const auto nodeMotion = this->nodeMotionType();
  const auto xsphMotion = (nodeMotion == NodeMotionType::XSPH);
  const auto fickianMotion = (nodeMotion == NodeMotionType::Fickian or
                              nodeMotion == NodeMotionType::EulerianFickian);
  const auto lagrangianMotion = nodeMotion == NodeMotionType::Lagrangian;
  const auto noMotion = (nodeMotion == NodeMotionType::Eulerian or
                         nodeMotion == NodeMotionType::EulerianFickian);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();

  // kernel
  const auto& W = this->kernel();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  const auto  W0 = W(0.0, 1.0);

  // Get the state and derivative FieldLists. 
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);

  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDx = derivs.fields(GSPHFieldNames::densityGradient, Vector::zero);
  auto  DepsDx = derivs.fields(GSPHFieldNames::specificThermalEnergyGradient, Vector::zero);
  auto  newRiemannDpDx = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannPressureGradient,Vector::zero);
  auto  newRiemannDvDx = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannVelocityGradient,Tensor::zero);
  auto  normalization = derivs.fields(HydroFieldNames::normalization, 0.0);
  
  CHECK(M.size() == numNodeLists);
  CHECK(DrhoDx.size() == numNodeLists);
  CHECK(DepsDx.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(newRiemannDpDx.size() == numNodeLists);
  CHECK(newRiemannDvDx.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists)

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, Wj, gWi, gWj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto M_thread = M.threadCopy(threadStack);
    auto DrhoDx_thread = DrhoDx.threadCopy(threadStack);
    auto DepsDx_thread = DepsDx.threadCopy(threadStack);
    auto newRiemannDpDx_thread = newRiemannDpDx.threadCopy(threadStack);
    auto newRiemannDvDx_thread = newRiemannDvDx.threadCopy(threadStack);
    auto DxDt_thread = DxDt.threadCopy(threadStack);
    auto normalization_thread = normalization.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;
      
      // Get the state for node i.
      const auto& vi = velocity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& voli = volume(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(voli > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DxDti = DxDt_thread(nodeListi,i);
      auto& normi = normalization(nodeListi,i);
      auto& DrhoDxi = DrhoDx_thread(nodeListi, i);
      auto& DepsDxi = DepsDx_thread(nodeListi, i);
      auto& newRiemannDpDxi = newRiemannDpDx_thread(nodeListi, i);
      auto& newRiemannDvDxi = newRiemannDvDx_thread(nodeListi, i); 
      auto& Mi = M_thread(nodeListi, i);

      // Get the state for node j
      const auto& vj = velocity(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& volj = volume(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(volj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DxDtj = DxDt_thread(nodeListj,j);
      auto& normj = normalization(nodeListj,j);
      auto& DrhoDxj = DrhoDx_thread(nodeListj, j);
      auto& DepsDxj = DepsDx_thread(nodeListj, j);
      auto& newRiemannDpDxj = newRiemannDpDx_thread(nodeListj, j);
      auto& newRiemannDvDxj = newRiemannDvDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);

      const auto rij = ri - rj;

      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      const auto psii = voli*Wi;
      const auto psij = volj*Wj;
      const auto gradPsii = voli*gradWi;
      const auto gradPsij = volj*gradWj;

      // gradients
      Mi -= rij.dyad(gradPsii);
      Mj -= rij.dyad(gradPsij);
      
      DrhoDxi -= (rhoi - rhoj) * gradPsii;
      DrhoDxj -= (rhoi - rhoj) * gradPsij;
      DepsDxi -= (epsi - epsj) * gradPsii;
      DepsDxj -= (epsi - epsj) * gradPsij;

      if (calcSpatialGradients){
        newRiemannDpDxi -= (Pi-Pj)*gradPsii;
        newRiemannDpDxj -= (Pi-Pj)*gradPsij;

        newRiemannDvDxi -= (vi-vj).dyad(gradPsii);
        newRiemannDvDxj -= (vi-vj).dyad(gradPsij);
      }

      // node motion relative to fluid  (MAYBE REMOVE SOUND SPEEED CHECK)
      //-----------------------------------------------------------
      if (xsphMotion) {
        const auto cij = 0.5*(ci+cj);
        const auto wij = cij * safeInv(max(10*(vi-vj).magnitude(),cij));
        DxDti -= wij*psii*(vi-vj);
        DxDtj -= wij*psij*(vj-vi);
      }
      if(fickianMotion){
        const auto fi = FastMath::pow4(Wi/(Hdeti*WnPerh));
        const auto fj = FastMath::pow4(Wj/(Hdetj*WnPerh));
        DxDti -= (1.0+fi)*gradPsii;
        DxDtj += (1.0+fj)*gradPsij;
      }

      normi += psii;
      normj += psij;

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  
  // Finish up the spatial gradient and nodal motion calculations
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = M[nodeListi]->nodeList();
    const auto  ni = nodeList.numInternalNodes();
    const auto& fluidNodeList = **(dataBase.fluidNodeListBegin() + nodeListi);
    const auto  allowALE = (fluidNodeList.allowALE() and !lagrangianMotion);

    // finish up spatial gradients
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      
      
      const auto& voli = volume(nodeListi,i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      auto& Mi = M(nodeListi, i);
      auto& normi = normalization(nodeListi, i);
      const auto Mdeti = std::abs(Mi.Determinant());

      normi += voli*Hdeti*W0;

      
      const auto enoughNeighbors =  numNeighborsi > Dimension::pownu(2);
      const auto goodM =  (Mdeti > 1e-2 and enoughNeighbors);                   

      Mi = ( goodM ? Mi.Inverse() : Tensor::one);

      if (correctSpatialGradients){
        auto& newRiemannDpDxi = newRiemannDpDx(nodeListi, i);
        auto& newRiemannDvDxi = newRiemannDvDx(nodeListi, i);
        auto& DrhoDxi = DrhoDx(nodeListi,i);
        auto& DepsDxi = DepsDx(nodeListi,i);

        DrhoDxi = Mi.Transpose()*DrhoDxi;
        DepsDxi = Mi.Transpose()*DepsDxi;
        newRiemannDpDxi = Mi.Transpose()*newRiemannDpDxi;
        newRiemannDvDxi = newRiemannDvDxi*Mi;
      }

      // if we are pure lagrangian do it here
      if (! allowALE) DxDti = vi;

    }

    // finish up the node motion
    if(allowALE){
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {

        const auto& ci = soundSpeed(nodeListi,i);
        const auto& vi = velocity(nodeListi,i);

        const auto& normi = normalization(nodeListi, i);
        auto& DxDti = DxDt(nodeListi,i);

        const auto nodeMotionLimiter = min(cfl*ci*safeInv(DxDti.magnitude()),1.0)
        if(xsphMotion) DxDti *= nodeMotionLimiter*nodeMotionCoeff/max(tiny, normi);
        if(fickianMotion) DxDti *= nodeMotionLimiter*nodeMotionCoeff*ci*ci*dt;
        if(!noMotion) DxDti += vi;
      }
    }
    
  }
  
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr){
          (*boundItr)->applyFieldListGhostBoundary(M);
          (*boundItr)->applyFieldListGhostBoundary(DxDt);
  }

  if (calcSpatialGradients){ 
    for (ConstBoundaryIterator boundItr = this->boundaryBegin();
          boundItr != this->boundaryEnd();
           ++boundItr){
      (*boundItr)->applyFieldListGhostBoundary(DepsDx);
      (*boundItr)->applyFieldListGhostBoundary(DrhoDx);
      (*boundItr)->applyFieldListGhostBoundary(newRiemannDpDx);
      (*boundItr)->applyFieldListGhostBoundary(newRiemannDvDx);
    }
  }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
 
}

} // spheral namespace
