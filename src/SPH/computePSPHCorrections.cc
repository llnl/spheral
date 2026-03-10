//------------------------------------------------------------------------------
// Compute the PSPH grad h corrections due to Hopkins 2013
//------------------------------------------------------------------------------
#include "computePSPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
computePSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& mass,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
                       const FieldList<Dimension, typename Dimension::Scalar>& gamma,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       const bool computeMassDensity,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHmassDensity,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHpbar,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHsoundSpeed,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHcorrection) {

  // Pre-conditions.
  const size_t numNodeLists = PSPHpbar.size();
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(specificThermalEnergy.size() == numNodeLists);
  REQUIRE(gamma.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(PSPHmassDensity.size() == numNodeLists);
  REQUIRE(PSPHsoundSpeed.size() == numNodeLists);
  REQUIRE(PSPHcorrection.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Zero out the result.
  if (computeMassDensity) PSPHmassDensity = 0.0;
  PSPHpbar = 0.0;
  PSPHcorrection = 0.0;

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);
  const double tiny = 1.0e-30;

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  FieldList<Dimension, Scalar> Nbar(FieldStorageType::CopyFields);      //Averaged Number of particles 
  FieldList<Dimension, Scalar> gradPbar(FieldStorageType::CopyFields);  //DpbarDh
  FieldList<Dimension, Scalar> gradNbar(FieldStorageType::CopyFields);  //DnbarDh
  for (const auto& fieldPtr: mass) {
    Nbar.appendNewField("Nbar", fieldPtr->nodeList(), 0.0);
    gradPbar.appendNewField("gradPbar", fieldPtr->nodeList(), 0.0);
    gradNbar.appendNewField("gradNbar", fieldPtr->nodeList(), 0.0);
  }

  // Sum the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    Scalar Wi, Wj, gWi, gWj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto PSPHmassDensity_thread = PSPHmassDensity.threadCopy(threadStack);
    auto PSPHpbar_thread = PSPHpbar.threadCopy(threadStack);
    auto PSPHcorrection_thread = PSPHcorrection.threadCopy(threadStack);
    auto Nbar_thread = Nbar.threadCopy(threadStack);
    auto gradPbar_thread = gradPbar.threadCopy(threadStack);
    auto gradNbar_thread = gradNbar.threadCopy(threadStack);

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  mi = mass(nodeListi, i);
      const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  gammai = gamma(nodeListi, i);
      const auto  invhi = (Hi.Trace()/Dimension::nDim);
      const auto  xi = (gammai-1.0)*mi*epsi;
     
      // Get the state for node j.
      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  mj = mass(nodeListj, j);
      const auto  invhj = (Hj.Trace()/Dimension::nDim);
      const auto  epsj = specificThermalEnergy(nodeListj, j);
      const auto  gammaj = gamma(nodeListj, j);
      const auto  xj = (gammaj-1.0)*mj*epsj;

      // Kernel weighting and gradient.
      const Vector rij = ri - rj;
      const Scalar etai = (Hi*rij).magnitude();
      const Scalar etaj = (Hj*rij).magnitude();
      W.kernelAndGradValue(etai, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaj, Hdetj, Wj, gWj);

      const auto gradhi = invhi*(Dimension::nDim*Wi+etai*gWi);
      const auto gradhj = invhj*(Dimension::nDim*Wj+etaj*gWj);

      if (computeMassDensity) {
        PSPHmassDensity_thread(nodeListi, i) += (nodeListi == nodeListj ? mj : mi)*Wj;
        PSPHmassDensity_thread(nodeListj, j) += (nodeListi == nodeListj ? mi : mj)*Wi;
      }

      PSPHpbar_thread(nodeListi, i) += xj*Wi;
      PSPHpbar_thread(nodeListj, j) += xi*Wj;

      Nbar_thread(nodeListi, i) += Wi;
      Nbar_thread(nodeListj, j) += Wj;
      
      gradPbar_thread(nodeListi, i) -= xj*gradhi;
      gradPbar_thread(nodeListj, j) -= xi*gradhj;

      gradNbar_thread(nodeListi, i) -= gradhi;
      gradNbar_thread(nodeListj, j) -= gradhj;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OMP parallel

  // Finish with the self contributions.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {

      // Get the state for node i.
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar mi = mass(nodeListi, i);
      const Scalar invhi = (Hi.Trace()/Dimension::nDim);
      const Scalar epsi = specificThermalEnergy(nodeListi, i);
      const Scalar gammai = gamma(nodeListi, i);
     
      const Scalar xi = (gammai-1.0)*mi*epsi;
      const Scalar gradh0 = invhi*(Dimension::nDim*Hdeti*W0);
      if (computeMassDensity) PSPHmassDensity(nodeListi, i) += mi*Hdeti*W0;
      PSPHpbar(nodeListi, i) += xi*Hdeti*W0;
      Nbar(nodeListi, i) += Hdeti*W0;
      gradPbar(nodeListi, i) += -xi*gradh0;
      gradNbar(nodeListi, i) += -gradh0;

      //const Scalar fi=1.0+gradNbari*safeInv(Dimension::nDim*Nbari*invhi);
      //PSPHcorrection(nodeListi, i)=gradPbari*safeInv(Dimension::nDim*(gammai-1.0)*Nbari*invhi*fi);
      const Scalar fi = 1.0 + gradNbar(nodeListi, i)/max(Dimension::nDim*Nbar(nodeListi, i)*invhi, tiny);
      PSPHcorrection(nodeListi, i) += gradPbar(nodeListi, i)/max(Dimension::nDim*(gammai-1.0)*Nbar(nodeListi, i)*invhi*fi, tiny);
      CHECK2((gammai-1.0)*epsi >= 0.0, i << " " << gammai << " " << epsi);
      PSPHsoundSpeed(nodeListi, i) = sqrt(std::max(0.0, gammai*(gammai - 1.0)*epsi));
    }
  }
}

}
