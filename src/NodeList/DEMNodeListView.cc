//---------------------------------Spheral++----------------------------------//
// DEMNodeListView -- allocation-free view/base state for DEMNodeList.
//----------------------------------------------------------------------------//
#include "NodeList/DEMNodeList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// when problem starts set our equilibrium overlap
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
void
DEMNodeListView<Dimension>::
setHfieldFromParticleRadius(const size_t startUniqueIndex,
                            const Scalar kernelExtent){
  
  const auto& radius = this->particleRadius();
  const auto& uniqId = this->uniqueIndex();
  auto& Hfield = this->Hfield();
  const auto ni = this->numInternalNodes();

#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    const auto ui = uniqId[i];
    if(ui >= startUniqueIndex) {
      const auto Ri = radius[i];
      auto& Hi = Hfield[i];
      const auto hInv = safeInv(2.0 * Ri * (1.0+mNeighborSearchBuffer)/kernelExtent);
      Hi = SymTensor::one() * hInv;
    }
  }
}

}
