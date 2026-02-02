//---------------------------------Spheral++----------------------------------//
// Compute the density from m/V
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/computeMFMDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

#include <limits.h>

namespace Spheral{

template<typename Dimension>
void
computeMFMDensity(const FieldList<Dimension, typename Dimension::Scalar>& mass,
                  const FieldList<Dimension, typename Dimension::Scalar>& volume,
                        FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = volume.size();
  REQUIRE(massDensity.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = volume[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      massDensity(nodeListi,i) =  mass(nodeListi,i)/std::max(volume(nodeListi,i),tiny);
    }
  }   
}     // function
}     // spheral namespace
