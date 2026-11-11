//------------------------------------------------------------------------------
// Append elements to a container. Grow resizable containers when necessary;
// fixed-capacity views report failure so their owning caller can resize and retry.
//------------------------------------------------------------------------------
#ifndef __Spheral_copy_n_conditional__
#define __Spheral_copy_n_conditional__

#include "config.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
#include <cstddef>

namespace Spheral {

template<typename InContainer,
         typename OutContainer>
SPHERAL_HOST_DEVICE
inline
bool
copy_n_conditional(const InContainer& source,
                   const size_t numSource,
                   OutContainer& dest,
                   const size_t sourceOffset = 0u,
                   const size_t destOffset = 0u) {
  REQUIRE(sourceOffset <= source.size());
  REQUIRE(numSource <= source.size() - sourceOffset);
  if (numSource == 0u) return true;
  const auto sizeAvail = destOffset < dest.size() ? dest.size() - destOffset : 0u;
  if (sizeAvail >= numSource) {
    std::copy_n(source.begin() + sourceOffset, numSource, dest.begin() + destOffset);
    return true;
  }

  if constexpr (requires { dest.resize(std::size_t{0}); }) {
    dest.resize(destOffset + numSource);
    std::copy_n(source.begin() + sourceOffset, numSource, dest.begin() + destOffset);
    return true;
  }

  return false;
}

}

#endif
