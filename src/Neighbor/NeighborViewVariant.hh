//---------------------------------Spheral++----------------------------------//
// NeighborViewVariant -- Non-owning handles to supported Neighbor views.
//
// Keep this header limited to forward declarations so DataBaseView can contain
// the variant without including the concrete NeighborView definitions.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborViewVariant__
#define __Spheral_NeighborViewVariant__

#include "Utilities/DBC.hh"

#include <variant>
#include <utility>

namespace Spheral {

template<typename Dimension> class TreeNeighborView;
template<typename Dimension> class PythonNeighbor;

template<typename Dimension>
using NeighborViewVariant = std::variant<TreeNeighborView<Dimension>*
#ifndef SPHERAL_ENABLE_GPU
                                       , PythonNeighbor<Dimension>*
#endif
                                        >;

//------------------------------------------------------------------------------
// Visitor pattern
//------------------------------------------------------------------------------
template<typename Variant, typename Visitor>
SPHERAL_HOST_DEVICE 
decltype(auto)
visitNeighborView(Variant&& neighborPtr,
                  Visitor&& visitor) {
  return std::visit(
    [&](const auto& neighborPtr) -> decltype(auto) {
      REQUIRE(neighborPtr);
      return std::forward<Visitor>(visitor)(*neighborPtr);
    },
    std::forward<Variant>(neighborPtr));
}

}

#endif
