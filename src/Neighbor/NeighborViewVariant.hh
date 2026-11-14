//---------------------------------Spheral++----------------------------------//
// NeighborViewVariant -- The closed set of Neighbor views.
//
// Visit this variant on the host, then capture the selected value view in a
// device kernel.  GPU builds contain only device-callable value views;
// CPU-only builds additionally contain the Python callback adapter.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborViewVariant__
#define __Spheral_NeighborViewVariant__

#include "Neighbor/TreeNeighborView.hh"
#include "Utilities/DBC.hh"
#ifndef SPHERAL_ENABLE_GPU
#include "Neighbor/PythonNeighborCallView.hh"
#endif

#include <type_traits>
#include <variant>
#include <utility>

namespace Spheral {

#ifdef SPHERAL_ENABLE_GPU
template<typename Dimension>
using NeighborViewVariant = std::variant<std::monostate,
                                         TreeNeighborView<Dimension>>;
#else
template<typename Dimension>
using NeighborViewVariant = std::variant<std::monostate,
                                         TreeNeighborView<Dimension>,
                                         PythonNeighborCallView<Dimension>>;
#endif

//------------------------------------------------------------------------------
// Visitor pattern
//------------------------------------------------------------------------------
template<typename Variant, typename Visitor>
decltype(auto)
visitNeighborView(Variant&& neighborVariant,
                  Visitor&& visitor) {
  return std::visit(
    [&](const auto& neighborView) -> decltype(auto) {
      using ViewType = std::decay_t<decltype(neighborView)>;
      if constexpr (std::is_same_v<ViewType, std::monostate>) {
        VERIFY2(false, "Attempt to visit an empty Neighbor view.");
      } else {
        return std::forward<Visitor>(visitor)(neighborView);
      }
    },
    std::forward<Variant>(neighborVariant));
}

}

#endif
