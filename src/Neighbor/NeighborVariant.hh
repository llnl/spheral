//---------------------------------Spheral++----------------------------------//
// NeighborVariants -- The variant type of possible Neighbor objects
//
// Created by J. Michael Owen, Mon Sep 28 16:53:54 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborVariant__
#define __Spheral_NeighborVariant__

#include <memory>
#include <variant>

namespace Spheral {

template<typename Dimension> class TreeNeighbor;
template<typename Dimension> class TreeNeighborView;

template<typename Dimension> class PythonNeighbor;

//------------------------------------------------------------------------------
// NeighborViews: used by clients to get potential Neighbors (possibly on device)
//------------------------------------------------------------------------------
template<typename Dimension>
using NeighborViewVariant = std::variant<std::unique_ptr<TreeNeighborView<Dimension>>,
                                         std::unique_ptr<PythonNeighbor<Dimension>>>;

//------------------------------------------------------------------------------
// NeighorValues: the full Neighbor objects that allocate and hold memory
//------------------------------------------------------------------------------
template<typename Dimension>
using NeighborVariant = std::variant<std::unique_ptr<TreeNeighbor<Dimension>>,
                                     std::unique_ptr<PythonNeighbor<Dimension>>>;

//------------------------------------------------------------------------------
// Visitor pattern
//------------------------------------------------------------------------------
template<typename Variant, typename Visitor>
decltype(auto)
visitNeighbor(Variant&& neighborPtr,
              Visitor&& visitor) {
  return std::visit(
    [&](const auto& neighborPtr) -> decltype(auto) {
      REQUIRE(neighborPtr);
      return std::forward<Visitor>(visitor)(*neighborPtr);
    },
    std::forward<Variant>(neighborPtr));
}

}

#include "Neighbor/PythonNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"

#endif
