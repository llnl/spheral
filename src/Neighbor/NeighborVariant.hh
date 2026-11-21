//---------------------------------Spheral++----------------------------------//
// NeighborVariant -- The variant types of possible Neighbor objects
//
// Created by J. Michael Owen, Mon Sep 28 16:53:54 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborVariant__
#define __Spheral_NeighborVariant__

#include <utility>
#include <variant>

#include "Neighbor/TreeNeighbor.hh"
#ifndef SPHERAL_ENABLE_GPU
#include "Neighbor/PythonNeighborCallView.hh"
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Neighbor dispatch objects.  Native neighbors are stored by value; CPU-only
// Python neighbors are represented by a lightweight, non-owning call view.
//------------------------------------------------------------------------------
#ifdef SPHERAL_ENABLE_GPU
template<typename Dimension>
using NeighborVariant = std::variant<TreeNeighbor<Dimension>>;
#else
template<typename Dimension>
using NeighborVariant = std::variant<TreeNeighbor<Dimension>,
                                     PythonNeighborCallView<Dimension>>;
#endif

//------------------------------------------------------------------------------
// Visitor pattern
//------------------------------------------------------------------------------
template<typename Variant, typename Visitor>
decltype(auto)
visitNeighbor(Variant&& neighborVariant,
              Visitor&& visitor) {
  return std::visit(
    [&](auto&& neighbor) -> decltype(auto) {
      return std::forward<Visitor>(visitor)(
        std::forward<decltype(neighbor)>(neighbor));
    },
    std::forward<Variant>(neighborVariant));
}

}

#endif
