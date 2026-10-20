//---------------------------------Spheral++----------------------------------//
// NeighborVariants -- The variant type of possible Neighbor objects
//
// Created by J. Michael Owen, Mon Sep 28 16:53:54 PDT 2026
//----------------------------------------------------------------------------//

#include "Neighbor/PythonNeighbor.hh"

#include <variant>

namespace Spheral {

// NeighborViews: used by clients to get potential Neighbors (possibly on device)
template<typename Dimension>
using NeighborViewVariant = std::variant<PythonNeighbor<Dimension>>;

// NeighorValues: the full Neighbor objects that allocate and hold memory
template<typename Dimension>
using NeighborValueVariant = std::variant<PythonNeighbor<Dimension>>;

}
