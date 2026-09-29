//---------------------------------Spheral++----------------------------------//
// PythonNeighbor
//
// PythonNeighbor just provides the virtual function hooks to support the
// NeighborView interface.
//
// NOTE: PythonNeighbor and its descendants are CPU-only (no GPU capability).
//
// Created by J. Michael Owen, Thu Sep 24 15:41:42 PDT 2026
//----------------------------------------------------------------------------//

#include "Neighbor/PythonNeighbor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
PythonNeighbor<Dimension>::
PythonNeighbor(NodeListView<Dimension>& nodeListView,
                   FieldView<Dimension, Vector>& nodeExtentView,
                   const NeighborSearchType searchType,
                   const Scalar kernelExtent):
  Neighbor<Dimension, PythonNeighbor<Dimension>>(nodeListView, nodeExtentView, searchType, kernelExtent) {
}

}
