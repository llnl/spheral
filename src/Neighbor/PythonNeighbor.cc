//---------------------------------Spheral++----------------------------------//
// PythonNeighbor
//
// PythonNeighbor just provides the virtual function hooks to support the
// NeighborView interface.  The virutal Neighbor interface is passed along
// for descendants to provide.
//
// NOTE: PythonNeighbor and its descendants are CPU-only (no GPU capability).
//----------------------------------------------------------------------------//

#include "Neighbor/PythonNeighbor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
PythonNeighbor<Dimension>::
PythonNeighbor(DataBase<Dimension>& dataBase,
               const NeighborSearchType searchType,
               const Scalar kernelExtent):
  Neighbor<Dimension>(dataBase, searchType, kernelExtent),
  NeighborView<Dimension, PythonNeighbor<Dimension>>(dataBase.view(), searchType, kernelExtent) {
}

}
