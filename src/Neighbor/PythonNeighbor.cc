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
#include "DataBase/DataBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
PythonNeighbor<Dimension>::
PythonNeighbor(DataBase<Dimension>& dataBase,
               const NeighborSearchType searchType,
               const Scalar kernelExtent):
  Neighbor<Dimension>(dataBase),
  NeighborView<Dimension, PythonNeighbor<Dimension>>(dataBase.view(), searchType, kernelExtent) {
}

//------------------------------------------------------------------------------
// Force the nodeExtent FieldList to be sized correctly
//------------------------------------------------------------------------------
template<typename Dimension>
void
PythonNeighbor<Dimension>::
resizeNodeExtents() {
  mDataBase.resizeGlobalFieldList(mNodeExtent,
                                  Vector::zero(),
                                  "node extents",
                                  false);
  mNodeExtentView = mNodeExtent.view();
}

//------------------------------------------------------------------------------
// Clear our state
//------------------------------------------------------------------------------
template<typename Dimension>
void
PythonNeighbor<Dimension>::
clear() {
  resizeNodeExtents();
}

}
