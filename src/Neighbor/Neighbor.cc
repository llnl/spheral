//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
// Split into Value/View pattern Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#include "Neighbor/Neighbor.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
Neighbor(DataBase<Dimension>& dataBase,
         const NeighborSearchType searchType,
         const Scalar kernelExtent):
  mDataBase(dataBase),
  mNodeExtent(FieldStorageType::CopyFields),
  mSearchType(searchType),
  mKernelExtent(kernelExtent) {
}

}
