//---------------------------------Spheral++----------------------------------//
// nodeBoundingBoxes
//
// Compute minimum volume bounding boxes for nodes and their extent.
//
// Created by JMO, Sun Jan 24 16:13:16 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_nodeBoundingBox__
#define __Spheral_nodeBoundingBox__

#include <utility>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class NodeList;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension> class DataBase;

//------------------------------------------------------------------------------
// The bounding box for a position and H.
//------------------------------------------------------------------------------
template<typename Dimension>
std::pair<typename Dimension::Vector, typename Dimension::Vector>
boundingBox(const typename Dimension::Vector& xi,
            const typename Dimension::SymTensor& Hi,
            const typename Dimension::Scalar& kernelExtent);

//------------------------------------------------------------------------------
// The bounding boxes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, std::pair<typename Dimension::Vector, typename Dimension::Vector> >
nodeBoundingBoxes(const NodeList<Dimension>& nodes,
                  const typename Dimension::Scalar kernelExtent);

//------------------------------------------------------------------------------
// The bounding boxes for all nodes in a DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, std::pair<typename Dimension::Vector, typename Dimension::Vector> >
nodeBoundingBoxes(const DataBase<Dimension>& dataBase);

}

#include "nodeBoundingBoxesInline.hh"

#endif
