//---------------------------------Spheral++----------------------------------//
// globalNodeIDs (contains numGlobalNodes and globalNodeIDs)
//
// Helper methods to assign unique global node IDs to all node in a NodeList,
// serial or parallel.  *Should* compute the same unique ID for each node no
// matter how the parallel domains are carved up.  
//
// Note this is accomplished by sorting based on position, so nodes occupying
// the same position are considered degenerate and may get different IDs for
// different decompositions.
//
// Created by JMO, Fri Oct 29 13:08:33 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_globalNodeIDs_hh__
#define __Spheral_globalNodeIDs_hh__

namespace Spheral {

template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;

//------------------------------------------------------------------------------
// Return the total global number of nodes in the NodeList (or DataBase).
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
numGlobalNodes(const NodeList<Dimension>& nodeList);
  
template<typename Dimension>
size_t
numGlobalNodes(const DataBase<Dimension>& dataBase);
  
template<typename Dimension, typename NodeListIterator>
size_t
numGlobalNodes(const NodeListIterator& begin,
               const NodeListIterator& end);
  
//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for the given NodeList, and return
// the set of them on this process.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, size_t>
globalNodeIDs(const NodeList<Dimension>& nodeList);

//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes across all NodeLists in
// a DataBase, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, size_t>
globalNodeIDs(const DataBase<Dimension>& dataBase);

template<typename Dimension, typename NodeListIterator>
FieldList<Dimension, size_t>
globalNodeIDs(const NodeListIterator& begin,
              const NodeListIterator& end);

}

#include "globalNodeIDsInline.hh"

#endif
