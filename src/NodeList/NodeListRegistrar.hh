//------------------------------------------------------------------------------
// NodeListRegistrar
// 
// A singleton class that maintains a running (sorted) record of NodeList names 
// that are currently in scope.  This is necessary to ensure that there is a
// globally available unique order to iterate over NodeLists (and Fields in 
// FieldLists) that will be uniform across processors.
//
// Created by JMO, Fri Aug  5 10:16:25 PDT 2005
//------------------------------------------------------------------------------
#ifndef __Spheral_NodeListRegistrar__
#define __Spheral_NodeListRegistrar__

#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Field/FieldBase.hh"
#include "Field/Field.hh"

#include <string>
#include <vector>

namespace Spheral {

// // Forward decalarations.
// template<typename Dimension> class NodeList;
// template<typename Dimension> class FluidNodeList;
// template<typename Dimension> class FieldBase;

template<typename Dimension>
class NodeListRegistrar {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef std::vector<NodeList<Dimension>*> ContainerType;
  typedef std::vector<FluidNodeList<Dimension>*> FluidContainerType;

  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;

  typedef typename FluidContainerType::iterator fluid_iterator;
  typedef typename FluidContainerType::const_iterator const_fluid_iterator;

  // Define a nested comparator class for comparing NodeLists by their names.
  struct NodeListComparator {
    int operator()(const NodeList<Dimension>* nodeListPtr1,
                   const NodeList<Dimension>* nodeListPtr2) const {
      return nodeListPtr1->name() < nodeListPtr2->name();
    }
    int operator()(const FieldBase<Dimension>* fieldPtr1,
                   const FieldBase<Dimension>* fieldPtr2) const {
      return fieldPtr1->nodeListPtr()->name() < fieldPtr2->nodeListPtr()->name();
    }
  };

  // Get the instance.
  static NodeListRegistrar& instance();

  // The number of NodeLists currently in play.
  int numNodeLists() const;
  int numFluidNodeLists() const;

  // Iterators over the current NodeLists (in proper order of course).
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  fluid_iterator fluidBegin();
  fluid_iterator fluidEnd();

  const_fluid_iterator fluidBegin() const;
  const_fluid_iterator fluidEnd() const;

  // The (sorted) set of registered NodeList names.
  std::vector<std::string> registeredNames() const;
  std::vector<std::string> registeredFluidNames() const;

  // Internal consistency checking.
  bool valid() const;

  // Flag (for use by others) indicating whether we want to run in 
  // strictly domain decomposition independent/reproducing mode (may have
  // a performance impact).
  bool domainDecompositionIndependent() const;
  void domainDecompositionIndependent(const bool x);

  //---------------------------------------------------------------------------
  // Static methods
  // Determine the proper place in a sequence of Fields that a given Field
  // should be inserted.
  // This is primarily a helper for building FieldLists properly.
  template<typename IteratorType, typename ThingyType>
  static
  IteratorType
  findInsertionPoint(const ThingyType& thingy,
                     const IteratorType begin,
                     const IteratorType end);

  // Sort the given iterator range into NodeList order
  template<typename IteratorType>
  static
  void
  sortInNodeListOrder(IteratorType begin, IteratorType end);

  // Helpers for findInsertionPoint that know how to extract a NodeList*
  // from the argument.
  static
  NodeList<Dimension>* 
  getNodeListPtr(NodeList<Dimension>* thingy) {
    return thingy;
  }

  static
  NodeList<Dimension>* 
  getNodeListPtr(const NodeList<Dimension>* thingy) {
    return const_cast<NodeList<Dimension>*>(thingy);
  }

  static
  NodeList<Dimension>* 
  getNodeListPtr(FluidNodeList<Dimension>* thingy) {
    return (NodeList<Dimension>*) thingy;
  }

  static
  NodeList<Dimension>* 
  getNodeListPtr(FieldBase<Dimension>* thingy) {
    return const_cast<NodeList<Dimension>*>(thingy->nodeListPtr());
  }

  static
  NodeList<Dimension>* 
  getNodeListPtr(const FieldBase<Dimension>* thingy) {
    return const_cast<NodeList<Dimension>*>(thingy->nodeListPtr());
  }

private:
  //--------------------------- Private Interface---------------------------//

  // The current set of NodeLists.
  ContainerType mNodeLists;
  FluidContainerType mFluidNodeLists;

  // Flag for the domain independent choice.
  bool mDomainDecompIndependent;

  // No public constructors, destructor, or assignment.
  NodeListRegistrar();
  ~NodeListRegistrar();

  // Grant friendship to NodeLists to register and unregister themselves.
  friend class NodeList<Dimension>;
  friend class FluidNodeList<Dimension>;
  void registerNodeList(NodeList<Dimension>& nodeList);
  void registerNodeList(FluidNodeList<Dimension>& nodeList);
  void unregisterNodeList(NodeList<Dimension>& nodeList);
  void unregisterNodeList(FluidNodeList<Dimension>& nodeList);

public:
  NodeListRegistrar(const NodeListRegistrar&) = delete;
  NodeListRegistrar& operator=(const NodeListRegistrar&) = delete;

};

}

#include "NodeListRegistrarInline.hh"

#endif
