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
#include "NodeListRegistrar.hh"
#include "Geometry/Dimension.hh"

#include <sstream>

namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
template<typename Dimension>
NodeListRegistrar<Dimension>&
NodeListRegistrar<Dimension>::
instance() {
  static NodeListRegistrar theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// The (sorted) set of registered NodeList names.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<std::string>
NodeListRegistrar<Dimension>::
registeredNames() const {
  std::vector<std::string> result;
  result.reserve(mNodeLists.size());
  for (const_iterator itr = mNodeLists.begin();
       itr != mNodeLists.end();
       ++itr) result.push_back((**itr).name());
  ENSURE(result.size() == mNodeLists.size());
  return result;
}

template<typename Dimension>
std::vector<std::string>
NodeListRegistrar<Dimension>::
registeredFluidNames() const {
  std::vector<std::string> result;
  result.reserve(mFluidNodeLists.size());
  for (const_fluid_iterator itr = mFluidNodeLists.begin();
       itr != mFluidNodeLists.end();
       ++itr) result.push_back((**itr).name());
  ENSURE(result.size() == mFluidNodeLists.size());
  return result;
}

//------------------------------------------------------------------------------
// Check that we're internally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NodeListRegistrar<Dimension>::
valid() const {
  // Pretty much just make sure that the registered NodeLists are properly
  // sorted by name.
  bool result = (mFluidNodeLists.size() <= mNodeLists.size());
  if (numNodeLists() > 0) {
    const_iterator itr = begin();
    while (result && itr < (end() - 1)) {
      result = ((**itr).name() < (**(itr + 1)).name());
      ++itr;
    }
  }
  if (numFluidNodeLists() > 0) {
    const_fluid_iterator itr = fluidBegin();
    while (result && itr < (fluidEnd() - 1)) {
      result = ((**itr).name() < (**(itr + 1)).name());
      ++itr;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Constructor (private).
//------------------------------------------------------------------------------
template<typename Dimension>
NodeListRegistrar<Dimension>::
NodeListRegistrar():
  mNodeLists(),
  mFluidNodeLists(),
  mDomainDecompIndependent(false) {
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
template<typename Dimension>
NodeListRegistrar<Dimension>::
~NodeListRegistrar() {
}

//------------------------------------------------------------------------------
// Ensure that both the name and its hashed identity are unique.  NodeListView
// comparisons use the hash, so admitting a collision would make two distinct
// NodeLists compare equal.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeListRegistrar<Dimension>::
verifyUniqueNameAndHash(const NodeList<Dimension>& nodeList) const {
  const auto newName = nodeList.name();
  const auto newHash = nodeList.hash();
  for (const auto* registeredNodeList: mNodeLists) {
    const auto registeredName = registeredNodeList->name();
    if (registeredName == newName) {
      std::stringstream msg;
      msg << "NodeListRegistrar ERROR: the name " << newName << '\n'
          << "is already registered.";
      VERIFY2(false, msg.str());
    }
    if (registeredNodeList->hash() == newHash) {
      std::stringstream msg;
      msg << "NodeListRegistrar ERROR: hash collision between the NodeList names "
          << newName << " and " << registeredName
          << " (hash " << newHash << ").";
      VERIFY2(false, msg.str());
    }
  }
}

//------------------------------------------------------------------------------
// Register a NodeList (private).
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeListRegistrar<Dimension>::
registerNodeList(NodeList<Dimension>& nodeList) {

  // Make sure this isn't a NodeList we already have registered.
  VERIFY2(std::find(mNodeLists.begin(), mNodeLists.end(), &nodeList) == mNodeLists.end(),
          "NodeListRegistrar ERROR: Attempt to register a NodeList we already have!");

  verifyUniqueNameAndHash(nodeList);

  // Find where we should insert the new NodeList in the ordering, and 
  // insert it.
  iterator itr = std::upper_bound(mNodeLists.begin(), mNodeLists.end(), &nodeList, NodeListComparator());
  mNodeLists.insert(itr, &nodeList);

  // Post-conditions.
  ENSURE(std::find(mNodeLists.begin(), mNodeLists.end(), &nodeList) != mNodeLists.end());
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Register a FluidNodeList (private).
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeListRegistrar<Dimension>::
registerNodeList(FluidNodeList<Dimension>& nodeList) {

  // Make sure this isn't a NodeList we already have registered.
  VERIFY2(std::find(mFluidNodeLists.begin(), mFluidNodeLists.end(), &nodeList) == mFluidNodeLists.end(),
          "NodeListRegistrar ERROR: Attempt to register a FluidNodeList we already have!");

  verifyUniqueNameAndHash(nodeList);

  // Find where we should insert the new NodeList in the ordering, and 
  // insert it.
  {
    iterator itr = std::upper_bound(mNodeLists.begin(), mNodeLists.end(), (NodeList<Dimension>*) &nodeList,
                                    NodeListComparator());
    mNodeLists.insert(itr, (NodeList<Dimension>*) &nodeList);
  }
  {
    fluid_iterator itr = std::upper_bound(mFluidNodeLists.begin(), mFluidNodeLists.end(), &nodeList,
                                          NodeListComparator());
    mFluidNodeLists.insert(itr, &nodeList);
  }

  // Post-conditions.
  ENSURE(std::find(mNodeLists.begin(), mNodeLists.end(), (NodeList<Dimension>*) &nodeList) != mNodeLists.end());
  ENSURE(std::find(mFluidNodeLists.begin(), mFluidNodeLists.end(), &nodeList) != mFluidNodeLists.end());
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Unregister a NodeList (private).
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeListRegistrar<Dimension>::
unregisterNodeList(NodeList<Dimension>& nodeList) {

  // Make sure this is a NodeList we have registered.
  iterator itr = std::find(mNodeLists.begin(), mNodeLists.end(), &nodeList);
  VERIFY2(itr != mNodeLists.end(),
          "NodeListRegistrar ERROR: Attempt to unregister a NodeList we don't have!");

  // Remove the thing and we're done.
  mNodeLists.erase(itr);

  // Post-conditions.
  ENSURE(std::find(mNodeLists.begin(), mNodeLists.end(), &nodeList) == mNodeLists.end());
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Unregister a FluidNodeList (private).
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeListRegistrar<Dimension>::
unregisterNodeList(FluidNodeList<Dimension>& nodeList) {

  // Make sure this is a NodeList we have registered.
  {
    iterator itr = std::find(mNodeLists.begin(), mNodeLists.end(), (NodeList<Dimension>*) &nodeList);
    VERIFY2(itr != mNodeLists.end(),
            "NodeListRegistrar ERROR: Attempt to unregister a NodeList we don't have!");
    mNodeLists.erase(itr);
  }

  // Remove it from the set of FluidNodeLists.
  {
    fluid_iterator itr = std::find(mFluidNodeLists.begin(), mFluidNodeLists.end(), &nodeList);
    VERIFY2(itr != mFluidNodeLists.end(),
            "NodeListRegistrar ERROR: Attempt to unregister a NodeList we don't have!");
    mFluidNodeLists.erase(itr);
  }

  // Post-conditions.
  ENSURE(std::find(mNodeLists.begin(), mNodeLists.end(), (NodeList<Dimension>*) &nodeList) == mNodeLists.end());
  ENSURE(std::find(mFluidNodeLists.begin(), mFluidNodeLists.end(), &nodeList) == mFluidNodeLists.end());
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Flag for whether we should try to run in a domain decomposition independent/
// reproducing mode.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NodeListRegistrar<Dimension>::
domainDecompositionIndependent() const {
  return mDomainDecompIndependent;
}

template<typename Dimension>
void
NodeListRegistrar<Dimension>::
domainDecompositionIndependent(const bool x) {
  mDomainDecompIndependent = x;
}

}
