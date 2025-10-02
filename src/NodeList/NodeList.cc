//---------------------------------Spheral++----------------------------------//
// NodeList -- An abstract base class for the NodeLists.
//
// We will define here the basic functionality we expect all NodeLists to 
// provide.
//
// Created by JMO, Wed Sep  8 21:54:50 PDT 1999
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "NodeListRegistrar.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Neighbor/Neighbor.hh"
#include "DataBase/State.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/DBC.hh"
#include "Utilities/packElement.hh"
#include "Utilities/SpheralMessage.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "NodeList.hh"

#include <algorithm>
using std::vector;
using std::list;
using std::string;
using std::cerr;
using std::endl;

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor with optional numInternal nodes, numGhost nodes, and name.
//------------------------------------------------------------------------------
template<typename Dimension>
NodeList<Dimension>::NodeList(std::string name,
                              const size_t numInternal,
                              const size_t numGhost,
                              const Scalar hmin,
                              const Scalar hmax,
                              const Scalar hminratio,
                              const Scalar nPerh,
                              const size_t maxNumNeighbors):
  mNumNodes(numInternal + numGhost),
  mFirstGhostNode(numInternal),
  mName(name),
  mMass(HydroFieldNames::mass),
  mPositions(HydroFieldNames::position),
  mVelocity(HydroFieldNames::velocity),
  mH(HydroFieldNames::H),
  mWork(HydroFieldNames::work),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mNodesPerSmoothingScale(nPerh),
  mMaxNumNeighbors(maxNumNeighbors),
  mFieldBases(),
  mNeighborPtr(nullptr),
  mDummyList(),
  mRestart(registerWithRestart(*this, 10)) {
  NodeListRegistrar<Dimension>::instance().registerNodeList(*this);
  mMass.setNodeList(*this);
  mPositions.setNodeList(*this);
  mVelocity.setNodeList(*this);
  mH.setNodeList(*this);
  mWork.setNodeList(*this);
  mDummyList.push_back(this);
  // It's never valid to have zero H's.
  mH = SymTensor::one;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NodeList<Dimension>::~NodeList() {
  DEBUG_LOG << "NodeList::~NodeList " << mName << " " << this;
  auto startingFields = mFieldBases;
  for (auto x: startingFields) x.get().unregisterNodeList();

  // Unregister ourselves from the NodeListRegistrar, freeing up our name.
  NodeListRegistrar<Dimension>::instance().unregisterNodeList(*this);

  // After we're done, all the field should have unregistered themselves
  // from the Node List.
  ENSURE(numFields() == 0u);
}

//------------------------------------------------------------------------------
// Set the number of internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::numInternalNodes(size_t size) {
  const auto numGhost = numGhostNodes();
  const auto oldFirstGhostNode = mFirstGhostNode;
  mFirstGhostNode = size;
  mNumNodes = size + numGhost;
  for (auto x: mFieldBases)     x.get().resizeFieldInternal(size, oldFirstGhostNode);
}

//------------------------------------------------------------------------------
// Set the number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::numGhostNodes(size_t size) {
  const auto numInternal = numInternalNodes();
  mNumNodes = numInternal + size;
  for (auto x: mFieldBases)     x.get().resizeFieldGhost(size);
}

//------------------------------------------------------------------------------
// Provide NodeIterators for all nodes in this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
AllNodeIterator<Dimension>
NodeList<Dimension>::nodeBegin() const {
  return AllNodeIterator<Dimension>(mDummyList.begin(),
                                    mDummyList.begin(),
                                    mDummyList.end());
}

template<typename Dimension>
AllNodeIterator<Dimension>
NodeList<Dimension>::nodeEnd() const {
  return AllNodeIterator<Dimension>(mDummyList.end(),
                                    mDummyList.begin(),
                                    mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
InternalNodeIterator<Dimension>
NodeList<Dimension>::internalNodeBegin() const {
  return InternalNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end());
}

template<typename Dimension>
InternalNodeIterator<Dimension>
NodeList<Dimension>::internalNodeEnd() const {
  return InternalNodeIterator<Dimension>(mDummyList.end(),
                                         mDummyList.begin(),
                                         mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
GhostNodeIterator<Dimension>
NodeList<Dimension>::ghostNodeBegin() const {
  return GhostNodeIterator<Dimension>(mDummyList.begin(),
                                      mDummyList.begin(),
                                      mDummyList.end(),
                                      firstGhostNode());
}

template<typename Dimension>
GhostNodeIterator<Dimension>
NodeList<Dimension>::ghostNodeEnd() const {
  return GhostNodeIterator<Dimension>(mDummyList.end(),
                                      mDummyList.begin(),
                                      mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for master neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>
NodeList<Dimension>::masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const {
  REQUIRE(mNeighborPtr != 0);
  REQUIRE(masterLists.size() == 1);
  if (masterLists[0].size() > 0) {
    return MasterNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         masterLists[0].begin(),
                                         masterLists);
  } else {
    return this->masterNodeEnd();
  }
}

template<typename Dimension>
MasterNodeIterator<Dimension>
NodeList<Dimension>::masterNodeEnd() const {
  return MasterNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node iterators for coarse neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>
NodeList<Dimension>::coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const {
  REQUIRE(mNeighborPtr != 0);
  REQUIRE(coarseNeighbors.size() == 1);
  if (coarseNeighbors[0].size() > 0) {
    return CoarseNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         coarseNeighbors[0].begin(),
                                         coarseNeighbors);
  } else {
    return this->coarseNodeEnd();
  }
}

template<typename Dimension>
CoarseNodeIterator<Dimension>
NodeList<Dimension>::coarseNodeEnd() const {
  return CoarseNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node iterators for refine neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>
NodeList<Dimension>::refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const {
  REQUIRE(mNeighborPtr != 0);
  REQUIRE(refineNeighbors.size() == 1);
  if (refineNeighbors[0].size() > 0) {
    return RefineNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         refineNeighbors[0].begin(),
                                         refineNeighbors);
  } else {
    return this->refineNodeEnd();
  }
}

template<typename Dimension>
RefineNodeIterator<Dimension>
NodeList<Dimension>::refineNodeEnd() const {
  return RefineNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Set the mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
mass(const Field<Dimension, typename Dimension::Scalar>& m) {
  mMass = m;
  mMass.name(HydroFieldNames::mass);
}

//------------------------------------------------------------------------------
// Set the position field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
positions(const Field<Dimension, typename Dimension::Vector>& r) {
  mPositions = r;
  mPositions.name(HydroFieldNames::position);
}

//------------------------------------------------------------------------------
// Set the velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
velocity(const Field<Dimension, typename Dimension::Vector>& v) {
  mVelocity = v;
  mVelocity.name(HydroFieldNames::velocity);
}

//------------------------------------------------------------------------------
// Set the H field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
Hfield(const Field<Dimension, typename Dimension::SymTensor>& H) {
  mH = H;
  mH.name(HydroFieldNames::H);
}

//------------------------------------------------------------------------------
// Set the work field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
work(const Field<Dimension, typename Dimension::Scalar>& w) {
  mWork = w;
  mWork.name(HydroFieldNames::work);
}

//------------------------------------------------------------------------------
// Compute and return the inverse field of H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
Hinverse(Field<Dimension, typename Dimension::SymTensor>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  for (auto i = 0u; i < numInternalNodes(); ++i) field(i) = mH(i).Inverse();
  field.name("H inverse");
}

//------------------------------------------------------------------------------
// Check if the given field is registered with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NodeList<Dimension>::haveField(const FieldBase<Dimension>& field) const {
  return (find_if(mFieldBases.begin(),
                  mFieldBases.end(),
                  [&](const std::reference_wrapper<FieldBase<Dimension>>& x) { return &(x.get()) == &field; }) != mFieldBases.end());
}

//------------------------------------------------------------------------------
// Return the type of node (internal or ghost).
//------------------------------------------------------------------------------
template<typename Dimension>
NodeType
NodeList<Dimension>::nodeType(size_t i) const {
  REQUIRE(i < mNumNodes);
  return (i < mFirstGhostNode ? 
          NodeType::InternalNode :
          NodeType::GhostNode);
}

//------------------------------------------------------------------------------
// Delete the given node indicies from the NodeList (including all Fields).
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
deleteNodes(const vector<size_t>& nodeIDs) {
  if (nodeIDs.size() > 0) {

    // First sort and make sure all node IDs are valid.
    vector<size_t> uniqueIDs(nodeIDs);
    std::sort(uniqueIDs.begin(), uniqueIDs.end());
    uniqueIDs.erase(std::unique(uniqueIDs.begin(), uniqueIDs.end()), uniqueIDs.end());
    CHECK(uniqueIDs.size() <= numNodes());
    CHECK(uniqueIDs.size() == 0 or uniqueIDs.back() < this->numNodes());

    // Determine how many internal, ghost, and total nodes we should end with.
    auto ghostDeleteItr = uniqueIDs.begin();
    while (ghostDeleteItr < uniqueIDs.end() &&
           *ghostDeleteItr < mFirstGhostNode) ++ghostDeleteItr;
    CHECK(ghostDeleteItr >= uniqueIDs.begin() && ghostDeleteItr <= uniqueIDs.end());
    const size_t numInternalNodesRemoved = std::distance(uniqueIDs.begin(), ghostDeleteItr);
    CHECK(numInternalNodesRemoved <= numInternalNodes());
    mNumNodes -= uniqueIDs.size();
    mFirstGhostNode -= numInternalNodesRemoved;
    CHECK(mFirstGhostNode <= mNumNodes);

    // Now iterate over the Fields defined on this NodeList, and remove the appropriate
    // elements from each.
    for (auto x: mFieldBases)     x.get().deleteElements(uniqueIDs);
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  for (auto x: mFieldBases) {
    CONTRACT_VAR(x) ;
    ENSURE(x.get().size() == mNumNodes);
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Pack the Field values for the given node indicies as appends onto the given
// vector of Scalars.
//------------------------------------------------------------------------------
template<typename Dimension>
list< vector<char> >
NodeList<Dimension>::
packNodeFieldValues(const vector<size_t>& nodeIDs) const {

  // Prepare the result.
  list<vector<char>> result;

  // Sort and make sure all node IDs are valid.
  vector<size_t> uniqueIDs(nodeIDs);
  std::sort(uniqueIDs.begin(), uniqueIDs.end());
  uniqueIDs.erase(unique(uniqueIDs.begin(), uniqueIDs.end()), uniqueIDs.end());
  CHECK(uniqueIDs.size() <= numNodes());
  CHECK(uniqueIDs.size() == 0 or uniqueIDs.back() < this->numNodes());

  // Iterate over all the Fields defined on this NodeList, and append it's packed 
  // field values to the stack.
  for (auto x: mFieldBases) {
    result.push_back(x.get().packValues(uniqueIDs));
  }

  ENSURE(result.size() == mFieldBases.size());
  return result;
}

//------------------------------------------------------------------------------
// Append the requested number of nodes onto this Field, and fill in the
// Field values using the provided packed data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
appendInternalNodes(const size_t numNewNodes,
                    const list<vector<char>>& packedFieldValues) {

  REQUIRE(packedFieldValues.size() == mFieldBases.size());

  // We only work if there are new nodes.
  if (numNewNodes > 0) {

    // Begin by resizing this NodeList appropriately.
    const auto beginInsertionIndex = numInternalNodes();
    numInternalNodes(beginInsertionIndex + numNewNodes);
    CHECK(numInternalNodes() == beginInsertionIndex + numNewNodes);

    // Loop over each Field, and have them fill in the new values from the
    // packed char buffers.
    vector<size_t> nodeIDs(numNewNodes);
    for (auto i = 0u; i < numNewNodes; ++i) nodeIDs[i] = beginInsertionIndex + i;
    auto bufItr = packedFieldValues.begin();
    for (auto x: mFieldBases) {
      CHECK(bufItr != packedFieldValues.end());
      x.get().unpackValues(nodeIDs, *bufItr);
      ++bufItr;
    }
    CHECK(bufItr == packedFieldValues.end());
  }
}

//------------------------------------------------------------------------------
// Reorder the nodes according to the requested ordering.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
reorderNodes(const vector<size_t>& newOrdering) {

  // The number of internal nodes.
  const auto n = this->numInternalNodes();

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(newOrdering.size() == n);
    vector<size_t> tmp(newOrdering);
    std::sort(tmp.begin(), tmp.end());
    for (auto i = 0u; i < n; ++i) REQUIRE(tmp[i] == i);
  }
  END_CONTRACT_SCOPE

  // Make sure we're not carting around ghost nodes.
  this->numGhostNodes(0);

  // The original ordering.
  vector<size_t> oldOrdering(n);
  for (auto i = 0u; i < n; ++i) oldOrdering[i] = i;

  // Pack up all the current nodal field values.
  list<vector<char>> packedFieldValues;
  for (auto x: mFieldBases) packedFieldValues.push_back(x.get().packValues(oldOrdering));
  CHECK(packedFieldValues.size() == mFieldBases.size());

  // Now unpack in the desired order.
  auto bufItr = packedFieldValues.begin();
  for (auto x: mFieldBases) {
    x.get().unpackValues(newOrdering, *bufItr);
    ++bufItr;
  }
  CHECK(bufItr == packedFieldValues.end());

  // Post-conditions.
  ENSURE(this->numInternalNodes() == n);
}

//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Dump the name of the NodeList.
  file.write(name(), pathName + "/name");

  // Dump the number of internal nodes (we assume that ghost information
  // does not need to be stored).
  file.write(numInternalNodes(), pathName + "/numNodes");

  // Dump each of the internal fields of the NodeList.
  file.write(mMass, pathName + "/mass");
  file.write(mPositions, pathName + "/positions");
  file.write(mVelocity, pathName + "/velocity");
  file.write(mH, pathName + "/H");
  file.write(mWork, pathName + "/work");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Restore the name of the NodeList.
  file.read(mName, pathName + "/name");

  // Read and reset the number of internal nodes.
  size_t numNodes;
  file.read(numNodes, pathName + "/numNodes");
  numInternalNodes(numNodes);

  // Now we can restore each of the internal fields of the NodeList.
  file.read(mMass, pathName + "/mass");
  file.read(mPositions, pathName + "/positions");
  file.read(mVelocity, pathName + "/velocity");
  file.read(mH, pathName + "/H");
  file.read(mWork, pathName + "/work");

  // The neighbor object doesn't actually write out state, but does need to be
  // reinitialized with the new NodeList state.
  mNeighborPtr->updateNodes();
}

//------------------------------------------------------------------------------
// Register a field with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::registerField(FieldBase<Dimension>& field) const {
  DEBUG_LOG << "NodeList::registerField : " << mName << " " << this << " : " << field.name() << " " << &field;
  if (haveField(field)) {
    SpheralMessage("WARNING: Attempt to register field " << &field << " (" << field.name()
                   << ") with NodeList " << this << " (" << this->name() << ") that already has it.");
  } else {
    mFieldBases.push_back(std::ref(field));
  }
}

//------------------------------------------------------------------------------
// Unregister a field that is listed with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::unregisterField(FieldBase<Dimension>& field) const {
  DEBUG_LOG << "NodeList::unregisterField : " << mName << " " << this << " : " << field.name() << " " << &field;
#pragma omp critical
  {
    auto itr = find_if(mFieldBases.begin(),
                       mFieldBases.end(),
                       [&](const std::reference_wrapper<FieldBase<Dimension>>& x) { return &(x.get()) == &field; });
    if (itr != mFieldBases.end()) mFieldBases.erase(itr);
  }
}

//------------------------------------------------------------------------------
// Register the given neighbor object with this node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
registerNeighbor(Neighbor<Dimension>& neighbor) {
  DEBUG_LOG << "NodeList::registerNeighbor : " << mName << " " << this << " : " << &neighbor;
  mNeighborPtr = &neighbor;
}

//------------------------------------------------------------------------------
// Unregister the current neighbor object from this node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::unregisterNeighbor() {
  DEBUG_LOG << "NodeList::unregisterNeighbor : " << mName << " " << this;
  mNeighborPtr = nullptr;
}

// //------------------------------------------------------------------------------
// // Notify all Fields registered on this NodeList to cache their coarse neighbor
// // values.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// NodeList<Dimension>::
// notifyFieldsCacheCoarseValues() const {
//   for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBases.begin();
//        fieldItr < mFieldBases.end();
//        ++fieldItr) {
//     (*fieldItr)->notifyNewCoarseNodes();
//   }
// }

// //------------------------------------------------------------------------------
// // Notify all Fields registered on this NodeList to cache their refine neighbor
// // values.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// NodeList<Dimension>::
// notifyFieldsCacheRefineValues() const {
//   for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBases.begin();
//        fieldItr < mFieldBases.end();
//        ++fieldItr) {
//     (*fieldItr)->notifyNewRefineNodes();
//   }
// }
}
