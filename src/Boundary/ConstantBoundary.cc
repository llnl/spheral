//---------------------------------Spheral++----------------------------------//
// ConstantBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Boundary.hh"
#include "Boundary/mapPositionThroughPlanes.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/DBC.hh"

#include "Boundary/ConstantBoundary.hh"
#include "Boundary/ConstantBoundaryUtilities.hh"

using std::vector;
using std::map;
using std::string;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantBoundary<Dimension>::
ConstantBoundary(DataBase<Dimension>& dataBase,
                 NodeList<Dimension>& nodeList,
                 const vector<size_t>& nodeIDs,
                 const GeomPlane<Dimension>& denialPlane):
  Boundary<Dimension>(),
  mDataBase(dataBase),
  mNodeListPtr(&nodeList),
  mBoundaryCount(nodeList.numFields()),
  mNodeFlags("ConstantBoundaryNodeFlags" + std::to_string(mBoundaryCount), nodeList, 0u),
  mNumConstantNodes(nodeIDs.size()),
  mDenialPlane(denialPlane),
  mReflectOperator(planarReflectingOperator(denialPlane)),
  mActive(false),
  mBufferedValues(),
  mRestart(registerWithRestart(*this)),
  mRedistribution(registerWithRedistribution(*this,
                                             &ConstantBoundary<Dimension>::notifyBeforeRedistribution,
                                             &ConstantBoundary<Dimension>::notifyAfterRedistribution)) {

  // Store the ids of the nodes we're watching.
  for (auto i: nodeIDs) {
    REQUIRE(i < nodeList.numInternalNodes());
    mNodeFlags[i] = 1u;
  }

  // Issue a big old warning!
  // if (Process::getRank() == 0) cerr << "WARNING: ConstantBoundary is currently not compatible with redistributing nodes!\nMake sure you don't allow redistribution with this Boundary condition." << endl;
}

//------------------------------------------------------------------------------
// setGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);

  if (mActive and &nodeList == mNodeListPtr) {
    auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    auto& cNodes = boundaryNodes.controlNodes;
    auto& gNodes = boundaryNodes.ghostNodes;
    
    auto currentNumGhostNodes = nodeList.numGhostNodes();
    auto firstNewGhostNode = nodeList.numNodes();
    nodeList.numGhostNodes(currentNumGhostNodes + mNumConstantNodes);
    cNodes.resize(mNumConstantNodes);
    gNodes.resize(mNumConstantNodes);
    for (auto i = 0u; i < mNumConstantNodes; ++i) {
      const auto j = firstNewGhostNode + i;
      mNodeFlags(j) = 1u;
      cNodes[i] = j;
      gNodes[i] = j;
    }
    
    this->updateGhostNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// updateGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  if (mActive and &nodeList == mNodeListPtr) {
    this->applyGhostBoundary(nodeList.positions());
    this->applyGhostBoundary(nodeList.Hfield());
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to the ghost nodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(FieldBase<Dimension>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mBufferedValues, false);
}

//------------------------------------------------------------------------------
// setViolationNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
  if (&nodeList == mNodeListPtr) {
    const auto& pos = nodeList.positions();
    auto&       boundaryNodes = this->accessBoundaryNodes(nodeList);
    auto&       vNodes = boundaryNodes.violationNodes;
    for (auto i = 0u; i < nodeList.numInternalNodes(); ++i) {
      if (mDenialPlane.compare(pos(i)) == -1) vNodes.push_back(i);
      // if (mNodeFlags[i] == 1 or mDenialPlane.compare(pos(i)) == -1) vNodes.push_back(i);
    }
    this->updateViolationNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// updateViolationNodes, correcting positions and H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
  if (&nodeList == mNodeListPtr) {
    auto& pos = nodeList.positions();
    auto& vel = nodeList.velocity();
    auto& H = nodeList.Hfield();

    // Look for any nodes in violation of the plane and reset their positions
    // and H's.
    for (auto i = 0u; i < nodeList.numInternalNodes(); ++i) {
      if (mDenialPlane.compare(pos(i)) == 1) {
        pos(i) = mapPositionThroughPlanes(pos(i), mDenialPlane, mDenialPlane);
        vel(i) = mReflectOperator*vel(i);
        H(i) = (mReflectOperator*H(i)*mReflectOperator).Symmetric();
      }
    }
  }
}

//------------------------------------------------------------------------------
// enforceBoundary
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(FieldBase<Dimension>&) const {
  // resetValues(field, this->nodeIndices(), mIntValues);
}

//------------------------------------------------------------------------------
// On problem startup we take our snapshot of the state.  We also then
// destroy the original internal nodes, as we will be replacing them
// with ghosts from now on.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::initializeProblemStartup(const bool final) {

  // Clear any existing data.
  mBufferedValues.clear();

  // Now take a snapshot of the Fields.
  const auto nodeIDs = this->nodeIndices();
  // cerr << "Node IDs: ";
  // std::copy(nodeIDs.begin(), nodeIDs.end(), std::ostream_iterator<int>(std::cerr, " "));
  // cerr << endl;
  storeFieldValues(*mNodeListPtr, nodeIDs, mBufferedValues);

  // If we're in cylindrical symmetry (RZ) we need to convert the mass to mass/(2*pi*r).
  if (GeometryRegistrar::coords() == CoordinateType::RZ){
    const auto& pos = mNodeListPtr->positions();
    const auto& mass = mNodeListPtr->mass();
    std::vector<char> buffer;
    for (auto i: nodeIDs) {
      const auto circi = 2.0*M_PI*abs(pos(i).y());
      CHECK(circi > 0.0);
      const auto mi = mass(i)/circi;
      packElement(mi, buffer);
    }
    const auto mkey = StateBase<Dimension>::key(mass);
    CHECK(mBufferedValues.find(mkey) != mBufferedValues.end());
    mBufferedValues[mkey] = buffer;
  }

  // Remove the origial internal nodes.
  if (final) {
    VERIFY2(not mActive, "ConstantBoundary::initializeProblemStartup ERROR -- called with final=True more than once");
    mNodeListPtr->deleteNodes(nodeIDs);
    mActive = true;
  }
}

//------------------------------------------------------------------------------
// Return the set of node IDs we're controlling.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<size_t>
ConstantBoundary<Dimension>::
nodeIndices() const {
  std::vector<size_t> result;
  if (mActive) {
    for (const auto& bnitr: this->accessBoundaryNodes()) {
      const auto& ghosts = bnitr.second.ghostNodes;
      std::copy(ghosts.begin(), ghosts.end(), std::back_inserter(result));
    }
  } else {
    for (auto i = 0u; i != mNodeListPtr->numNodes(); ++i) {
      if (mNodeFlags(i) == 1u) result.push_back(i);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Return a unique label for restart.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
ConstantBoundary<Dimension>::label() const {
  return "ConstantBoundary" + std::to_string(mBoundaryCount);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mNumConstantNodes, pathName + "/numConstantNodes");
  file.write(mActive, pathName + "/active");
  file.write(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  for (const auto& p: mBufferedValues) {
    keys.push_back(p.first);
    std::string val(p.second.begin(), p.second.end());
    file.write(val, pathName + "/BufferedValues/" + p.first);
  }
  file.write(keys, pathName + "/keys");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName)  {
  file.read(mNumConstantNodes, pathName + "/numConstantNodes");
  file.read(mActive, pathName + "/active");
  file.read(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  file.read(keys, pathName + "/keys");
  mBufferedValues.clear();
  for (const auto& key: keys) {
    std::string val;
    file.read(val, pathName + "/BufferedValues/" + key);
    mBufferedValues[key] = vector<char>(val.begin(), val.end());
  }
}

//------------------------------------------------------------------------------
// Redistribution methods
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
notifyBeforeRedistribution() {
  VERIFY2(false, "ConstantBoundary ERROR: node redistribution not allowed with constant boundaries.");
}

template<typename Dimension>
void
ConstantBoundary<Dimension>::
notifyAfterRedistribution() {
  VERIFY2(false, "ConstantBoundary ERROR: node redistribution not allowed with constant boundaries.");
}

}
