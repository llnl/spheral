//---------------------------------Spheral++----------------------------------//
// NodeListBase -- host-side base class for NodeLists.
//
// Provides identity, registration, and common host configuration for NodeLists.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeListBase_hh__
#define __Spheral_NodeListBase_hh__

#include "Utilities/span.hh"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace Spheral {

template<typename Dimension> class FieldBase;
template<typename Dimension> class Neighbor;
class FileIO;

enum class NodeType {
  InternalNode = 0,
  GhostNode = 1
};

template<typename Dimension>
class NodeListBase {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;

  using FieldBaseSpan = SPHERAL_SPAN_TYPE<std::reference_wrapper<FieldBase<Dimension>>>;

  using FieldBaseIterator = typename std::vector<std::reference_wrapper<FieldBase<Dimension>>>::iterator;
  using const_FieldBaseIterator = typename std::vector<std::reference_wrapper<FieldBase<Dimension>>>::const_iterator;

  // Constructors and destructor.
  NodeListBase(std::string name,
               const Scalar hmin,
               const Scalar hmax,
               const Scalar hminratio,
               const Scalar nPerh,
               const size_t maxNumNeighbors);
  virtual ~NodeListBase() = default;

  // Access the name of the NodeList.
  std::string name() const;

  // Node count contract.
  virtual size_t numNodes() const = 0;
  virtual size_t numInternalNodes() const = 0;
  virtual size_t numGhostNodes() const = 0;
  virtual size_t firstGhostNode() const = 0;
  virtual void numInternalNodes(size_t size) = 0;
  virtual void numGhostNodes(size_t size) = 0;

  // Provide iterators over the FieldBases defined on this NodeList.
  FieldBaseIterator registeredFieldsBegin();
  FieldBaseIterator registeredFieldsEnd();

  const_FieldBaseIterator registeredFieldsBegin() const;
  const_FieldBaseIterator registeredFieldsEnd() const;

  FieldBaseSpan registeredFields() const;

  // Provide methods to query Fields defined over this NodeList.
  size_t numFields() const;
  bool haveField(const FieldBase<Dimension>& field) const;

  // Access the neighbor object.
  Neighbor<Dimension>& neighbor() const;

  // Neighbor object registration.
  void registerNeighbor(Neighbor<Dimension>& neighbor);
  void unregisterNeighbor();

  // The target number of nodes per smoothing scale (for calculating the ideal H).
  Scalar nodesPerSmoothingScale() const;
  void nodesPerSmoothingScale(Scalar val);

  // The maximum number of neighbors we want to have (for calculating the ideal H).
  size_t maxNumNeighbors() const;
  void maxNumNeighbors(size_t val);

  // Allowed range of smoothing scales for use in calculating H.
  Scalar hmin() const;
  void hmin(Scalar val);

  Scalar hmax() const;
  void hmax(Scalar val);

  Scalar hminratio() const;
  void hminratio(Scalar val);

  // Methods required for restarting.
  virtual std::string label() const = 0;
  virtual void dumpState(FileIO& file, const std::string& pathName) const = 0;
  virtual void restoreState(const FileIO& file, const std::string& pathName) = 0;

  // Some operators.
  bool operator==(const NodeListBase& rhs) const;
  bool operator!=(const NodeListBase& rhs) const;

  // No default constructor, copying, or assignment.
  NodeListBase() = delete;
  NodeListBase(const NodeListBase& nodes) = delete;
  NodeListBase& operator=(const NodeListBase& rhs) = delete;

protected:
  //--------------------------- Protected Interface ---------------------------//
  void name(std::string name);

  // Methods to handle registering Fields and Neighbors.
  void registerField(FieldBase<Dimension>& field) const;
  void unregisterField(FieldBase<Dimension>& field) const;

  friend class FieldBase<Dimension>;

  std::string mName;

  // List of Fields that are defined over this NodeList.
  mutable std::vector<std::reference_wrapper<FieldBase<Dimension>>> mFieldBases;

  Neighbor<Dimension>* mNeighborPtr;

  // Stuff for how H is handled.
  Scalar mhmin, mhmax, mhminratio;
  Scalar mNodesPerSmoothingScale;
  size_t mMaxNumNeighbors;
};

}

#include "NodeListBaseInline.hh"

#endif
