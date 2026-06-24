//---------------------------------Spheral++----------------------------------//
// NodeList -- An abstract base class for the NodeLists.
//
// We will define here the basic functionality we expect all NodeLists to 
// provide.
//
// Created by JMO, Tue Jun  8 23:58:36 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeList__
#define __Spheral_NodeList__

#include "DataOutput/registerWithRestart.hh"
#include "NodeList/NodeListBase.hh"
#include "NodeList/NodeListView.hh"
#include "Utilities/DBC.hh"

#include <string>
#include <list>
#include <vector>

namespace Spheral {

template<typename Dimension> class AllNodeIterator;
template<typename Dimension> class InternalNodeIterator;
template<typename Dimension> class GhostNodeIterator;
template<typename Dimension> class MasterNodeIterator;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class State;
template<typename Dimension> class Neighbor;
template<typename Dimension> class DataBase;
template<typename Dimension> class TableKernel;
template<typename Dimension> class FieldBase;
template<typename Dimension, typename Value> class Field;
class FileIO;

template<typename Dimension>
class NodeList:
    public NodeListBase<Dimension>,
    public NodeListView<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using BaseType = NodeListBase<Dimension>;
  using ViewType = NodeListView<Dimension>;

  using FieldBaseSpan = typename BaseType::FieldBaseSpan;
  using FieldBaseIterator = typename BaseType::FieldBaseIterator;
  using const_FieldBaseIterator = typename BaseType::const_FieldBaseIterator;

  using BaseType::registeredFieldsBegin;
  using BaseType::registeredFieldsEnd;
  using BaseType::registeredFields;
  using BaseType::numFields;
  using BaseType::haveField;
  using BaseType::neighbor;
  using BaseType::registerNeighbor;
  using BaseType::unregisterNeighbor;
  using BaseType::nodesPerSmoothingScale;
  using BaseType::maxNumNeighbors;
  using BaseType::hmin;
  using BaseType::hmax;
  using BaseType::hminratio;

  // Constructors
  explicit NodeList(std::string name,
                    const size_t numInternal,
                    const size_t numGhost,
                    const Scalar hmin = 1.0e-20,
                    const Scalar hmax = 1.0e20,
                    const Scalar hminratio = 0.1,
                    const Scalar nPerh = 2.01,
                    const size_t maxNumNeighbors = 500);

  // Destructor
  virtual ~NodeList();

  // Access the name of the NodeList.
  std::string name() const { return BaseType::name(); }

  // Get or set the number of Nodes.
  virtual size_t numNodes()                           const override { return ViewType::numNodes(); }
  virtual size_t numInternalNodes()                   const override { return ViewType::numInternalNodes(); }
  virtual size_t numGhostNodes()                      const override { CHECK(ViewType::firstGhostNode() <= ViewType::numNodes()); return ViewType::numNodes() - ViewType::firstGhostNode(); }
  virtual size_t firstGhostNode()                     const override { CHECK(ViewType::firstGhostNode() <= ViewType::numNodes()); return ViewType::firstGhostNode(); }
  virtual void numInternalNodes(size_t size) override;
  virtual void numGhostNodes(size_t size) override;

  // Provide the standard NodeIterators over the nodes of this NodeList.
  AllNodeIterator<Dimension> nodeBegin() const;
  AllNodeIterator<Dimension> nodeEnd() const;
          
  InternalNodeIterator<Dimension> internalNodeBegin() const;
  InternalNodeIterator<Dimension> internalNodeEnd() const;
          
  GhostNodeIterator<Dimension> ghostNodeBegin() const;
  GhostNodeIterator<Dimension> ghostNodeEnd() const;
          
  MasterNodeIterator<Dimension> masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const;
  MasterNodeIterator<Dimension> masterNodeEnd() const;
          
  CoarseNodeIterator<Dimension> coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const;
  CoarseNodeIterator<Dimension> coarseNodeEnd() const;

  RefineNodeIterator<Dimension> refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const;
  RefineNodeIterator<Dimension> refineNodeEnd() const;

  // The NodeList state Fields.
  Field<Dimension, Scalar>& mass()                          { return mMass; }
  Field<Dimension, Vector>& positions()                     { return mPositions; }
  Field<Dimension, Vector>& velocity()                      { return mVelocity; } 
  Field<Dimension, SymTensor>& Hfield()                     { return mH; }

  const Field<Dimension, Scalar>& mass()              const { return mMass; }
  const Field<Dimension, Vector>& positions()         const { return mPositions; }
  const Field<Dimension, Vector>& velocity()          const { return mVelocity; } 
  const Field<Dimension, SymTensor>& Hfield()         const { return mH; }

  void mass(const Field<Dimension, Scalar>& m);
  void positions(const Field<Dimension, Vector>& r);
  void velocity(const Field<Dimension, Vector>& v);
  void Hfield(const Field<Dimension, SymTensor>& H);

  // Anyone can acces the work field as a mutable field.
  Field<Dimension, Scalar>& work()                    const { return mWork; }
  void work(const Field<Dimension, Scalar>& w);

  // These are quantities which are not stored, but can be computed.
  void Hinverse(Field<Dimension, SymTensor>& field) const;

  // NodeLists can contain ghost nodes (either communicated from neighbor
  // processors, or simply created for boundary conditions).
  NodeType nodeType(size_t i) const;

  //****************************************************************************
  // Methods for adding/removing individual nodes to/from the NodeList
  // (including all Field information.  These methods are primarily useful
  // for redistributing Nodes between parallel domains.
  virtual void deleteNodes(const std::vector<size_t>& nodeIDs);
  virtual std::list<std::vector<char>>  packNodeFieldValues(const std::vector<size_t>& nodeIDs) const;
  virtual void appendInternalNodes(const size_t numNewNodes,
                                   const std::list<std::vector<char>>& packedFieldValues);

  // A related method for reordering the nodes.
  virtual void reorderNodes(const std::vector<size_t>& newOrdering);
  //****************************************************************************

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label()                              const override { return "NodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

  // Some operators.
  bool operator==(const NodeList& rhs)                     const { return this == &rhs; }
  bool operator!=(const NodeList& rhs)                     const { return !(*this == rhs); }

  // Get the view of NodeList-owned data needed by ConnectivityMap.
  ViewType view() const;

  // No default constructor, copying, or assignment.
  NodeList() = delete;
  NodeList(const NodeList& nodes) = delete;
  NodeList& operator=(const NodeList& rhs) = delete;

protected:
  //--------------------------- Protected Interface ---------------------------//
  void refreshView();

  void registerField(FieldBase<Dimension>& field) const { BaseType::registerField(field); }
  void unregisterField(FieldBase<Dimension>& field) const { BaseType::unregisterField(field); }

  friend class FieldBase<Dimension>;

private:
  //--------------------------- Private Interface ---------------------------//
  // State fields.
  Field<Dimension, Scalar> mMass;
  Field<Dimension, Vector> mPositions;
  Field<Dimension, Vector> mVelocity;
  Field<Dimension, SymTensor> mH;

  // The work field is mutable.
  mutable Field<Dimension, Scalar> mWork;

  // Provide a dummy vector to buid NodeIterators against.
  std::vector<NodeList<Dimension>*> mDummyList;

  // The restart registration.
  RestartRegistrationType mRestart;

  using BaseType::mName;
  using BaseType::mFieldBases;
  using BaseType::mNeighborPtr;
  using ViewType::mNumNodes;
  using ViewType::mFirstGhostNode;
};

}

#include "Field/Field.hh"

#endif
