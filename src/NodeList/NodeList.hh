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
#include "NodeList/NodeListView.hh"
#include "Utilities/span.hh"

#include <string>
#include <list>
#include <vector>
#include <functional>
#include <memory>

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
class NodeList: public NodeListView<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ViewType = NodeListView<Dimension>;
  using ViewType::numInternalNodes;
  using ViewType::numGhostNodes;

  using FieldBaseSpan = SPHERAL_SPAN_TYPE<std::reference_wrapper<FieldBase<Dimension>>>;

  using FieldBaseIterator = typename std::vector<std::reference_wrapper<FieldBase<Dimension>>>::iterator;
  using const_FieldBaseIterator = typename std::vector<std::reference_wrapper<FieldBase<Dimension>>>::const_iterator;

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
  std::string name()                                  const { return mName; }

  // Set the number of Nodes.
  void numInternalNodes(size_t size);
  void numGhostNodes(size_t size);

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
  Field<Dimension, Scalar>& massRZ()                        { VERIFY2(mMassRZptr, "NodeList::massRZ ERROR: not allocated"); return *mMassRZptr; }
  Field<Dimension, Vector>& positions()                     { return mPositions; }
  Field<Dimension, Vector>& velocity()                      { return mVelocity; } 
  Field<Dimension, SymTensor>& Hfield()                     { return mH; }

  const Field<Dimension, Scalar>& mass()              const { return mMass; }
  const Field<Dimension, Scalar>& massRZ()            const { VERIFY2(mMassRZptr, "NodeList::massRZ ERROR: not allocated"); return *mMassRZptr; }
  const Field<Dimension, Vector>& positions()         const { return mPositions; }
  const Field<Dimension, Vector>& velocity()          const { return mVelocity; } 
  const Field<Dimension, SymTensor>& Hfield()         const { return mH; }

  void mass(const Field<Dimension, Scalar>& m);
  void massRZ(const Field<Dimension, Scalar>& m);
  void positions(const Field<Dimension, Vector>& r);
  void velocity(const Field<Dimension, Vector>& v);
  void Hfield(const Field<Dimension, SymTensor>& H);

  // Anyone can acces the work field as a mutable field.
  Field<Dimension, Scalar>& work()                    const { return mWork; }
  void work(const Field<Dimension, Scalar>& w);

  // These are quantities which are not stored, but can be computed.
  void Hinverse(Field<Dimension, SymTensor>& field) const;

  // Provide iterators over the FieldBases defined on this NodeList.
  FieldBaseIterator registeredFieldsBegin()                     { return mFieldBases.begin(); }    
  FieldBaseIterator registeredFieldsEnd()                       { return mFieldBases.end(); }      

  const_FieldBaseIterator registeredFieldsBegin()         const { return mFieldBases.begin(); }    
  const_FieldBaseIterator registeredFieldsEnd()           const { return mFieldBases.end(); }       

  FieldBaseSpan registeredFields()                        const { return FieldBaseSpan(mFieldBases.data(), mFieldBases.size()); }

  // Provide methods to add and subtract Fields which are defined over a
  // NodeList.
  size_t numFields() const { return mFieldBases.size(); }
  bool haveField(const FieldBase<Dimension>& field) const;

  // Access the neighbor object.
  Neighbor<Dimension>& neighbor()                          const { CHECK(mNeighborPtr != nullptr); return *mNeighborPtr; }

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
  virtual std::string label()                              const { return "NodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Some operators.
  bool operator==(const NodeList& rhs)                     const { return this == &rhs; }
  bool operator!=(const NodeList& rhs)                     const { return !(*this == rhs); }

  // Get the view of NodeList-owned data.
  ViewType view();

  // Neighbor object registration
  void registerNeighbor(Neighbor<Dimension>& neighbor);
  void unregisterNeighbor();

  // No default constructor, copying, or assignment.
  NodeList() = delete;
  NodeList(const NodeList& nodes) = delete;
  NodeList& operator=(const NodeList& rhs) = delete;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Methods to handle registering Fields and Neighbors
  void registerField(FieldBase<Dimension>& field) const;
  void unregisterField(FieldBase<Dimension>& field) const;
  void refreshView();
  void refreshMassView();
  void refreshMassRZView();
  void refreshPositionsView();
  void refreshVelocityView();
  void refreshHfieldView();
  void refreshWorkView();

  friend class FieldBase<Dimension>;

private:
  //--------------------------- Private Interface ---------------------------//
  std::string mName;

  // State fields.
  Field<Dimension, Scalar> mMass;
  Field<Dimension, Vector> mPositions;
  Field<Dimension, Vector> mVelocity;
  Field<Dimension, SymTensor> mH;

  // Optional fields.
  std::unique_ptr<Field<Dimension, Scalar>> mMassRZptr;

  // The work field is mutable.
  mutable Field<Dimension, Scalar> mWork;

  // List of Fields that are defined over this NodeList.
  mutable std::vector<std::reference_wrapper<FieldBase<Dimension>>> mFieldBases;
  Neighbor<Dimension>* mNeighborPtr;

  // Provide a dummy vector to buid NodeIterators against.
  std::vector<NodeList<Dimension>*> mDummyList;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "Field/Field.hh"

#endif
