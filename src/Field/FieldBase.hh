//---------------------------------Spheral++----------------------------------//
// FieldBase -- An abstract base class to provide a generic handle on Fields.
//
// Created by JMO, Tue Sep 14 22:36:19 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldBase_hh__
#define __Spheral_FieldBase_hh__

#include <string>
#include <vector>
#include <memory>

namespace Spheral {

// Forward declarations
template<typename Dimension> class NodeList;
template<typename Dimension> class FieldListBase;

template<typename Dimension>
class FieldBase {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename std::string FieldName;
  typedef typename Dimension::Scalar Scalar;

  // Constructors.
  FieldBase(FieldName name);
  FieldBase(FieldName name,
            const NodeList<Dimension>& nodeList);
  FieldBase(const FieldBase& fieldBase);
  virtual std::shared_ptr<FieldBase> clone() const = 0;

  // Destructor.
  virtual ~FieldBase();

  // Assignment operator.
  virtual FieldBase& operator=(const FieldBase& rhs);

  // Require descendent fields be able to test equivalence.
  virtual bool operator==(const FieldBase& rhs) const = 0;
  bool operator!=(const FieldBase& rhs) const;

  // Access the name.
  FieldName name() const;
  void name(FieldName name);

  // Provide methods to access and set the NodeList.
  const NodeList<Dimension>& nodeList() const;
  const NodeList<Dimension>* nodeListPtr() const;
  void unregisterNodeList();

  // Methods every field must provide.
  virtual unsigned size() const = 0;
  virtual void Zero() = 0;
  virtual void setNodeList(const NodeList<Dimension>& nodeList) = 0;
  virtual std::vector<char> packValues(const std::vector<int>& nodeIDs) const = 0;
  virtual void unpackValues(const std::vector<int>& nodeIDs,
                            const std::vector<char>& buffer) = 0;
  virtual void copyElements(const std::vector<int>& fromIndices,
                            const std::vector<int>& toIndices) = 0;
  virtual bool fixedSizeDataType() const = 0;
  virtual int numValsInDataType() const = 0;
  virtual int sizeofDataType() const = 0;
  virtual int computeCommBufferSize(const std::vector<int>& packIndices,
                                    const int sendProc,
                                    const int recvProc) const = 0;

//   // Methods to support cacheing of coarse and refine neighbor values.
//   void notifyNewCoarseNodes() const;
//   void notifyNewRefineNodes() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  void setFieldBaseNodeList(const NodeList<Dimension>& nodeListPtr);

  friend class NodeList<Dimension>;
  virtual void resizeField(unsigned size) = 0;
  virtual void resizeFieldInternal(unsigned size, unsigned oldFirstGhostNode) = 0;
  virtual void resizeFieldGhost(unsigned size) = 0;
  virtual void deleteElement(int nodeID) = 0;
  virtual void deleteElements(const std::vector<int>& nodeIDs) = 0;

  // Make the FieldListBase a friend, so that it can use the registration
  // methods.
  friend class FieldListBase<Dimension>;
  void registerFieldList(const FieldListBase<Dimension>& fieldList) const;
  void unregisterFieldList(const FieldListBase<Dimension>& fieldList) const;
  bool haveFieldList(const FieldListBase<Dimension>& fieldList) const;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldName mName;
  const NodeList<Dimension>* mNodeListPtr;

  // The set of FieldLists currently referencing this Field.
  mutable std::vector<const FieldListBase<Dimension>*> mFieldListBaseList;

  // Disallow the default constructor.
  FieldBase();
};

}

#include "FieldBaseInline.hh"

#endif
