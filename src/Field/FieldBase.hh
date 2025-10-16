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

template<typename Dimension>
class FieldBase {

public:
  //--------------------------- Public Interface ---------------------------//
  using FieldName = std::string;
  using Scalar = typename Dimension::Scalar;

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
  virtual size_t size() const = 0;
  virtual void Zero() = 0;
  virtual void setNodeList(const NodeList<Dimension>& nodeList) = 0;
  virtual std::vector<char> packValues(const std::vector<size_t>& nodeIDs) const = 0;
  virtual void unpackValues(const std::vector<size_t>& nodeIDs,
                            const std::vector<char>& buffer) = 0;
  virtual void copyElements(const std::vector<size_t>& fromIndices,
                            const std::vector<size_t>& toIndices) = 0;
  virtual bool fixedSizeDataType() const = 0;
  virtual size_t numValsInDataType() const = 0;
  virtual size_t sizeofDataType() const = 0;
  virtual size_t computeCommBufferSize(const std::vector<size_t>& packIndices,
                                       const int sendProc,
                                       const int recvProc) const = 0;

//   // Methods to support cacheing of coarse and refine neighbor values.
//   void notifyNewCoarseNodes() const;
//   void notifyNewRefineNodes() const;

  // Disallow the default constructor.
  FieldBase() = delete;

protected:
  //--------------------------- Protected Interface ---------------------------//
  void setFieldBaseNodeList(const NodeList<Dimension>& nodeListPtr);

  friend class NodeList<Dimension>;
  virtual void resizeField(size_t size) = 0;
  virtual void resizeFieldInternal(size_t size, size_t oldFirstGhostNode) = 0;
  virtual void resizeFieldGhost(size_t size) = 0;
  virtual void deleteElement(size_t nodeID) = 0;
  virtual void deleteElements(const std::vector<size_t>& nodeIDs) = 0;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldName mName;
  const NodeList<Dimension>* mNodeListPtr;
};

}

#include "FieldBaseInline.hh"

#endif
