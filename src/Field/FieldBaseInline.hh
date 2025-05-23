#include "Utilities/DBC.hh"
#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given name.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::
FieldBase(typename FieldBase<Dimension>::FieldName name):
  mName(name),
  mNodeListPtr(0),
  mFieldListBaseList() {
}

//------------------------------------------------------------------------------
// Construct for a given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::
FieldBase(typename FieldBase<Dimension>::FieldName name,
          const NodeList<Dimension>& nodeList):
  mName(name),
  mNodeListPtr(&nodeList),
  mFieldListBaseList() {
  mNodeListPtr->registerField(*this);
}

//------------------------------------------------------------------------------
// Copy Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::FieldBase(const FieldBase& fieldBase):
  mName(fieldBase.name()),
  mNodeListPtr(fieldBase.nodeListPtr()),
  mFieldListBaseList() {
  mNodeListPtr->registerField(*this);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::~FieldBase() {
  if (mNodeListPtr) mNodeListPtr->unregisterField(*this);
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>&
FieldBase<Dimension>::operator=(const FieldBase<Dimension>& rhs) {
  if (this != &rhs) {
    if (mNodeListPtr) 
      mNodeListPtr->unregisterField(*this);
    mNodeListPtr = rhs.mNodeListPtr;
    mFieldListBaseList = std::vector<const FieldListBase<Dimension>*>();
    mNodeListPtr->registerField(*this);
  }
  return *this;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FieldBase<Dimension>::operator!=(const FieldBase<Dimension>& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Get the name string.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename FieldBase<Dimension>::FieldName
FieldBase<Dimension>::name() const {
  return mName;
}

//------------------------------------------------------------------------------
// Set the name string.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::name(typename FieldBase<Dimension>::FieldName name) {
  mName = name;
}

//------------------------------------------------------------------------------
// Return a reference to this field's NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
FieldBase<Dimension>::nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// Return a pointer to this field's NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>*
FieldBase<Dimension>::nodeListPtr() const {
  return mNodeListPtr;
}

//------------------------------------------------------------------------------
// Set the node list pointer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::unregisterNodeList() {
  mNodeListPtr->unregisterField(*this);
  mNodeListPtr = 0;
}

//------------------------------------------------------------------------------
// Set the node list pointer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::setFieldBaseNodeList(const NodeList<Dimension>& nodeList) {
  if (mNodeListPtr != 0) unregisterNodeList();
  mNodeListPtr = &nodeList;
  nodeList.registerField(*this);
}

//------------------------------------------------------------------------------
// Register a new FieldList as containing this Field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::
registerFieldList(const FieldListBase<Dimension>& fieldListBase) const {
  REQUIRE(!haveFieldList(fieldListBase));
  mFieldListBaseList.push_back(&fieldListBase);
}

//------------------------------------------------------------------------------
// Unregister a FieldList from this Field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::
unregisterFieldList(const FieldListBase<Dimension>& fieldListBase) const {
  typename std::vector<const FieldListBase<Dimension>*>::iterator itr = std::find(mFieldListBaseList.begin(),
                                                                                  mFieldListBaseList.end(),
                                                                                  &fieldListBase);
  REQUIRE(itr != mFieldListBaseList.end());
  mFieldListBaseList.erase(itr);
  ENSURE(!haveFieldList(fieldListBase));
}

//------------------------------------------------------------------------------
// Test if a given FieldList is registered with this Field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FieldBase<Dimension>::
haveFieldList(const FieldListBase<Dimension>& fieldListBase) const {
  return std::find(mFieldListBaseList.begin(),
                   mFieldListBaseList.end(),
                   &fieldListBase) != mFieldListBaseList.end();
}

}
