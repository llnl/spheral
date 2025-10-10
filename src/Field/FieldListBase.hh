//---------------------------------Spheral++----------------------------------//
// FieldListBase -- A base class to provide a generic handle on FieldLists.
//
// Created by JMO, Sun Apr 11 14:28:59 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListBase_hh__
#define __Spheral_FieldListBase_hh__

#include <vector>

namespace Spheral {

template<typename Dimension> class FieldBase;

template<typename Dimension>
class FieldListBase {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef FieldBase<Dimension>* ElementType;
  typedef FieldBase<Dimension>* value_type;    // STL compatibility
  typedef std::vector<ElementType> StorageType;

  typedef typename StorageType::iterator iterator;
  typedef typename StorageType::const_iterator const_iterator;
  typedef typename StorageType::reverse_iterator reverse_iterator;
  typedef typename StorageType::const_reverse_iterator const_reverse_iterator;

  // Constructors.
  FieldListBase();
  FieldListBase(const FieldListBase& rhs);

  // Destructor.
  virtual ~FieldListBase() = default;

  // Assignment operator.
  FieldListBase& operator=(const FieldListBase& rhs);

  // Require descendent types to fill in our iterators.
  virtual iterator begin_base() = 0;
  virtual iterator end_base() = 0;
  virtual reverse_iterator rbegin_base() = 0;
  virtual reverse_iterator rend_base() = 0;

  virtual const_iterator begin_base() const = 0;
  virtual const_iterator end_base() const = 0;
  virtual const_reverse_iterator rbegin_base() const = 0;
  virtual const_reverse_iterator rend_base() const = 0;

protected:
  //--------------------------- Protected Interface ---------------------------//
  mutable bool mNewCoarseNodes;
  mutable bool mNewRefineNodes;

  // Provide methods for the FieldList to register with its member Fields.
  void registerWithField(const FieldBase<Dimension>& fieldBase) const;
  void unregisterFromField(const FieldBase<Dimension>& fieldBase) const;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#include "FieldListBaseInline.hh"

#endif
