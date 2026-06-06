//---------------------------------Spheral++----------------------------------//
// StateBase -- The base class for State and StateDerivatives method, providing
// the common methods for registering and storing Fields, as well as serving
// up FieldLists based on these Fields.
//
// Important usage notes!
// Copying StateBase objects (either by copy constructor or assignment) always
// copies by reference!  In other words, even if you have a StateBase object
// that is using internal storage to hold it's Fields, copying or assigning 
// another StateBase object with this one will only create new references to the
// original's Fields.  You have to explicitly call "copyFields()" if you want
// to do the expensive copy.
// 
// Created by JMO, Wed Aug 25 22:23:35 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StateBase_hh__
#define __Spheral_StateBase_hh__

#include "Field/FieldBase.hh"

#include <any>
#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <map>
#include <list>
#include <set>

namespace Spheral {

// Forward declaration.
template<typename Dimension> class NodeList;
template<typename Dimension> class FieldBase;
template<typename Dimension> class FieldListBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class Mesh;
template<typename Dimension> class ReproducingKernel;
enum class RKOrder : int;

template<typename Dimension>
class StateBase {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;
  using ConnectivityMapType = typename Spheral::ConnectivityMap<Dimension>;
  using MeshType = typename Spheral::Mesh<Dimension>;

  using ConnectivityMapPtr = std::shared_ptr<ConnectivityMapType>;
  using MeshPtr = std::shared_ptr<MeshType>;

  using KeyType = std::string;
  using FieldName = typename FieldBase<Dimension>::FieldName;

  // Constructors, destructor.
  StateBase();
  StateBase(const StateBase& rhs) = default;
  StateBase& operator=(const StateBase& rhs) = default;
  virtual ~StateBase() = default;

  // Test if two StateBases have equivalent fields.
  virtual bool operator==(const StateBase& rhs) const;

  //............................................................................
  // Enroll state
  virtual              void enroll(FieldBase<Dimension>& field);
  virtual              void enroll(std::shared_ptr<FieldBase<Dimension>>& fieldPtr);
  virtual              void enroll(FieldListBase<Dimension>& fieldList);
  template<typename T> void enroll(const KeyType& key, T& thing);
  template<typename T> void enrollCopy(const KeyType& key, T& thing);

  //............................................................................
  // Access Fields
  template<typename Value> Field<Dimension, Value>& field(const KeyType& key) const;
  template<typename Value> Field<Dimension, Value>& field(const KeyType& key,
                                                          const Value& dummy) const;

  // Get all registered fields of the given data type
  template<typename Value> std::vector<Field<Dimension, Value>*> allFields() const;
  template<typename Value> std::vector<Field<Dimension, Value>*> allFields(const Value& dummy) const;

  //............................................................................
  // Access FieldLists
  template<typename Value> FieldList<Dimension, Value> fields(const std::string& name,
                                                              bool allowNone = false) const;
  template<typename Value> FieldList<Dimension, Value> fields(const std::string& name, 
                                                              const Value& dummy,
                                                              bool allowNone = false) const;

  //............................................................................
  // Access an arbitrary type
  template<typename Value> Value& get(const KeyType& key) const;
  template<typename Value> Value& get(const KeyType& key, const Value& dummy) const;

  //............................................................................
  // Access an arbitrary type as a pointer.
  // Returns nullptr on failure
  template<typename Value> Value* getPtr(const KeyType& key) const;
  template<typename Value> Value* getPtr(const KeyType& key, const Value& dummy) const;

  //............................................................................
  // Test if the specified Field or key is currently registered.
  bool registered(const KeyType& key) const;
  bool registered(const FieldBase<Dimension>& field) const;
  bool registered(const FieldListBase<Dimension>& fieldList) const;
  bool fieldNameRegistered(const FieldName& fieldName) const;

  //............................................................................
  // Return the complete set of keys registered
  std::vector<KeyType> keys() const;

  // The field keys including mangling with NodeList names
  std::vector<KeyType> fullFieldKeys() const;

  // The non-field (miscellaneous) keys
  std::vector<KeyType> miscKeys() const;

  // Return the set of known field names (unencoded from our internal mangling
  // convention with the NodeList name).
  std::vector<FieldName> fieldNames() const;

  //............................................................................
  // A state object can carry around a reference to a ConnectivityMap.
  void enrollConnectivityMap(ConnectivityMapPtr connectivityMapPtr);
  const ConnectivityMapType& connectivityMap() const;
  ConnectivityMapType& connectivityMap();

  //............................................................................
  // We also provide support for registering a mesh (though only one per StateBase).
  void enrollMesh(MeshPtr meshPtr);
  bool meshRegistered() const;
  const MeshType& mesh() const;
  MeshType& mesh();

  //............................................................................
  // Set the Fields equal to those in another State object.
  void assign(const StateBase<Dimension>& rhs);

  // Assign just the fields with the given name to those in another State object.
  template<typename Value>
  void assignFields(const StateBase<Dimension>& rhs, const std::string name);

  // Force the StateBase to create new internally owned copies of all state.
  virtual void copyState();

  //............................................................................
  // Construct the lookup key for the given field.
  static KeyType key(const FieldBase<Dimension>& field);

  // Construct the lookup key for the given FieldList.
  static KeyType key(const FieldListBase<Dimension>& fieldList);

  // Encode our underlying convention for combining Field and NodeList names into a 
  // single Key.
  static KeyType buildFieldKey(const KeyType& fieldName, const KeyType& nodeListName);
  static void splitFieldKey(const KeyType& key, KeyType& fieldName, KeyType& nodeListName);

protected:
  //--------------------------- Protected Interface ---------------------------//
  using StorageType = std::map<KeyType, std::any>;
  using CacheType = std::list<std::any>;

  // Protected data.
  StorageType     mStorage;
  CacheType       mCache;
  std::set<const NodeList<Dimension>*> mNodeListPtrs;
  ConnectivityMapPtr mConnectivityMapPtr;
  MeshPtr mMeshPtr;
};

}

#include "StateBaseInline.hh"

#endif

