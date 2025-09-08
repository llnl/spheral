#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Utilities/DBC.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Field/FieldList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Write Field as binary blob
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::write(const Field<Dimension, DataType>& value,
              const std::string path) {
  const auto buf = value.serialize();
  this->write(buf, path);
}

//------------------------------------------------------------------------------
// Read Field as binary blob
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::read(Field<Dimension, DataType>& value,
             const std::string path) const {
  vector<char> buf;
  this->read(buf, path);
  value.deserialize(buf);
}

//------------------------------------------------------------------------------
// Write a FieldList 
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::write(const FieldList<Dimension, DataType>& fieldList,
              const std::string path) {

  // Is the FieldList responsible for it's own memory?  If so, we have to 
  // provide additional information so it can properly restore itself.
  if (fieldList.storageType() == FieldStorageType::CopyFields) {
    std::vector<std::string> nodeListNames;
    for (const auto fieldPtr: fieldList) nodeListNames.push_back(fieldPtr->nodeList().name());
    this->write(nodeListNames, path + "/NodeListNames");
  }

  // Loop over each Field, and write each one using the descendent method.
  auto ifield = 0u;
  for (const auto fieldPtr: fieldList) {
    this->write(*fieldPtr, path + "/Field" + std::to_string(ifield++));
  }
}

//------------------------------------------------------------------------------
// Read a FieldList
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::read(FieldList<Dimension, DataType>& fieldList,
             const std::string path) const {

  // Is the FieldList responsible for it's own memory?  If so, we have to 
  // first make sure it has memory for each of the NodeLists it's defined against.
  if (fieldList.storageType() == FieldStorageType::CopyFields) {
    // We need the NodeListRegistrar.
    const NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
    const size_t numNodeLists = registrar.numNodeLists();
    CONTRACT_VAR(numNodeLists);
    const std::vector<std::string> registeredNames = registrar.registeredNames();

    // Read the set of NodeLists this FieldList is associated with.
    std::vector<std::string> nodeListNames;
    this->read(nodeListNames, path + "/NodeListNames");
    for (const auto& name: nodeListNames) {
      
      // Find the NodeList this is.
      const size_t nodeListi = std::distance(registeredNames.begin(), 
                                             find(registeredNames.begin(), registeredNames.end(), name));
      VERIFY(nodeListi < numNodeLists);
      const NodeList<Dimension>& nodeList = **(registrar.begin() + nodeListi);

      // If necessary, insert a field into this FieldList for this NodeList.
      if (not fieldList.haveNodeList(nodeList)) {
        fieldList.appendNewField("Unnamed Field", nodeList, DataTypeTraits<DataType>::zero());
      }
    }
  }

  // Loop over each Field, and read each one using the descendent method.
  auto ifield = 0u;
  for (auto fieldPtr: fieldList) {
    read(*fieldPtr, path + "/Field" + std::to_string(ifield++));
  }
}

//------------------------------------------------------------------------------
// Write a std::vector<DataType>.
//------------------------------------------------------------------------------
// Generic method
template<typename DataType>
inline
void
FileIO::write(const std::vector<DataType>& x, const std::string path) {
  std::vector<char> buf;
  packElement(x, buf);
  this->write_vector_char(buf, path);
}

// Specialize for some types that can be treated as arrays of doubles more efficiently/portably
//..............................................................................
template<typename Value>
inline
void
FileIO::writeVector(const std::vector<Value>& x, const std::string path) {
  const auto n = x.size();
  const auto ne = Value::numElements;
  std::vector<double> buf(n*ne);
  for (auto i = 0u; i < n; ++i) std::copy(x[i].begin(), x[i].end(), &buf[i*ne]);
  this->write(buf, path);
}
template<> inline void FileIO::write<Dim<1>::Vector>         (const std::vector<Dim<1>::Vector>& x,          const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<1>::Tensor>         (const std::vector<Dim<1>::Tensor>& x,          const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<1>::SymTensor>      (const std::vector<Dim<1>::SymTensor>& x,       const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<1>::ThirdRankTensor>(const std::vector<Dim<1>::ThirdRankTensor>& x, const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<2>::Vector>         (const std::vector<Dim<2>::Vector>& x,          const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<2>::Tensor>         (const std::vector<Dim<2>::Tensor>& x,          const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<2>::SymTensor>      (const std::vector<Dim<2>::SymTensor>& x,       const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<2>::ThirdRankTensor>(const std::vector<Dim<2>::ThirdRankTensor>& x, const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<3>::Vector>         (const std::vector<Dim<3>::Vector>& x,          const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<3>::Tensor>         (const std::vector<Dim<3>::Tensor>& x,          const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<3>::SymTensor>      (const std::vector<Dim<3>::SymTensor>& x,       const std::string path) { this->writeVector(x, path); }
template<> inline void FileIO::write<Dim<3>::ThirdRankTensor>(const std::vector<Dim<3>::ThirdRankTensor>& x, const std::string path) { this->writeVector(x, path); }

//------------------------------------------------------------------------------
// Read a std::vector<DataType>.
//------------------------------------------------------------------------------
// Generic method
template<typename DataType>
inline
void
FileIO::read(std::vector<DataType>& x, const std::string path) const {
  const auto buf = this->read_vector_char(path);
  auto itr = buf.begin();
  unpackElement(x, itr, buf.end());
  ENSURE(itr == buf.end());
}

// Specialize for some types that can be treated as arrays of doubles more efficiently/portably
//..............................................................................
template<typename Value>
inline
void
FileIO::readVector(std::vector<Value>& x, const std::string path) const {
  std::vector<double> buf;
  this->read(buf, path);
  const auto ne = Value::numElements;
  const auto n = buf.size()/ne;
  x.resize(n);
  for (auto i = 0u; i < n; ++i) std::copy(&buf[i*ne], &buf[i*ne] + ne, x[i].begin());
}
template<> inline void FileIO::read<Dim<1>::Vector>          (std::vector<Dim<1>::Vector>& x,          const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<1>::Tensor>          (std::vector<Dim<1>::Tensor>& x,          const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<1>::SymTensor>       (std::vector<Dim<1>::SymTensor>& x,       const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<1>::ThirdRankTensor> (std::vector<Dim<1>::ThirdRankTensor>& x, const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<2>::Vector>          (std::vector<Dim<2>::Vector>& x,          const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<2>::Tensor>          (std::vector<Dim<2>::Tensor>& x,          const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<2>::SymTensor>       (std::vector<Dim<2>::SymTensor>& x,       const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<2>::ThirdRankTensor> (std::vector<Dim<2>::ThirdRankTensor>& x, const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<3>::Vector>          (std::vector<Dim<3>::Vector>& x,          const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<3>::Tensor>          (std::vector<Dim<3>::Tensor>& x,          const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<3>::SymTensor>       (std::vector<Dim<3>::SymTensor>& x,       const std::string path) const { this->readVector(x, path); }
template<> inline void FileIO::read<Dim<3>::ThirdRankTensor> (std::vector<Dim<3>::ThirdRankTensor>& x, const std::string path) const { this->readVector(x, path); }

//------------------------------------------------------------------------------
// Safe method to try and read from a path if it exists.
// Returns: 0 => successful
//          1 => path does not exist
//          2 => unable to read value
//------------------------------------------------------------------------------
template<typename T>
inline
int
FileIO::readIfAvailable(T& value, const std::string path) const {
  if (not this->pathExists(path)) return 1;
  try {
    this->read(value, path);
  } catch(...) {
    return 2;
  }
  return 0;
}

}
