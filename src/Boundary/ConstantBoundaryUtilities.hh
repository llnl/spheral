#ifndef __Spheral_ConstantBoundaryUtilities__
#define __Spheral_ConstantBoundaryUtilities__

namespace Spheral {

//------------------------------------------------------------------------------
// Store the Field values for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
storeFieldValues(const NodeList<Dimension>& nodeList,
                 const std::vector<size_t>& nodeIDs,
                 std::map<std::string, std::vector<char>>& values) {
  // std::cerr << "storeFieldValues starting size: " << values.size() << std::endl;
  for (auto fieldItr = nodeList.registeredFieldsBegin();
       fieldItr != nodeList.registeredFieldsEnd();
       ++fieldItr) {
    const auto buffer = (**fieldItr).packValues(nodeIDs);
    const auto key = StateBase<Dimension>::key(**fieldItr);
    // std::cerr << "Storing key " << key << std::endl;
    // if (values.find(key) != values.end()) std::cerr << "ConstantBoundaryUtilities::storeFieldValues collision for key " << key << std::endl;
    values[key] = buffer;
    // std::cerr << "    Stored " << buffer.size() << " values for " << key << std::endl;
  }
}
  
//------------------------------------------------------------------------------
// Set the Field values in the given field using the given map.
//------------------------------------------------------------------------------
template<typename Dimension>
void
resetValues(FieldBase<Dimension>& field,
            const std::vector<size_t>& nodeIDs,
            const std::map<std::string, std::vector<char>>& values,
            const bool dieOnMissingField) {
  // Find this Field in the set of stored values.
  const auto key = StateBase<Dimension>::key(field);
  auto itr = values.find(key);
  VERIFY2(itr != values.end() or not dieOnMissingField,
          "ConstantBoundary error: " << key << " not found in stored field values.");

  // Now set the values.
  if (itr != values.end()) {
    const auto& buffer = itr->second;
    field.unpackValues(nodeIDs, buffer);
  }
}

//------------------------------------------------------------------------------
// Copy all fields from -> to by IDs
//------------------------------------------------------------------------------
template<typename Dimension>
void
copyFieldValues(const NodeList<Dimension>& nodeList,
                const std::vector<size_t>& fromIDs,
                const std::vector<size_t>& toIDs) {
  for (auto fieldItr = nodeList.registeredFieldsBegin();
       fieldItr != nodeList.registeredFieldsEnd();
       ++fieldItr) {
    (**fieldItr).copyElements(fromIDs, toIDs);
  }
}

//------------------------------------------------------------------------------
// Extract the encoded values
//------------------------------------------------------------------------------
template<typename Value>
std::vector<Value>
extractBufferedValues(const std::vector<char>& buffer) {
  std::vector<Value> result;
  auto itr = buffer.begin();
  auto endbuf = buffer.end();
  auto n = 0;
  while (itr < endbuf) {
    result.resize(++n);
    unpackElement(result.back(), itr, endbuf);
  }
  return result;
}

}

#endif
