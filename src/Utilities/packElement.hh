//---------------------------------Spheral++----------------------------------//
// packElement
// Methods to pack/unpack sets of values into vector<char>.  Useful for 
// constructing MPI buffers of data.
//
// Created by J. Michael Owen, Tue Oct 10 22:16:29 PDT 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_packElement__
#define __Spheral_packElement__

#include "Geometry/Dimension.hh"
#include "Geometry/PolyClipperUtilities.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/DomainNode.hh"
#include "Utilities/uniform_random.hh"
#include "RK/RKCorrectionParams.hh"
#include "RK/RKCoefficients.hh"

#include <stdint.h>
#include <vector>
#include <map>
#include <set>
#include <iterator>
#include <string>
#include <tuple>

#ifdef USE_MPI
#include <mpi.h>
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

//------------------------------------------------------------------------------
// A standalone method to pack a given DataType into a vector
// The generic version assumes the DataType has a begin/end iterator pair, and
// a fixed size.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
packElement(const DataType& value,
            std::vector<char>& buffer) {
  const int packSize = sizeof(typename DataTypeTraits<DataType>::ElementType);
  for (typename DataType::const_iterator valueItr = value.begin();
       valueItr != value.end();
       ++valueItr) {
    char* data = reinterpret_cast<char*>(const_cast<typename DataTypeTraits<DataType>::ElementType*>(&(*valueItr)));
    for (int i = 0; i != packSize; ++i) {
      buffer.push_back(*(data + i));
    }
  }
}

// Specialization for an char type.
template<>
inline
void
packElement<char>(const char& value, 
                 std::vector<char>& buffer) {
  const int packSize = sizeof(char);
  char* data = reinterpret_cast<char*>(const_cast<char*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}

// Specialization for an int type.
template<>
inline
void
packElement<int>(const int& value, 
                 std::vector<char>& buffer) {
  const int packSize = sizeof(int);
  char* data = reinterpret_cast<char*>(const_cast<int*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}

// // Specialization for size_t.
// template<>
// inline
// void
// packElement<size_t>(const size_t& value, 
//                     std::vector<char>& buffer) {
//   const int packSize = sizeof(size_t);
//   char* data = reinterpret_cast<char*>(const_cast<size_t*>(&value));
//   for (int i = 0; i != packSize; ++i) {
//     buffer.push_back(*(data + i));
//   }
// }

// // Specialization for an unsigned int type.
// template<>
// inline
// void
// packElement<unsigned>(const unsigned& value, 
//                       std::vector<char>& buffer) {
//   const int packSize = sizeof(unsigned);
//   char* data = reinterpret_cast<char*>(const_cast<unsigned*>(&value));
//   for (int i = 0; i != packSize; ++i) {
//     buffer.push_back(*(data + i));
//   }
// }

// // Specialization for an unsigned long int type.
#if __APPLE__
template<>
inline
void
packElement<unsigned long>(const unsigned long& value, 
                           std::vector<char>& buffer) {
  const int packSize = sizeof(unsigned long);
  char* data = reinterpret_cast<char*>(const_cast<unsigned long*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}
#endif

// // Specialization for an unsigned long long int type.
// template<>
// inline
// void
// packElement<unsigned long long>(const unsigned& value, 
//                                 std::vector<char>& buffer) {
//   const int packSize = sizeof(unsigned long long);
//   char* data = reinterpret_cast<char*>(const_cast<unsigned long long*>(&value));
//   for (int i = 0; i != packSize; ++i) {
//     buffer.push_back(*(data + i));
//   }
// }

// Specialization for an uint32_t int type.
template<>
inline
void
packElement<uint32_t>(const uint32_t& value, 
                      std::vector<char>& buffer) {
  const int packSize = sizeof(uint32_t);
  char* data = reinterpret_cast<char*>(const_cast<uint32_t*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}

// Specialization for an uint64_t int type.
template<>
inline
void
packElement<uint64_t>(const uint64_t& value, 
                                std::vector<char>& buffer) {
  const int packSize = sizeof(uint64_t);
  char* data = reinterpret_cast<char*>(const_cast<uint64_t*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}

// Specialization for a scalar (float) type.
template<>
inline
void
packElement<float>(const float& value, 
                   std::vector<char>& buffer) {
  const int packSize = sizeof(float);
  char* data = reinterpret_cast<char*>(const_cast<float*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}

// Specialization for a scalar (double) type.
template<>
inline
void
packElement<double>(const double& value, 
                    std::vector<char>& buffer) {
  const int packSize = sizeof(double);
  char* data = reinterpret_cast<char*>(const_cast<double*>(&value));
  for (int i = 0; i != packSize; ++i) {
    buffer.push_back(*(data + i));
  }
}

// Specialization for a std::string.
template<>
inline
void
packElement<std::string>(const std::string& value, 
                         std::vector<char>& buffer) {
  const size_t size = value.size();
  packElement(size, buffer);
  buffer.insert(buffer.end(), value.begin(), value.begin() + size);
}

// Specialization for a std::pair of known types.
template<typename T1, typename T2>
inline
void
packElement(const std::pair<T1, T2>& value,
            std::vector<char>& buffer) {
  packElement(value.first, buffer);
  packElement(value.second, buffer);
}

// Specialization for a std::tuple of three common elements.
template<typename T>
inline
void
packElement(const std::tuple<T, T, T>& value,
            std::vector<char>& buffer) {
  packElement(std::get<0>(value), buffer);
  packElement(std::get<1>(value), buffer);
  packElement(std::get<2>(value), buffer);
}

// Specialization for a std::tuple of four common elements.
template<typename T>
inline
void
packElement(const std::tuple<T, T, T, T>& value,
            std::vector<char>& buffer) {
  packElement(std::get<0>(value), buffer);
  packElement(std::get<1>(value), buffer);
  packElement(std::get<2>(value), buffer);
  packElement(std::get<3>(value), buffer);
}

// Specialization for a std::tuple of five common elements.
template<typename T>
inline
void
packElement(const std::tuple<T, T, T, T, T>& value,
            std::vector<char>& buffer) {
  packElement(std::get<0>(value), buffer);
  packElement(std::get<1>(value), buffer);
  packElement(std::get<2>(value), buffer);
  packElement(std::get<3>(value), buffer);
  packElement(std::get<4>(value), buffer);
}

// Specialize for a std::vector<DataType>.
// Assumes the elements of the vector<> are of a type we already know how to pack.
template<typename DataType>
inline
void
packElement(const std::vector<DataType>& value, 
            std::vector<char>& buffer) {

  // First push the size of the vector onto the buffer.
  const unsigned size = value.size();
  packElement(size, buffer);

  // Now push each of the elements on.
  for (typename std::vector<DataType>::const_iterator itr = value.begin();
       itr != value.end();
       ++itr) {
    packElement(*itr, buffer);
  }
}

// Specialize for std::array<DataType, size>
template<typename DataType, std::size_t size>
inline
void
packElement(const std::array<DataType, size>& value, 
            std::vector<char>& buffer) {
  // We don't need to store size information like a vector
  // Push the elements on
  for (typename std::array<DataType, size>::const_iterator itr = value.begin();
       itr != value.end();
       ++itr) {
    packElement(*itr, buffer);
  }
}

// Specialize for a std::unordered_map<T1, T2, T3>
// Assumes that we know how to pack T1 and T2
template<typename T1, typename T2, typename T3>
inline
void
packElement(const std::unordered_map<T1, T2, T3>& value,
            std::vector<char>& buffer) {
  // Push the size of the unordered_map onto the buffer.
  const unsigned size = value.size();
  packElement(size, buffer);

  // Push each of the elements on
  for (typename std::unordered_map<T1, T2, T3>::const_iterator itr = value.begin();
       itr != value.end();
       ++itr) {
    packElement(itr->first, buffer);
    packElement(itr->second, buffer);
  }
}

// Specialize for a std::set<DataType>.
// Assumes the elements of the vector<> are of a type we already know how to pack.
template<typename DataType>
inline
void
packElement(const std::set<DataType>& value, 
            std::vector<char>& buffer) {

  // First push the size of the set onto the buffer.
  const unsigned size = value.size();
  packElement(size, buffer);

  // Now push each of the elements on.
  for (const auto& x: value) packElement(x, buffer);
}

// RKOrder
template<>
inline
void
packElement(const RKOrder& value,
            std::vector<char>& buffer) {
  packElement(static_cast<int>(value), buffer);
}

// RKCoefficients<Dimension>
template<typename Dimension>
inline
void
packElement(const RKCoefficients<Dimension>& value,
            std::vector<char>& buffer) {
  packElement(value.correctionOrder, buffer);
  packElement(value.coeffs, buffer);
}

// uniform_random
inline
void
packElement(const uniform_random& value,
            std::vector<char>& buffer) {
  value.serialize(buffer);
}

//------------------------------------------------------------------------------
// A standalone method to unpack an encoded value from the given
// vector<*>::iterator to the dataType.
// Note these methods increment the buffer iterator as they go, which alters 
// the value in the calling method.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
unpackElement(DataType& value, 
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(typename DataTypeTraits<DataType>::ElementType);
  for (typename DataType::iterator valueItr = value.begin();
       valueItr != value.end();
       ++valueItr) {
    volatile char* data = reinterpret_cast<char*>(&(*valueItr));
    for (int i = 0; i != packSize; ++i, ++itr) {
      CHECK(itr < endPackedVector);
      *(data + i) = *itr;
    }
  }
  ENSURE(itr <= endPackedVector);
}

// Specialization for an char type.
template<>
inline
void
unpackElement<char>(char& value,
                   std::vector<char>::const_iterator& itr,
                   const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(char);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}

// Specialization for an int type.
template<>
inline
void
unpackElement<int>(int& value,
                   std::vector<char>::const_iterator& itr,
                   const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(int);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}

// // Specialization for a size_t.
// template<>
// inline
// void
// unpackElement<size_t>(size_t& value,
//                       std::vector<char>::const_iterator& itr,
//                       const std::vector<char>::const_iterator& endPackedVector) {
//   const int packSize = sizeof(size_t);
//   char* data = reinterpret_cast<char*>(&value);
//   for (int i = 0; i != packSize; ++i, ++itr) {
//     CHECK(itr < endPackedVector);
//     *(data + i) = *itr;
//   }
//   ENSURE(itr <= endPackedVector);
// }

// // Specialization for an unsigned int type.
// template<>
// inline
// void
// unpackElement<unsigned>(unsigned& value,
//                         std::vector<char>::const_iterator& itr,
//                         const std::vector<char>::const_iterator& endPackedVector) {
//   const int packSize = sizeof(unsigned);
//   char* data = reinterpret_cast<char*>(&value);
//   for (int i = 0; i != packSize; ++i, ++itr) {
//     CHECK(itr < endPackedVector);
//     *(data + i) = *itr;
//   }
//   ENSURE(itr <= endPackedVector);
// }

// // Specialization for an unsigned long int type.
#if __APPLE__
template<>
inline
void
unpackElement<unsigned long>(unsigned long& value,
                        std::vector<char>::const_iterator& itr,
                        const std::vector<char>::const_iterator& endPackedVector) {
  const int packSize = sizeof(unsigned long);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}
#endif

// Specialization for an uint32_t int type.
template<>
inline
void
unpackElement<uint32_t>(uint32_t& value,
                        std::vector<char>::const_iterator& itr,
                        const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(uint32_t);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}

// Specialization for an uint64_t int type.
template<>
inline
void
unpackElement<uint64_t>(uint64_t& value,
                        std::vector<char>::const_iterator& itr,
                        const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(uint64_t);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}

// Specialization for a scalar (float) type.
template<>
inline
void
unpackElement<float>(float& value,
                     std::vector<char>::const_iterator& itr,
                     const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(float);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}

// Specialization for a scalar (double) type.
template<>
inline
void
unpackElement<double>(double& value,
                      std::vector<char>::const_iterator& itr,
                      const std::vector<char>::const_iterator& endPackedVector) {
  CONTRACT_VAR(endPackedVector);
  const int packSize = sizeof(double);
  char* data = reinterpret_cast<char*>(&value);
  for (int i = 0; i != packSize; ++i, ++itr) {
    CHECK(itr < endPackedVector);
    *(data + i) = *itr;
  }
  ENSURE(itr <= endPackedVector);
}

// Specialization for a std::string
template<>
inline
void
unpackElement<std::string>(std::string& value,
                           std::vector<char>::const_iterator& itr,
                           const std::vector<char>::const_iterator& endPackedVector) {
  size_t size;
  unpackElement(size, itr, endPackedVector);
  CHECK(itr + size <= endPackedVector);
  std::string result(itr, itr + size);
  value = result;
  itr += size;
}

// std::pair<T1,T2>
template<typename T1, typename T2>
inline
void
unpackElement(std::pair<T1, T2>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.first, itr, endPackedVector);
  unpackElement(value.second, itr, endPackedVector);
}

// std::tuple<T,T,T>
template<typename DataType>
inline
void
unpackElement(std::tuple<DataType, DataType, DataType>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  DataType x, y, z;
  unpackElement(x, itr, endPackedVector);
  unpackElement(y, itr, endPackedVector);
  unpackElement(z, itr, endPackedVector);
  std::get<0>(value) = x;
  std::get<1>(value) = y;
  std::get<2>(value) = z;
}

// std::tuple<T,T,T,T>
template<typename DataType>
inline
void
unpackElement(std::tuple<DataType, DataType, DataType, DataType>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  DataType x, y, z, a;
  unpackElement(x, itr, endPackedVector);
  unpackElement(y, itr, endPackedVector);
  unpackElement(z, itr, endPackedVector);
  unpackElement(a, itr, endPackedVector);
  std::get<0>(value) = x;
  std::get<1>(value) = y;
  std::get<2>(value) = z;
  std::get<3>(value) = a;
}

// std::tuple<T,T,T>
template<typename DataType>
inline
void
unpackElement(std::tuple<DataType, DataType, DataType, DataType, DataType>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  DataType x, y, z, a, b;
  unpackElement(x, itr, endPackedVector);
  unpackElement(y, itr, endPackedVector);
  unpackElement(z, itr, endPackedVector);
  unpackElement(a, itr, endPackedVector);
  unpackElement(b, itr, endPackedVector);
  std::get<0>(value) = x;
  std::get<1>(value) = y;
  std::get<2>(value) = z;
  std::get<3>(value) = a;
  std::get<4>(value) = b;
}

// Handle the vector<DataType> case, so long as DataType is one of the types
// we know how to unpack.
template<typename DataType>
inline
void
unpackElement(std::vector<DataType>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {

  // Read the size of the vector.
  unsigned int size;
  unpackElement(size, itr, endPackedVector);
  CHECK2(size <= std::distance(itr, endPackedVector),
         "Crazy buffer size:  " << size << " " << std::distance(itr, endPackedVector));

  // Now iterate over the number of elements we will unpack, and push them onto
  // the value.
  value.clear();
  for (unsigned int i = 0; i != size; ++i) {
    DataType element;
    unpackElement(element, itr, endPackedVector);
    value.push_back(element);
  }

  ENSURE(itr <= endPackedVector);
}

// std::array<DataType, size>
template<typename DataType, std::size_t size>
inline
void
unpackElement(std::array<DataType, size>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  // Unlike std::vector, we don't need the size data
  // Unpack the elements
  for (int i = 0; i != size; ++i) {
    DataType element;
    unpackElement(element, itr, endPackedVector);
    value[i] = element;
  }
}

// std::unordered_map<T1, T2, T3>, as long as we know how to unpack T1, T2
template<typename T1, typename T2, typename T3>
inline
void
unpackElement(std::unordered_map<T1, T2, T3>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  // Read the size of the map
  unsigned size;
  unpackElement(size, itr, endPackedVector);
  CHECK(2 * size <= std::distance(itr, endPackedVector));

  // Iterate over the number of elements to unpack and add them to the map
  value.clear();
  for (int i = 0; i < (int)size; ++i) {
    T1 key;
    T2 val;
    unpackElement(key, itr, endPackedVector);
    unpackElement(val, itr, endPackedVector);
    value[key] = val;
  }

  ENSURE(itr <= endPackedVector);
}

// Handle the set<DataType> case, so long as DataType is one of the types
// we know how to unpack.
template<typename DataType>
inline
void
unpackElement(std::set<DataType>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {

  // Read the size of the vector.
  unsigned int size;
  unpackElement(size, itr, endPackedVector);
  CHECK2(size <= std::distance(itr, endPackedVector),
         "Crazy buffer size:  " << size << " " << std::distance(itr, endPackedVector));

  // Now iterate over the number of elements we will unpack, and push them onto
  // the value.
  value.clear();
  for (unsigned int i = 0; i != size; ++i) {
    DataType element;
    unpackElement(element, itr, endPackedVector);
    value.insert(element);
  }

  ENSURE(itr <= endPackedVector);
}

// RKOrder
template<>
inline
void
unpackElement(RKOrder& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  int iorder;
  unpackElement(iorder, itr, endPackedVector);
  value = static_cast<RKOrder>(iorder);
  ENSURE(itr <= endPackedVector);
}

// RKCoeffients
template<typename Dimension>
inline
void
unpackElement(RKCoefficients<Dimension>& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.correctionOrder, itr, endPackedVector);
  unpackElement(value.coeffs, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

// uniform_random
inline
void
unpackElement(uniform_random& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  value.deserialize(itr, endPackedVector);
}

//------------------------------------------------------------------------------
// Compute the size of a buffer necessary to hold the given field data.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
int
computeBufferSize(const Field<Dimension, DataType>& /*field*/,
                  const std::vector<size_t>& packIndices,
                  const int /*sendProc*/,
                  const int /*recvProc*/) {
  return (packIndices.size() * 
          DataTypeTraits<DataType>::numElements(DataType()) * 
          sizeof(typename DataTypeTraits<DataType>::ElementType));
}

// Specialization for a std::vector
template<typename Dimension, typename DataType>
inline
int
computeBufferSize(const Field<Dimension, std::vector<DataType> >& field,
                  const std::vector<size_t>& packIndices,
                  const int sendProc,
                  const int recvProc) {

  // Byte size of the primitive type of the vector.
  const int elementSize = (DataTypeTraits<DataType>::numElements(DataType()) * 
                           sizeof(typename DataTypeTraits<DataType>::ElementType));

  // Find the rank of this processor.
  int rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(Communicator::communicator(), &rank);
#else
  CONTRACT_VAR(recvProc);
#endif
  REQUIRE(rank == sendProc || rank == recvProc);

  // The send proc can compute the necessary size.
  int bufSize = 0;
  if (rank == sendProc) {
    for (auto i: packIndices) {
      bufSize += field(i).size();
    }
    bufSize *= elementSize;
  }

  // Communicate the result to the receiving processor.
#ifdef USE_MPI
  if (rank == recvProc) {
    MPI_Status status;
    MPI_Recv(&bufSize, 1, MPI_INT, sendProc, 103, Communicator::communicator(), &status);
  } else if (rank == sendProc) {
    MPI_Send(&bufSize, 1, MPI_INT, recvProc, 103, Communicator::communicator());
  }
#endif
  return bufSize;
}

//------------------------------------------------------------------------------
// Pack the given Field values into a 1-D array for transmission.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<char>
packFieldValues(const Field<Dimension, DataType>& field,
                const std::vector<size_t>& packIndices) {

  // Prepare the return vector.
  std::vector<char> result;

  // Loop over the elements of the Field we are packing, and push the 
  // packed elements onto the result.
  for (auto i: packIndices) {
    CHECK(i < field.numElements());
    packElement(field(i), result);
  }

  return result;
}

// Same thing for all internal values
template<typename Dimension, typename DataType>
inline
std::vector<char>
packFieldValues(const Field<Dimension, DataType>& field) {
  std::vector<char> result;
  const auto n = field.numInternalElements();
  for (auto i = 0u; i < n; ++i) packElement(field[i], result);
  return result;
}

//------------------------------------------------------------------------------
// Unpack encoded values from the given vector to the Field at the indicated
// indices.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
unpackFieldValues(Field<Dimension, DataType>& field,
                  const std::vector<size_t>& packIndices,
                  const std::vector<char>& packedValues) {

  // Loop over the elements of the Field we are unpacking.
  auto bufItr = packedValues.begin();
  for (auto i: packIndices) {
    CHECK(i < field.numElements());
    CHECK(bufItr < packedValues.end());
    unpackElement(field(i), bufItr, packedValues.end());
    CHECK(bufItr <= packedValues.end());
  }

  ENSURE(bufItr == packedValues.end());
}

// Same thing for all internal values
template<typename Dimension, typename DataType>
inline
void
unpackFieldValues(Field<Dimension, DataType>& field,
                  const std::vector<char>& packedValues) {
  auto bufItr = packedValues.begin();
  auto endItr = packedValues.end();
  const auto n = field.numInternalElements();
  for (auto i = 0u; i < n; ++i) {
    unpackElement(field[i], bufItr, endItr);
    CHECK(bufItr <= endItr);
  }
  CHECK(bufItr == endItr);
}

//------------------------------------------------------------------------------
// Pack the Volume types.
//------------------------------------------------------------------------------
template<>
inline
void
packElement(const Dim<1>::Box& value,
            std::vector<char>& buffer) {
  packElement(value.center(), buffer);
  packElement(value.extent(), buffer);
}

template<>
inline
void
packElement(const GeomPolygon& value,
            std::vector<char>& buffer) {
  packElement(value.vertices(), buffer);
  packElement(value.facetVertices(), buffer);
}

template<>
inline
void
packElement(const GeomPolyhedron& value,
            std::vector<char>& buffer) {
  packElement(value.vertices(), buffer);
  packElement(value.facetVertices(), buffer);
}

//------------------------------------------------------------------------------
// Unpack the Volume types.
//------------------------------------------------------------------------------
template<>
inline
void
unpackElement(Dim<1>::Box& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  Dim<1>::Vector center;
  double extent;
  unpackElement(center, itr, endPackedVector);
  unpackElement(extent, itr, endPackedVector);
  value.center(center);
  value.extent(extent);
}              

template<>
inline
void
unpackElement(GeomPolygon& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  std::vector<Dim<2>::Vector> vertices;
  std::vector<std::vector<unsigned> > facetVertices;
  unpackElement(vertices, itr, endPackedVector);
  unpackElement(facetVertices, itr, endPackedVector);
  value.reconstruct(vertices, facetVertices);
}

template<>
inline
void
unpackElement(GeomPolyhedron& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  std::vector<Dim<3>::Vector> vertices;
  std::vector<std::vector<unsigned> > facetVertices;
  unpackElement(vertices, itr, endPackedVector);
  unpackElement(facetVertices, itr, endPackedVector);
  value.reconstruct(vertices, facetVertices);
}

//------------------------------------------------------------------------------
// We also provide the ability to pack/unpack std::map's of known types.
//------------------------------------------------------------------------------
template<typename Key, typename Value>
inline
void
packElement(const std::map<Key, Value>& mapvalue,
            std::vector<char>& buffer) {

  // Pack the size first.
  packElement(unsigned(mapvalue.size()), buffer);

  // Now walk the key, value pairs and pack them up.
  for (typename std::map<Key, Value>::const_iterator itr = mapvalue.begin();
       itr != mapvalue.end();
       ++itr) {
    packElement(itr->first, buffer);
    packElement(itr->second, buffer);
  }
}

template<typename Key, typename Value>
inline
void
unpackElement(std::map<Key, Value>& mapvalue,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {

  // Get the size.
  unsigned numElements;
  unpackElement(numElements, itr, endPackedVector);

  // Unpack the individual (key, value) pairs and put them in.
  size_t i;
  Key key;
  Value value;
  for (i = 0; i != numElements; ++i) {
    unpackElement(key, itr, endPackedVector);
    unpackElement(value, itr, endPackedVector);
    mapvalue[key] = value;
  }
  ENSURE(itr <= endPackedVector);
}

//------------------------------------------------------------------------------
// Handle the PolyClipper types.
//------------------------------------------------------------------------------
// PolyClipper::Vertex2d
template<>
inline
void
packElement<PolyClipperVertex2d>(const PolyClipperVertex2d& value, 
                                 std::vector<char>& buffer) {
  packElement(value.position, buffer);
  packElement(value.neighbors, buffer);
  packElement(value.comp, buffer);
  packElement(value.ID, buffer);
  packElement(value.clips, buffer);
}

template<>
inline
void
unpackElement<PolyClipperVertex2d>(PolyClipperVertex2d& value, 
                                   std::vector<char>::const_iterator& itr,
                                   const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.position, itr, endPackedVector);
  unpackElement(value.neighbors, itr, endPackedVector);
  unpackElement(value.comp, itr, endPackedVector);
  unpackElement(value.ID, itr, endPackedVector);
  unpackElement(value.clips, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

//..............................................................................
// PolyClipper::Vertex3d
template<>
inline
void
packElement<PolyClipperVertex3d>(const PolyClipperVertex3d& value, 
                                 std::vector<char>& buffer) {
  packElement(value.position, buffer);
  packElement(value.neighbors, buffer);
  packElement(value.comp, buffer);
  packElement(value.ID, buffer);
  packElement(value.clips, buffer);
}

template<>
inline
void
unpackElement<PolyClipperVertex3d>(PolyClipperVertex3d& value, 
                                   std::vector<char>::const_iterator& itr,
                                   const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.position, itr, endPackedVector);
  unpackElement(value.neighbors, itr, endPackedVector);
  unpackElement(value.comp, itr, endPackedVector);
  unpackElement(value.ID, itr, endPackedVector);
  unpackElement(value.clips, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

//..............................................................................
// PolyClipper::Plane2d
template<>
inline
void
packElement<PolyClipperPlane2d>(const PolyClipperPlane2d& value, 
                                std::vector<char>& buffer) {
  packElement(value.dist, buffer);
  packElement(value.normal, buffer);
  packElement(value.ID, buffer);
}

template<>
inline
void
unpackElement<PolyClipperPlane2d>(PolyClipperPlane2d& value, 
                                  std::vector<char>::const_iterator& itr,
                                  const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.dist, itr, endPackedVector);
  unpackElement(value.normal, itr, endPackedVector);
  unpackElement(value.ID, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

//..............................................................................
// PolyClipper::Plane3d
template<>
inline
void
packElement<PolyClipperPlane3d>(const PolyClipperPlane3d& value, 
                                std::vector<char>& buffer) {
  packElement(value.dist, buffer);
  packElement(value.normal, buffer);
  packElement(value.ID, buffer);
}

template<>
inline
void
unpackElement<PolyClipperPlane3d>(PolyClipperPlane3d& value, 
                                  std::vector<char>::const_iterator& itr,
                                  const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.dist, itr, endPackedVector);
  unpackElement(value.normal, itr, endPackedVector);
  unpackElement(value.ID, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

//------------------------------------------------------------------------------
// DomainNode
//------------------------------------------------------------------------------
template<>
inline
void
packElement<DomainNode<Dim<1>>>(const DomainNode<Dim<1>>& value, 
                                std::vector<char>& buffer) {
  packElement(value.localNodeID, buffer);
  packElement(value.uniqueLocalNodeID, buffer);
  packElement(value.globalNodeID, buffer);
  packElement(value.nodeListID, buffer);
  packElement(value.domainID, buffer);
  packElement(value.work, buffer);
  packElement(value.position, buffer);
}

template<>
inline
void
packElement<DomainNode<Dim<2>>>(const DomainNode<Dim<2>>& value, 
                                std::vector<char>& buffer) {
  packElement(value.localNodeID, buffer);
  packElement(value.uniqueLocalNodeID, buffer);
  packElement(value.globalNodeID, buffer);
  packElement(value.nodeListID, buffer);
  packElement(value.domainID, buffer);
  packElement(value.work, buffer);
  packElement(value.position, buffer);
}

template<>
inline
void
packElement<DomainNode<Dim<3>>>(const DomainNode<Dim<3>>& value, 
                                std::vector<char>& buffer) {
  packElement(value.localNodeID, buffer);
  packElement(value.uniqueLocalNodeID, buffer);
  packElement(value.globalNodeID, buffer);
  packElement(value.nodeListID, buffer);
  packElement(value.domainID, buffer);
  packElement(value.work, buffer);
  packElement(value.position, buffer);
}

template<>
inline
void
unpackElement<DomainNode<Dim<1>>>(DomainNode<Dim<1>>& value, 
                                  std::vector<char>::const_iterator& itr,
                                  const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.localNodeID, itr, endPackedVector);
  unpackElement(value.uniqueLocalNodeID, itr, endPackedVector);
  unpackElement(value.globalNodeID, itr, endPackedVector);
  unpackElement(value.nodeListID, itr, endPackedVector);
  unpackElement(value.domainID, itr, endPackedVector);
  unpackElement(value.work, itr, endPackedVector);
  unpackElement(value.position, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}
template<>
inline
void
unpackElement<DomainNode<Dim<2>>>(DomainNode<Dim<2>>& value, 
                                  std::vector<char>::const_iterator& itr,
                                  const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.localNodeID, itr, endPackedVector);
  unpackElement(value.uniqueLocalNodeID, itr, endPackedVector);
  unpackElement(value.globalNodeID, itr, endPackedVector);
  unpackElement(value.nodeListID, itr, endPackedVector);
  unpackElement(value.domainID, itr, endPackedVector);
  unpackElement(value.work, itr, endPackedVector);
  unpackElement(value.position, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}
template<>
inline
void
unpackElement<DomainNode<Dim<3>>>(DomainNode<Dim<3>>& value, 
                                  std::vector<char>::const_iterator& itr,
                                  const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.localNodeID, itr, endPackedVector);
  unpackElement(value.uniqueLocalNodeID, itr, endPackedVector);
  unpackElement(value.globalNodeID, itr, endPackedVector);
  unpackElement(value.nodeListID, itr, endPackedVector);
  unpackElement(value.domainID, itr, endPackedVector);
  unpackElement(value.work, itr, endPackedVector);
  unpackElement(value.position, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

}

#endif
