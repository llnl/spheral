//---------------------------------Spheral++----------------------------------//
// DataTypeTraits.  A simple trait class to define the number of independent
// elements in the given DataType.
//
// Created by J. Michael Owen, Thu Oct 18 16:36:21 PDT 2001
//----------------------------------------------------------------------------//
#ifndef DataTypeTraits_HH
#define DataTypeTraits_HH

#include <climits>
#include <stdint.h>
#include <vector>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/PolyClipperUtilities.hh"
#include "Distributed/RegisterMPIDataTypes.hh"
#include "Utilities/DomainNode.hh"
#include "RK/RKCorrectionParams.hh"
#include "RK/RKCoefficients.hh"
#include "axom/sidre.hpp"
#include "Utilities/uniform_random.hh"

#ifdef SPHERAL_ENABLE_MPI
extern "C" {
#include <mpi.h>
#include <stdint.h>
#include <limits.h>

// Macro to help us figure out the best type for size_t
#if SIZE_MAX == UCHAR_MAX
   #define SPHERAL_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define SPHERAL_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define SPHERAL_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define SPHERAL_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define SPHERAL_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "Unable to figure out MPI type for size_t"
#endif
}
#endif

namespace Spheral {

template<typename DataType> struct DataTypeTraits {using AxomType = double;};

// template <typename T> struct is_rank_n_tensor : std::false_type {};
// template <int U> struct is_rank_n_tensor<GeomThirdRankTensor<U>> : std::true_type {};
// template <int U> struct is_rank_n_tensor<GeomFourthRankTensor<U>> : std::true_type {};
// template <int U> struct is_rank_n_tensor<GeomFifthRankTensor<U>> : std::true_type {};

//------------------------------------------------------------------------------
// Provide a trait test if type provides operator<
namespace TypeTraits {
  template<typename T1, typename T2, typename = void>
  struct has_less_than : std::false_type {};

  template<typename T1, typename T2>
  struct has_less_than<T1, T2, std::void_t<decltype(std::declval<T1>() < std::declval<T2>())>> : std::true_type {};
}

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<bool> {
  typedef bool ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  static bool zero() { return false; }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return MPI_C_BOOL; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::INT8_ID; }
  using AxomType = bool;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<char> {
  typedef char ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  SPHERAL_HOST_DEVICE static int zero() { return '\0'; }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return MPI_CHAR; }
#endif

#if (CHAR_MIN==0)
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::UINT8_ID; }
#else
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::INT8_ID; }
#endif

  using AxomType = char;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<int> {
  typedef int ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  SPHERAL_HOST_DEVICE static int zero() { return 0; }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return MPI_INT; }
  static MPI_Datatype MpiLocDataType() { return MPI_2INT; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::INT_ID; }
  using AxomType = int;
};

//------------------------------------------------------------------------------
// template<>
// struct DataTypeTraits<size_t> {
//   typedef size_t ElementType;
//   static bool fixedSize() { return true; }
//   static size_t numElements(const ElementType& x) { return 1; }
//   static size_t zero() { return 0U; }
// #ifdef SPHERAL_ENABLE_MPI
//   static MPI_Datatype MpiDataType() { return SPHERAL_MPI_SIZE_T; }
// #endif
//   static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::UINT64_ID; }
//   using AxomType = uint64_t;
// };

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<uint32_t> {
  typedef uint32_t ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  static uint32_t zero() { return 0UL; }
#ifdef SPHERAL_ENABLE_MPI
#ifdef MPI_UINT32_T
  static MPI_Datatype MpiDataType() { return MPI_UINT32_T; }
#else
  static MPI_Datatype MpiDataType() { return MPI_UNSIGNED; }
#endif
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::UINT32_ID; }
  using AxomType = uint32_t;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<uint64_t> {
  typedef uint64_t ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  SPHERAL_HOST_DEVICE static uint64_t zero() { return 0ULL; }
#ifdef SPHERAL_ENABLE_MPI
#ifdef MPI_UINT64_T
  static MPI_Datatype MpiDataType() { return MPI_UINT64_T; }
#else
  static MPI_Datatype MpiDataType() { return MPI_UNSIGNED_LONG_LONG; }
#endif
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::UINT64_ID; }
  using AxomType = uint64_t;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<float> {
  typedef float ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  SPHERAL_HOST_DEVICE static float zero() { return 0.0; }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return MPI_FLOAT; }
  static MPI_Datatype MpiLocDataType() { return MPI_FLOAT_INT; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::FLOAT_ID; }
  using AxomType = float;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<double> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType&) { return 1; }
  SPHERAL_HOST_DEVICE static double zero() { return 0.0; }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return MPI_DOUBLE; }
  static MPI_Datatype MpiLocDataType() { return MPI_DOUBLE_INT; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<std::string> {
  typedef std::string ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const ElementType& x) { return x.size(); }
  static ElementType zero() { return ""; }

#if (CHAR_MIN==0)
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::UINT8_ID; }
#else
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::INT8_ID; }
#endif

  using AxomType = char;
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<std::vector<Value> > {
  typedef std::vector<Value> ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const std::vector<Value>& x) { return x.size(); }
  static std::vector<Value> zero() { return std::vector<Value>(); }

  static axom::sidre::DataTypeId axomTypeID() { return DataTypeTraits<Value>::axomTypeID(); }
  using AxomType = typename DataTypeTraits<Value>::AxomType;
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<std::set<Value> > {
  typedef std::set<Value> ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const std::set<Value>& x) { return x.size(); }
  static std::set<Value> zero() { return std::set<Value>(); }
  using AxomType = typename DataTypeTraits<Value>::AxomType;
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<std::tuple<Value, Value, Value> > {
  typedef Value ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const std::tuple<Value, Value, Value>&) { return 3; }
  static std::tuple<Value, Value, Value> zero() { return std::make_tuple(DataTypeTraits<Value>::zero(), 
                                                                         DataTypeTraits<Value>::zero(),
                                                                         DataTypeTraits<Value>::zero()); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<Value>::MpiDataType(); }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return DataTypeTraits<Value>::axomTypeID(); }
  using AxomType = typename DataTypeTraits<Value>::AxomType;
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<std::tuple<Value, Value, Value, Value> > {
  typedef std::tuple<Value, Value, Value, Value> ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const std::tuple<Value, Value, Value, Value>&) { return 4; }
  static std::tuple<Value, Value, Value, Value> zero() { return std::make_tuple(DataTypeTraits<Value>::zero(),
                                                                                DataTypeTraits<Value>::zero(),
                                                                                DataTypeTraits<Value>::zero(),
                                                                                DataTypeTraits<Value>::zero()); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<Value>::MpiDataType(); }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return DataTypeTraits<Value>::axomTypeID(); }
  using AxomType = typename DataTypeTraits<Value>::AxomType;
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<std::tuple<Value, Value, Value, Value, Value> > {
  typedef std::tuple<Value, Value, Value, Value, Value> ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const std::tuple<Value, Value, Value, Value, Value>&) { return 5; }
  static std::tuple<Value, Value, Value, Value, Value> zero() { return std::make_tuple(DataTypeTraits<Value>::zero(),
                                                                                       DataTypeTraits<Value>::zero(),
                                                                                       DataTypeTraits<Value>::zero(),
                                                                                       DataTypeTraits<Value>::zero(),
                                                                                       DataTypeTraits<Value>::zero()); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<Value>::MpiDataType(); }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return DataTypeTraits<Value>::axomTypeID(); }
  using AxomType = typename DataTypeTraits<Value>::AxomType;
};

//------------------------------------------------------------------------------
template<typename Value1, typename Value2>
struct DataTypeTraits<std::pair<Value1, Value2> > {
  typedef std::pair<Value1, Value2> ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const std::pair<Value1, Value2>&) { return 2; }
  static std::pair<Value1, Value2> zero() { return std::make_pair(DataTypeTraits<Value1>::zero(), DataTypeTraits<Value2>::zero()); }
  static axom::sidre::DataTypeId axomTypeID() { VERIFY2(false, "axom interface not checked for std::pair<T1,T2>"); return DataTypeTraits<Value1>::axomTypeID(); }
  using AxomType = typename DataTypeTraits<Value1>::AxomType;
};

//------------------------------------------------------------------------------
template<typename Value1, typename Value2, typename Hash>
struct DataTypeTraits<std::unordered_map<Value1, Value2, Hash> > {
  typedef std::pair<Value1, Value2> ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const std::unordered_map<Value1, Value2, Hash>& x) { return x.size(); }
  static std::unordered_map<Value1, Value2, Hash> zero() { return std::unordered_map<Value1, Value2, Hash>(); }
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<1>::Vector> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::Vector&) { return Dim<1>::Vector::numElements(); }
  static Dim<1>::Vector zero() { return Dim<1>::Vector::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector1d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<1>::Tensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::Tensor&) { return Dim<1>::Tensor::numElements(); }
  static Dim<1>::Tensor zero() { return Dim<1>::Tensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor1d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<1>::SymTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::SymTensor&) { return Dim<1>::SymTensor::numElements(); }
  static Dim<1>::SymTensor zero() { return Dim<1>::SymTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor1d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<1>::ThirdRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::ThirdRankTensor&) { return Dim<1>::ThirdRankTensor::numElements(); }
  static Dim<1>::ThirdRankTensor zero() { return Dim<1>::ThirdRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor1d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<1>::FourthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::FourthRankTensor&) { return Dim<1>::FourthRankTensor::numElements(); }
  static Dim<1>::FourthRankTensor zero() { return Dim<1>::FourthRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FourthRankTensor1d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<1>::FifthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::FifthRankTensor&) { return Dim<1>::FifthRankTensor::numElements(); }
  static Dim<1>::FifthRankTensor zero() { return Dim<1>::FifthRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FifthRankTensor1d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<1>::FacetedVolume> {
  typedef double ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const Dim<1>::FacetedVolume) { return 0; }
  static Dim<1>::FacetedVolume zero() { return Dim<1>::FacetedVolume(); }

  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<2>::Vector> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::Vector&) { return Dim<2>::Vector::numElements(); }
  static Dim<2>::Vector zero() { return Dim<2>::Vector::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector2d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<2>::Tensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::Tensor&) { return Dim<2>::Tensor::numElements(); }
  static Dim<2>::Tensor zero() { return Dim<2>::Tensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor2d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<2>::SymTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::SymTensor&) { return Dim<2>::SymTensor::numElements(); }
  static Dim<2>::SymTensor zero() { return Dim<2>::SymTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor2d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<2>::ThirdRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::ThirdRankTensor&) { return Dim<2>::ThirdRankTensor::numElements(); }
  static Dim<2>::ThirdRankTensor zero() { return Dim<2>::ThirdRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor2d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<2>::FourthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::FourthRankTensor&) { return Dim<2>::FourthRankTensor::numElements(); }
  static Dim<2>::FourthRankTensor zero() { return Dim<2>::FourthRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FourthRankTensor2d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<2>::FifthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::FifthRankTensor&) { return Dim<2>::FifthRankTensor::numElements(); }
  static Dim<2>::FifthRankTensor zero() { return Dim<2>::FifthRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FifthRankTensor2d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<2>::FacetedVolume> {
  typedef double ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const Dim<2>::FacetedVolume) { return 0; }
  static Dim<2>::FacetedVolume zero() { return Dim<2>::FacetedVolume(); }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<3>::Vector> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::Vector&) { return Dim<3>::Vector::numElements(); }
  static Dim<3>::Vector zero() { return Dim<3>::Vector::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector3d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<3>::Tensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::Tensor&) { return Dim<3>::Tensor::numElements(); }
  static Dim<3>::Tensor zero() { return Dim<3>::Tensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor3d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<3>::SymTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::SymTensor&) { return Dim<3>::SymTensor::numElements(); }
  static Dim<3>::SymTensor zero() { return Dim<3>::SymTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor3d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<3>::ThirdRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::ThirdRankTensor&) { return Dim<3>::ThirdRankTensor::numElements(); }
  static Dim<3>::ThirdRankTensor zero() { return Dim<3>::ThirdRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor3d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<3>::FourthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::FourthRankTensor&) { return Dim<3>::FourthRankTensor::numElements(); }
  static Dim<3>::FourthRankTensor zero() { return Dim<3>::FourthRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FourthRankTensor3d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<3>::FifthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::FifthRankTensor&) { return Dim<3>::FifthRankTensor::numElements(); }
  static Dim<3>::FifthRankTensor zero() { return Dim<3>::FifthRankTensor::zero(); }
#ifdef SPHERAL_ENABLE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FifthRankTensor3d; }
#endif
  static axom::sidre::DataTypeId axomTypeID() { return axom::sidre::DOUBLE_ID; }
  using AxomType = double;
};

template<>
struct DataTypeTraits<Dim<3>::FacetedVolume> {
  typedef double ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const Dim<3>::FacetedVolume) { return 0; }
  static Dim<3>::FacetedVolume zero() { return Dim<3>::FacetedVolume(); }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<PolyClipperVertex2d> {
  typedef PolyClipperVertex2d ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const ElementType& x) { return (DataTypeTraits<Dim<2>::Vector>::numElements(Dim<2>::Vector::zero()) +
                                                         2u + 
                                                         2u +
                                                         x.clips.size()); }
  static ElementType zero() { return PolyClipperVertex2d(); }
  using AxomType = double;
};

template<>
struct DataTypeTraits<PolyClipperVertex3d> {
  typedef PolyClipperVertex3d ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const ElementType& x) { return (DataTypeTraits<Dim<3>::Vector>::numElements(Dim<3>::Vector::zero()) +
                                                         x.neighbors.size() + 
                                                         2u +
                                                         x.clips.size()); }
  static ElementType zero() { return PolyClipperVertex3d(); }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<PolyClipperPlane2d> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const PolyClipperPlane2d&) { return 4; }
  static PolyClipperPlane2d zero() { return PolyClipperPlane2d(); }
  using AxomType = double;
};

template<>
struct DataTypeTraits<PolyClipperPlane3d> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const PolyClipperPlane3d&) { return 5; }
  static PolyClipperPlane3d zero() { return PolyClipperPlane3d(); }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<int ndim>
struct DataTypeTraits<DomainNode<Dim<ndim>>> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const DomainNode<Dim<ndim>>&) { return DomainNode<Dim<ndim>>::packSize(); }
  static DomainNode<Dim<ndim>> zero() { return DomainNode<Dim<ndim>>({0, 0, 0, 0, 0, 0.0, Dim<ndim>::Vector::zero()}); }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<RKOrder> {
  typedef int ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const RKOrder&) { return 1; }
  static RKOrder zero() { return RKOrder::ZerothOrder; }
  using AxomType = int;
};

//------------------------------------------------------------------------------
template<int ndim>
struct DataTypeTraits<RKCoefficients<Dim<ndim>>> {
  typedef double ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const RKCoefficients<Dim<ndim>>& x) { return x.size() + 1; }
  static RKCoefficients<Dim<ndim>> zero() { return RKCoefficients<Dim<ndim>>(); }
  using AxomType = double;
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<uniform_random> {
  typedef char ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const uniform_random&) { return 2*sizeof(size_t) + 2*sizeof(double); }
  static uniform_random zero() { return uniform_random(); }
  using AxomType = char;
};

} // namespace Spheral

#endif
