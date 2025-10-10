#ifndef SPHERAL_BASIC_DATATYPES_HH
#define SPHERAL_BASIC_DATATYPES_HH

#include "Utilities/DataTypeTraits.hh"
#include "config.hh"

// clang-format off
using TEST_FIELD_DATATYPES = camp::list<
  //bool
  char
  ,int
  ,size_t
  ,float
  ,double
  ,Spheral::Dim<3>::Vector
  ,Spheral::Dim<3>::Vector3d
  ,Spheral::Dim<3>::Tensor
  ,Spheral::Dim<3>::SymTensor
  ,Spheral::Dim<3>::ThirdRankTensor
  ,Spheral::Dim<3>::FourthRankTensor
  ,Spheral::Dim<3>::FifthRankTensor
>;
// clang-format on

namespace Spheral {
namespace TestUtils {

/**
 * Specializable functions to be used in datatype tests. The default are set up to perform
 * a scalar operation on a given type. If scalar multiplication is not available specializations
 * for the type to be tested should be defined below.
 */
template <typename T> SPHERAL_HOST_DEVICE T initT(T const&) { return T(4); }
template <typename T> SPHERAL_HOST_DEVICE void mutateT(T &val) { val *= 2; }
template <typename T> SPHERAL_HOST_DEVICE T resultT(T const&) { return T(4) * 2; }

} // namespace impl
} // namespace spheral

#endif // SPHERAL_BASIC_DATATYPES_HH
