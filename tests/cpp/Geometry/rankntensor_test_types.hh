#ifndef SPHERAL_RANKNTENSOR_TEST_TYPES_HH
#define SPHERAL_RANKNTENSOR_TEST_TYPES_HH

#include "test-base.hh"
#include "test-basic-exec-policies.hh"

#include "Geometry/GeomFifthRankTensor.hh"
#include "Geometry/GeomFourthRankTensor.hh"
#include "Geometry/GeomThirdRankTensor.hh"

/*
 * Instantiate the typed test suite for different Rank-N tensor types.
 */
using RANK_N_TENSOR_TYPES = camp::list<
  Spheral::GeomThirdRankTensor<3>
  ,Spheral::GeomFourthRankTensor<3>
  ,Spheral::GeomFifthRankTensor<3>
>;

using RANK_N_TENSOR_TEST_TYPES = Spheral::Test<
  camp::cartesian_product<EXEC_TYPES, RANK_N_TENSOR_TYPES>>::Types;

#endif //  SPHERAL_RANKNTENSOR_TEST_TYPES_HH
