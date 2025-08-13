#ifndef SPHERAL_FIELD_TEST_TYPES_HH
#define SPHERAL_FIELD_TEST_TYPES_HH

#include "test-basic-datatypes.hh"
#include "test-basic-exec-policies.hh"

using FIELD_TEST_TYPES = Spheral::Test<
    camp::cartesian_product<EXEC_TYPES, TEST_FIELD_DATATYPES>>::Types;

#endif //  SPHERAL_FIELD_TEST_TYPES_HH
