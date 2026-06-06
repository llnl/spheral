#ifndef SPHERAL_TEST_BASE_HH
#define SPHERAL_TEST_BASE_HH

#include "gtest/gtest.h"
#include "camp/camp.hpp"

namespace Spheral {

template <class T> struct Test;

template <class... T> struct Test<camp::list<T...>> {
  using Types = ::testing::Types<T...>;
};

} // namespace Spheral

#endif //  SPHERAL_TEST_BASE_HH
