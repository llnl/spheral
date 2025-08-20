// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Utilities/FlexPointer.hh"
#include <Utilities/Logger.hh>
#include "Geometry/Geom3Vector.hh"

// Base class with a virtual function
class Base {
public:
  SPHERAL_HOST_DEVICE virtual double operate(double a_val) const = 0;
  SPHERAL_HOST_DEVICE virtual ~Base() = default;
  SPHERAL_HOST_DEVICE Base() = default;
  SPHERAL_HOST_DEVICE Base(const double a_val1, const double a_val2) :
    m_val1(a_val1),
    m_val2(a_val2)
  {}
  SPHERAL_HOST_DEVICE double getVal1() { return m_val1; }
  SPHERAL_HOST_DEVICE double getVal2() { return m_val2; }
protected:
  double m_val1 = 0.;
  double m_val2 = 0.;
};

// Derived class 1
class DerivedA final : public Base {
public:
  SPHERAL_HOST_DEVICE DerivedA(const double a_val1, const double a_val2, const double a_val3) :
    Base(a_val1, a_val2), m_val3(a_val3) {}
  SPHERAL_HOST_DEVICE double operate(double a_val) const override {
    return a_val*(m_val1 + m_val2 + m_val3);
  }
protected:
  double m_val3 = 10.;
};

// Derived class 2
class DerivedB final : public Base {
public:
  SPHERAL_HOST_DEVICE DerivedB(const double a_val1, const double a_val2) :
    Base(a_val1, a_val2) {}
  SPHERAL_HOST_DEVICE double operate(double a_val) const override {
    return -a_val*(m_val1 + m_val2);
  }
};

class FlexPointerTest : public ::testing::Test {
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(FlexPointerTypedTest);
template <typename T> class FlexPointerTypedTest : public FlexPointerTest {};

GPU_TYPED_TEST_P(FlexPointerTypedTest, Start) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    Spheral::GeomVector<3> vec(1., 2., 3.);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(FlexPointerTypedTest, BasicCapture) {
  const double val0 = 2.;
  const double val1 = 1.;
  const double val2 = 3.;
  const double val3 = 4.;
  const double ref_valA = val0*(val1 + val2 + val3);
  const double ref_valB = -val0*(val1 + val2);
  chai::ExecutionSpace space = chai::GPU;
  if (typeid(RAJA::seq_exec) == typeid(TypeParam)) {
    space = chai::CPU;
  }
  for (int i = 0; i < 10; ++i) {
    Spheral::FlexPointer<Base> base_ptr;
    double ref_val = ref_valA;
    if (i%2 == 0) {
      base_ptr.initialize<DerivedA>(space, val1, val2, val3);
    } else {
      base_ptr.initialize<DerivedB>(space, val1, val2);
      ref_val = ref_valB;
    }
    Base* d_ptr = base_ptr.getPointer();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      double val = d_ptr->operate(val0);
      SPHERAL_ASSERT_FLOAT_EQ(val, ref_val);
    EXEC_IN_SPACE_END()
  }
}

REGISTER_TYPED_TEST_SUITE_P(FlexPointerTypedTest, Start, BasicCapture);

INSTANTIATE_TYPED_TEST_SUITE_P(FlexPointer, FlexPointerTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
