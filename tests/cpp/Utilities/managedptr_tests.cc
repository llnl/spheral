// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"
#include "chai/managed_ptr.hpp"

#include <Utilities/Logger.hh>
#include "Geometry/GeomVector.hh"

// Base class with a virtual function
class Base {
public:
  SPHERAL_HOST_DEVICE virtual double operate(double a_val) const = 0;
  SPHERAL_HOST_DEVICE virtual ~Base() = default;
  SPHERAL_HOST_DEVICE Base() = default;
  SPHERAL_HOST_DEVICE Base(double a_val1, double a_val2) :
    m_val1(a_val1),
    m_val2(a_val2)
  {}
  SPHERAL_HOST_DEVICE double getVal1() { return m_val1; }
  SPHERAL_HOST_DEVICE double getVal2() { return m_val2; }
  SPHERAL_HOST_DEVICE void setVal1(const double a_val) { m_val1 = a_val; }
protected:
  double m_val1 = 0.;
  double m_val2 = 0.;
};

// Derived class 1
class DerivedA final : public Base {
public:
  SPHERAL_HOST_DEVICE DerivedA(double a_val1, double a_val2, double a_val3) :
    Base(a_val1, a_val2), m_val3(a_val3) {}
  SPHERAL_HOST_DEVICE virtual double operate(double a_val) const override {
    return a_val*(m_val1 + m_val2 + m_val3);
  }
protected:
  double m_val3 = 10.;
  using Base::m_val1;
  using Base::m_val2;
};

// Derived class 2
class DerivedB final : public Base {
public:
  SPHERAL_HOST_DEVICE DerivedB(double a_val1, double a_val2) :
    Base(a_val1, a_val2) {}
  SPHERAL_HOST_DEVICE virtual double operate(double a_val) const override {
    return -a_val*(m_val1 + m_val2);
  }
protected:
  using Base::m_val1;
  using Base::m_val2;
};

class ManagedPointerTest : public ::testing::Test {
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(ManagedPointerTypedTest);
template <typename T> class ManagedPointerTypedTest : public ManagedPointerTest {};

GPU_TYPED_TEST_P(ManagedPointerTypedTest, Start) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    Spheral::GeomVector<3> vec(1., 2., 3.);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(ManagedPointerTypedTest, BasicCapture) {
  double val0 = 2.;
  double val1 = 1.;
  double val2 = 3.;
  double val3 = 4.;
  double ref_valA = val0*(val1 + val2 + val3);
  double ref_valB = -val0*(val1 + val2);
  for (int i = 0; i < 10; ++i) {
    double ref_val = ref_valA;
    chai::managed_ptr<Base> d_ptr;
    if (i%2 == 0) {
      d_ptr = chai::make_managed<DerivedA>(val1, val2, val3);
    } else {
      d_ptr = chai::make_managed<DerivedB>(val1, val2);
      ref_val = ref_valB;
    }
    EXEC_IN_SPACE_BEGIN(TypeParam)
      double val = d_ptr->operate(val0);
      SPHERAL_ASSERT_FLOAT_EQ(val, ref_val);
    EXEC_IN_SPACE_END()
    d_ptr.free();
  }
}

GPU_TYPED_TEST_P(ManagedPointerTypedTest, ModifyClass) {
  double val0 = 2.;
  double val1 = 1.;
  double val2 = 3.;
  double val3 = 4.;
  double val12 = 10.;
  double ref_valA = val0*(val1 + val2 + val3);
  double ref_valA2 = val0*(val12 + val2 + val3);
  chai::managed_ptr<Base> d_ptr = chai::make_managed<DerivedA>(val1, val2, val3);
  for (int i = 0; i < 10; ++i) {
    double ref_val = ref_valA;
    if (i%2 == 0) {
      EXEC_IN_SPACE_BEGIN(TypeParam)
        d_ptr->setVal1(val1);
      EXEC_IN_SPACE_END()
    } else {
      EXEC_IN_SPACE_BEGIN(TypeParam)
        d_ptr->setVal1(val12);
      EXEC_IN_SPACE_END()
      ref_val = ref_valA2;
    }
    EXEC_IN_SPACE_BEGIN(TypeParam)
      double val = d_ptr->operate(val0);
      SPHERAL_ASSERT_FLOAT_EQ(val, ref_val);
    EXEC_IN_SPACE_END()
  }
  d_ptr.free();
}

REGISTER_TYPED_TEST_SUITE_P(ManagedPointerTypedTest, Start, BasicCapture, ModifyClass);

INSTANTIATE_TYPED_TEST_SUITE_P(ManagedPointer, ManagedPointerTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
