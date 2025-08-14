#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

class OffloadTest : public ::testing::Test {};

// Setting up G Test for Field
TYPED_TEST_SUITE_P(OffloadTypedTest);
template <typename T> class OffloadTypedTest : public OffloadTest {};

class TestType {
public:
  int data[3];
  SPHERAL_HOST_DEVICE
  TestType() : data{0, 0, 0} {}

  SPHERAL_HOST_DEVICE
  TestType(int d) : data{d, d, d} {}

  SPHERAL_HOST_DEVICE
  void operator+=(TestType const &rhs) {
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
  }

  SPHERAL_HOST_DEVICE
  bool operator==(int i) {
    return (data[0] == i && data[1] == i && data[2] == i);
  }
};

#define TEST_SIZE 50

/**
 * Testing execution of a RAJA kernel for Host & any compiled offload
 * platforms.
 */
GPU_TYPED_TEST_P(OffloadTypedTest, RajaLoop) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using WORK_RESOURCE = typename camp::at<TypeParam, camp::num<1>>::type;

  using T = TestType;

  std::vector<T> h_data(TEST_SIZE);
  for (size_t i = 0; i < TEST_SIZE; i++) {
    h_data[i] = T(i);
  }
  auto resource = WORK_RESOURCE::get_default();
  T *o_data = resource.template allocate<T>(TEST_SIZE);

  resource.memcpy(o_data, h_data.data(), TEST_SIZE * sizeof(T));

  RAJA::forall<WORK_EXEC_POLICY>(
      resource, TRS_UINT(0, TEST_SIZE),
      [=] SPHERAL_HOST_DEVICE(int i) { o_data[i] += T(i + 1); });

  resource.memcpy(h_data.data(), o_data, TEST_SIZE * sizeof(T));

  for (int i = 0; i < TEST_SIZE; ++i) {
    ASSERT_TRUE(h_data[i] == i + i + 1);
  }

  resource.deallocate(o_data);
}

REGISTER_TYPED_TEST_SUITE_P(OffloadTypedTest, RajaLoop);

INSTANTIATE_TYPED_TEST_SUITE_P(Offload, OffloadTypedTest,
                               EXEC_RESOURCE_TYPES, );
