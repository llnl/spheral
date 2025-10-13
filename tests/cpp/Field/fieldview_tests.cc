// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/FieldView.hh"

/**
 * These are unit tests for Spehral::FieldView with a basic double datatype.
 * Spheral::FieldView is a host/device capable. It is tested using typed
 * tests to check for correct execution on both host and device.
 */

using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldViewDouble = Spheral::FieldView<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

// Default Testing Size.
static constexpr int N = 100;

// FieldViewTest is Constructed at the start of each unit test.
class FieldViewTest : public ::testing::Test {
public:
  GPUCounters gcounts;
  NodeList_t nl = NodeList_t("DataNodeList", N, 0);

  // Increment variables for each action and space
  auto callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
               chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? gcounts.DToHCopies++ : gcounts.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? gcounts.HNumAlloc++ : gcounts.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? gcounts.HNumFree++ : gcounts.DNumFree++;
    };
  }
};

// Setting up Templated Test for FieldView
TYPED_TEST_SUITE_P(FieldViewTypedTest);
template <typename T> class FieldViewTypedTest : public FieldViewTest {};

/**
 * Host/Device test for the FieldView being captured in a RAJA execution space.
 * GPU execution spaces should trigger an allocation on the device, a copy from
 * the host to the device, and a deallocation when the Field Dtor is triggered.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, ExecutionSpaceCapture) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    FieldDouble field("ExecSpaceCapture", gpu_this->nl, 4.0);
    field.setCallback(gpu_this->callback());
    SPHERAL_ASSERT_EQ(field.numElements(), N);

    auto field_v = field.view();
    SPHERAL_ASSERT_EQ(field_v.numElements(), N);

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, field.numElements()),
       [=] SPHERAL_HOST_DEVICE (size_t i) {
         SPHERAL_ASSERT_EQ(field_v[i], 4.0);
       });

  } // field and any GPU allocation should be released here.

  GPUCounters ref_count;
  if (typeid(WORK_EXEC_POLICY) != typeid(SEQ_EXEC_POLICY)) {
    ref_count.HToDCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }

  COMP_COUNTERS(gpu_this->gcounts, ref_count);
}

/**
 * This test ensures the FieldView Data is migrated back and forth between
 * RAJA execution spaces through implicit capture.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, MultiSpaceCapture) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    FieldDouble field("MultiSpaceCapture", gpu_this->nl, 4.0);
    field.setCallback(gpu_this->callback());

    // Setup Field Data w/ iota values.
    std::vector<double> data(N);
    std::iota(data.begin(), data.end(), 0.0);
    field = data;

    auto field_v = field.view();

    // Execute in working execution space.
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, field.numElements()),
       [=] SPHERAL_HOST_DEVICE (int i) {
         field_v[i] *= 2;
       });

    // Execute in a CPU execution space.
    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, field.numElements()),
       [=, &field](int i) {
         SPHERAL_ASSERT_EQ(field_v[i], i * 2);
         SPHERAL_ASSERT_EQ(field[i], i * 2);
       });

  } // field and any GPU allocation should be released here.

  GPUCounters ref_count;
  if (typeid(WORK_EXEC_POLICY) != typeid(SEQ_EXEC_POLICY)) {
    ref_count.HToDCopies = 1;
    ref_count.DToHCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }

  COMP_COUNTERS(gpu_this->gcounts, ref_count);
}

/**
 * Test the multi-view semantics for a copy. If multiple views are made from a
 * single Field then only one copy should be performed as both views will reference
 * the same data.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, MultiViewSemantics) {
  const double val = 4.;
  using WORK_EXEC_POLICY = TypeParam;
  {
    FieldDouble field("MultiViewSemantics", gpu_this->nl, val);
    field.setCallback(gpu_this->callback());
    SPHERAL_ASSERT_EQ(field.numElements(), N);

    // Retreive multiple FieldViews from a Single Field.
    auto field_v0 = field.view();
    auto field_v1 = field.view();
    auto field_v2 = field.view();
    auto field_v3 = field.view();
    auto field_v4 = field.view();
    auto field_v5 = field.view();
    auto field_v6 = field.view();
    auto field_v7 = field.view();
    auto field_v8 = field.view();
    auto field_v9 = field.view();

    // Capture and execute on all FieldView objs in the working space.
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, field.numElements()),
       [=] SPHERAL_HOST_DEVICE (int i) {
         SPHERAL_ASSERT_EQ(field_v0[i], val);
         SPHERAL_ASSERT_EQ(field_v1[i], val);
         SPHERAL_ASSERT_EQ(field_v2[i], val);
         SPHERAL_ASSERT_EQ(field_v3[i], val);
         SPHERAL_ASSERT_EQ(field_v4[i], val);
         SPHERAL_ASSERT_EQ(field_v5[i], val);
         SPHERAL_ASSERT_EQ(field_v6[i], val);
         SPHERAL_ASSERT_EQ(field_v7[i], val);
         SPHERAL_ASSERT_EQ(field_v8[i], val);
         SPHERAL_ASSERT_EQ(field_v9[i], val);
       });

  } // field and any GPU allocation should be released here.

  GPUCounters ref_count;
  if (typeid(WORK_EXEC_POLICY) != typeid(SEQ_EXEC_POLICY)) {
    ref_count.HToDCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }

  COMP_COUNTERS(gpu_this->gcounts, ref_count);
}

/**
 * Resize the field after a copy to the execution space. The Second view
 * Call should trigger a free of any GPU memory and reassign the FieldView
 * CPU pointer to the underlying vectors new address. This test should expect
 * two allocations, two copies to the device and two deallocations on the
 * device.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, ResizeField) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    const double val = 4.0;
    FieldDouble field("ResizeField", gpu_this->nl, val);
    field.setCallback(gpu_this->callback());
    SPHERAL_ASSERT_EQ(field.numElements(), N);

    auto field_v = field.view();
    SPHERAL_ASSERT_EQ(field_v.numElements(), N);

    // Capture the FieldView in the working execution space.
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, field.numElements()),
       [=] SPHERAL_HOST_DEVICE (int i) {
         SPHERAL_ASSERT_EQ(field_v[i], val);
       });

    // We shouldn't have any ghost nodes, but we double check here.
    SPHERAL_ASSERT_EQ(field.numInternalElements(), field.numElements());

    // Resize the NodeList to 10x the original size.
    gpu_this->nl.numInternalNodes(N * 10);

    // Assign field_v again. This should trigger a deallocation of the original
    // GPU data.
    field_v = field.view();
    SPHERAL_ASSERT_EQ(field_v.numElements(), N * 10);

    // Capture field_v in the working execution space again. This should trigger
    // a new copy to the Device if executing on the GPU.
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, field.numElements()),
       [=] SPHERAL_HOST_DEVICE(int i) {
         if (i < N) {SPHERAL_ASSERT_EQ(field_v[i], val);}
         else { SPHERAL_ASSERT_EQ(field_v[i], 0); }
       });

    SPHERAL_ASSERT_EQ(field.numElements(), N * 10);

  } // field and any GPU allocation should be released here.

  GPUCounters ref_count;
  if (typeid(WORK_EXEC_POLICY) != typeid(SEQ_EXEC_POLICY)) {
    ref_count.HToDCopies = 2;
    ref_count.DNumAlloc = 2;
    ref_count.DNumFree = 2;
  }

  COMP_COUNTERS(gpu_this->gcounts, ref_count);
}

REGISTER_TYPED_TEST_SUITE_P(FieldViewTypedTest, ExecutionSpaceCapture, MultiSpaceCapture,
                            MultiViewSemantics, ResizeField);

INSTANTIATE_TYPED_TEST_SUITE_P(Field, FieldViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );

