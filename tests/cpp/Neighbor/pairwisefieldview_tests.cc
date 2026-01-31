// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "Neighbor/NodePairList.hh"
#include "Neighbor/PairwiseField.hh"

#include <random>
#include <memory>

//------------------------------------------------------------------------------
// These are unit tests for Spheral::PairwiseFieldView with a basic double datatype.
// Spheral::PairwiseFieldView is a host/device capable. It is tested using typed
// tests to check for correct execution on both host and device.
//------------------------------------------------------------------------------
using DIM3 = Spheral::Dim<3>;
using NPIT = Spheral::NodePairIdxType;
using NPLV = Spheral::NodePairListView;
using NPL = Spheral::NodePairList;
using NPLVec = std::vector<Spheral::NodePairIdxType>;
using PairwiseFieldDouble = Spheral::PairwiseField<DIM3, double, 1u>;
using PairwiseFieldViewDouble = Spheral::PairwiseFieldView<double, 1u>;
using PairwiseFieldDoubleDouble = Spheral::PairwiseField<DIM3, double, 2u>;
using PairwiseFieldViewDoubleDouble = Spheral::PairwiseFieldView<double, 2u>;

// Default Testing Size.
static constexpr int N = 10000;

// PairwiseFieldViewTest is constructed at the start of each unit test.
class PairwiseFieldViewTest : public ::testing::Test {
public:
  GPUCounters gcounts;

  // Build a NodePairList
  std::shared_ptr<NPL> createNPL(const size_t count = N) {
    std::random_device rd;
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<size_t> dist(0u, count);
    NPLVec vals(count);
    {
      size_t i = 0u;
      while (i < count) {
        vals[i] = NPIT(dist(gen), dist(gen), dist(gen), dist(gen));
        if ((vals[i].i_list != vals[i].j_list) or
            (vals[i].i_node != vals[i].j_node)) ++i;
      }
    }
    return std::make_shared<NPL>(vals);
  }
  
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

// Setting up Templated Test for PairwiseFieldView
TYPED_TEST_SUITE_P(PairwiseFieldViewTypedTest);
template <typename T> class PairwiseFieldViewTypedTest : public PairwiseFieldViewTest {};


//------------------------------------------------------------------------------
// Host/Device test for the PairwiseFieldView being captured in a RAJA execution space.
// GPU execution spaces should trigger an allocation on the device, a copy from
// the host to the device, and a deallocation when the Field Dtor is triggered.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(PairwiseFieldViewTypedTest, ExecutionSpaceCapture) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    auto pairs = gpu_this->createNPL();
    SPHERAL_ASSERT_EQ(pairs->size(), N);

    PairwiseFieldDouble dfield(pairs);
    PairwiseFieldDoubleDouble ddfield(pairs);
    SPHERAL_ASSERT_EQ(dfield.size(), N);
    SPHERAL_ASSERT_EQ(ddfield.size(), N);

    auto dfield_v = dfield.view();
    auto ddfield_v = ddfield.view();
    SPHERAL_ASSERT_EQ(dfield_v.size(), N);
    SPHERAL_ASSERT_EQ(ddfield_v.size(), N);

    // Fill in known values for the pairwise fields
    for (auto i = 0u; i < N; ++i) {
      dfield_v[i] = 0.0;
      ddfield_v[i][0] = 0.0;
      ddfield_v[i][1] = 1.0;
    }

    // Check the values in a RAJA loop
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, N),
       [=] SPHERAL_HOST_DEVICE (size_t i) {
         SPHERAL_ASSERT_EQ(dfield_v[i], 0.0);
         SPHERAL_ASSERT_EQ(ddfield_v[i][0], 0.0);
         SPHERAL_ASSERT_EQ(ddfield_v[i][1], 1.0);
       });

  } // field and any GPU allocation should be released here.
}

//------------------------------------------------------------------------------
// This test ensures the PairwiseFieldView Data is migrated back and forth between
// RAJA execution spaces through implicit capture.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(PairwiseFieldViewTypedTest, MultiSpaceCapture) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    auto pairs = gpu_this->createNPL();
    SPHERAL_ASSERT_EQ(pairs->size(), N);

    PairwiseFieldDouble dfield(pairs);
    SPHERAL_ASSERT_EQ(dfield.size(), N);

    // Fill in known values for the pairwise fields
    for (auto i = 0u; i < N; ++i) {
      dfield[i] = 10.0 * i;
    }

    auto dfield_v = dfield.view();
    SPHERAL_ASSERT_EQ(dfield_v.size(), N);

    // Execute in working execution space.
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, dfield.size()),
       [=] SPHERAL_HOST_DEVICE (int i) {
         dfield_v[i] *= 2;
       });

    // Execute in a CPU execution space.
    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, dfield.size()),
       [=, &dfield](int i) {
         SPHERAL_ASSERT_EQ(dfield_v[i], i * 20.0);
         SPHERAL_ASSERT_EQ(dfield[i], i * 20.0);
       });

  } // field and any GPU allocation should be released here.
}

//------------------------------------------------------------------------------
// Test the multi-view semantics for a copy. If multiple views are made from a
// single PairwiseField then only one copy should be performed as both views will reference
// the same data.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(PairwiseFieldViewTypedTest, MultiViewSemantics) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    auto pairs = gpu_this->createNPL();
    SPHERAL_ASSERT_EQ(pairs->size(), N);

    PairwiseFieldDouble field(pairs);
    SPHERAL_ASSERT_EQ(field.size(), N);

    // Fill in known values for the pairwise fields
    for (auto i = 0u; i < N; ++i) {
      field[i] = 10.0 * i;
    }

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

    // Capture and execute on all PairwiseFieldView objs in the working space.
    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, field.size()),
       [=] SPHERAL_HOST_DEVICE (int i) {
         SPHERAL_ASSERT_EQ(field_v0[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v1[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v2[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v3[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v4[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v5[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v6[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v7[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v8[i], 10.0 * i);
         SPHERAL_ASSERT_EQ(field_v9[i], 10.0 * i);
       });

  } // field and any GPU allocation should be released here.
}


REGISTER_TYPED_TEST_SUITE_P(PairwiseFieldViewTypedTest,
                            ExecutionSpaceCapture, MultiSpaceCapture, MultiViewSemantics);

INSTANTIATE_TYPED_TEST_SUITE_P(Field, PairwiseFieldViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );

