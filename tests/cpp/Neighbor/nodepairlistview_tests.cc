// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"
#include "Neighbor/NodePairList.hh"

#include <typeinfo>

using NPIT = Spheral::NodePairIdxType;
using NPLV = Spheral::NodePairListView;
using NPL = Spheral::NodePairList;
using NPLVec = std::vector<Spheral::NodePairIdxType>;

static constexpr size_t N = 5;

class NPLVTest : public ::testing::Test {
public:
  GPUCounters n_count;
  // Helper to create a ContainerType with values [start, start+count)
  NPLVec createVec(size_t count = N) {
    NPLVec vals;
    for (size_t i = 0; i < count; ++i) {
      NPIT nit(i, i+1, 2*i, 2*i+1, (double)i);
      vals.push_back(nit);
    }
    return vals;
  }
  NPL createContainer(size_t count = N) {
    NPL npl(createVec(N));
    return npl;
  }
  // Increment variables for each action and space
  auto callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
               chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? n_count.DToHCopies++ : n_count.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? n_count.HNumAlloc++ : n_count.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? n_count.HNumFree++ : n_count.DNumFree++;
    };
  }
};

// Setting up Typed Test Suite for NodePairListView
TYPED_TEST_SUITE_P(NPLViewTypedTest);
template <typename T> class NPLViewTypedTest : public NPLVTest {};

// Test default constructor
GPU_TYPED_TEST_P(NPLViewTypedTest, DefaultConstructor) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    NPLV npl_v;
    SPHERAL_ASSERT_EQ(npl_v.size(), 0);
  EXEC_IN_SPACE_END()
}

// Test copy and assignment
GPU_TYPED_TEST_P(NPLViewTypedTest, CopyAssign) {
  NPL npl = gpu_this->createContainer();
  SPHERAL_ASSERT_EQ(npl.size(), N);
  npl.setUserCallback(gpu_this->callback());
  {
    NPL npl_2(npl);
    NPL npl_3 = npl;
    // TODO: Make callbacks copied
    npl_2.setUserCallback(gpu_this->callback());
    npl_3.setUserCallback(gpu_this->callback());
    SPHERAL_ASSERT_NE(npl.data(), npl_2.data());
    SPHERAL_ASSERT_EQ(npl_2.size(), npl.size());
    SPHERAL_ASSERT_EQ(npl_3.size(), npl.size());
    NPLV npl2_view = npl_2.view();
    NPLV npl3_view = npl_3.view();
    NPLV npl_v = npl.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_NE(npl2_view.data(), npl_v.data());
      SPHERAL_ASSERT_NE(npl3_view.data(), npl_v.data());
    EXEC_IN_SPACE_END()
    RAJA::forall<TypeParam>(TRS_UINT(0, N),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
        SPHERAL_ASSERT_EQ(npl2_view[i], npl_v[i]);
        SPHERAL_ASSERT_EQ(npl3_view[i], npl_v[i]);
      });
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    // npl is not destroyed yet so only 2 frees on device
    ref_count = {3, 0, 0, 3, 0, 2};
  }
  COMP_COUNTERS(gpu_this->n_count, ref_count);
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, ConstructorFromContainer) {
  {
    NPL npl = gpu_this->createContainer();
    npl.setUserCallback(gpu_this->callback());
    NPLV npl_v = npl.view();
    SPHERAL_ASSERT_EQ(npl_v.size(), N);
    SPHERAL_ASSERT_EQ(npl_v.data(), npl.data());

    RAJA::forall<TypeParam>(TRS_UINT(0, N),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
        NPIT nit_ref(i, i+1, 2*i, 2*i+1, (double)i);
        SPHERAL_ASSERT_EQ(npl_v[i], nit_ref);
        npl_v[i].i_node *= 2;
        npl_v[i].i_list -= 1;
        npl_v[i].f_couple *= 2.;
      });

    npl_v.move(chai::CPU);

    for (size_t i = 0; i < N; ++i) {
      NPIT nit_ref(2*i, i, 2*i, 2*i+1, 2*(double)i);
      SPHERAL_ASSERT_EQ(npl[i], nit_ref);
    }
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count = {1, 1, 0, 1, 0, 1};
  }
  COMP_COUNTERS(gpu_this->n_count, ref_count);
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, Touch) {
  {
    NPL npl = gpu_this->createContainer();
    npl.setUserCallback(gpu_this->callback());
    NPLV npl_v = npl.view();

    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(npl_v.size(), N);
    EXEC_IN_SPACE_END()
    npl[0].i_list = 4; // Modify the data on the host
    npl_v.touch(chai::CPU); // Change the execution space for the MA
    npl_v = npl.view(); // Create a new view

    RAJA::forall<TypeParam>(TRS_UINT(0, N),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
         if (i == 0) {
           SPHERAL_ASSERT_EQ(npl_v[i].i_list, 4);
         }
       });
    npl_v.touch(chai::CPU);
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count = {2, 0, 0, 1, 0, 1};
  }
  COMP_COUNTERS(gpu_this->n_count, ref_count);
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, Resize) {
  {
    NPLVec npl_vec = gpu_this->createVec();
    NPL npl(npl_vec);
    npl.setUserCallback(gpu_this->callback());
    NPLV npl_v = npl.view();

    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(npl_v.size(), N);
    EXEC_IN_SPACE_END()

    NPIT nit(4, 4, 4, 4, 4.);
    npl_vec.push_back(nit);
    npl.fill(npl_vec);
    // TODO: Make it so callbacks persist
    npl.setUserCallback(gpu_this->callback());
    npl_v = npl.view();

    RAJA::forall<TypeParam>(TRS_UINT(0, npl.size()),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
        SPHERAL_ASSERT_EQ(npl_v.size(), N+1);
        if (i == N) {
          SPHERAL_ASSERT_EQ(npl_v[i].i_list, 4);
          npl_v[i].i_node = 6;
        }
      });
    npl_v.move(chai::CPU);
    SPHERAL_ASSERT_EQ(npl_v[N].i_node, npl[N].i_node);
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count = {2, 1, 0, 2, 0, 2};
  }
  COMP_COUNTERS(gpu_this->n_count, ref_count);
}

GPU_TYPED_TEST_P(NPLViewTypedTest, ScopeChanges) {
  {
    NPLVec npl_vec = gpu_this->createVec();
    NPL npl(std::move(npl_vec));
    npl.setUserCallback(gpu_this->callback());
    {
      NPLV npl2_v = npl.view();
      SPHERAL_ASSERT_EQ(npl2_v.size(), npl.size());
      RAJA::forall<TypeParam>(TRS_UINT(0, npl.size()),
        [=] SPHERAL_HOST_DEVICE(size_t i) {
          SPHERAL_ASSERT_EQ(npl2_v[i].i_node, i);
        });
    }
    NPLV npl3_v = npl.view();
    RAJA::forall<TypeParam>(TRS_UINT(0, npl.size()),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
        SPHERAL_ASSERT_EQ(npl3_v[i].i_node, i);
      });
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count = {1, 0, 0, 1, 0, 1};
  }
  COMP_COUNTERS(gpu_this->n_count, ref_count);
}

GPU_TYPED_TEST_P(NPLViewTypedTest, EmptyInit) {
  NPLVec npl_vec = gpu_this->createVec();
  NPL npl(npl_vec);
  NPLV npl_v = npl.view();
  SPHERAL_ASSERT_EQ(npl_v.size(), N);
  SPHERAL_ASSERT_EQ(npl.size(), N);
  npl_vec = gpu_this->createVec();
  npl_vec.clear();
  npl.fill(npl_vec);
  npl_v = npl.view();
  SPHERAL_ASSERT_EQ(npl_v.size(), 0u);
  SPHERAL_ASSERT_EQ(npl.size(), 0u);
}

GPU_TYPED_TEST_P(NPLViewTypedTest, ViewVector) {
  NPLVec npl_vec = gpu_this->createVec();
  NPL npl(std::move(npl_vec));
  npl_vec = gpu_this->createVec();
  NPL npl2(std::move(npl_vec));
  {
    std::vector<NPLV> view_vec;
    for (size_t i = 0; i < 4; ++i) {
      view_vec.push_back(npl.view());
    }
    NPLV& ref_view = view_vec[3];
    SPHERAL_ASSERT_EQ(ref_view[0].i_node, 0u);
    view_vec.clear();
    for (size_t i = 0; i < 4; ++i) {
      view_vec.push_back(npl.view());
    }
    NPLV& ref_view2 = view_vec[3];
    SPHERAL_ASSERT_EQ(ref_view2[0].i_node, 0u);
    NPLV later_view = npl.view();
    SPHERAL_ASSERT_EQ(later_view[0].i_node, 0u);
  }
  SPHERAL_ASSERT_EQ(npl.size(), N);
}

REGISTER_TYPED_TEST_SUITE_P(NPLViewTypedTest, DefaultConstructor, CopyAssign,
                            ConstructorFromContainer, Touch, Resize, ScopeChanges,
                            EmptyInit, ViewVector);

INSTANTIATE_TYPED_TEST_SUITE_P(NodePairListView, NPLViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
