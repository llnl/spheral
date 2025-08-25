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

static constexpr size_t N = 5;

class NodePairListViewTest : public ::testing::Test {
public:
  GPUCounters n_count;
  // Helper to create a ContainerType with values [start, start+count)
  NPL createContainer(size_t count = N) {
    NPL npl ;
    for (size_t i = 0; i < count; ++i) {
      NPIT nit(i, i+1, 2*i, 2*i+1, (double)i);
      npl.push_back(nit);
    }
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

// Setting up Typed Test Suite for GeomVector
TYPED_TEST_SUITE_P(NPLViewTypedTest);
template <typename T> class NPLViewTypedTest : public NodePairListViewTest {};

// Test default constructor
GPU_TYPED_TEST_P(NPLViewTypedTest, DefaultConstructor) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    NPLV npl_v;
    SPHERAL_ASSERT_EQ(npl_v.size(), 0);
  EXEC_IN_SPACE_END()
}

// Test default constructor
GPU_TYPED_TEST_P(NPLViewTypedTest, CopyAssign) {
  NPL npl = gpu_this->createContainer();
  {
    NPL npl_2(npl);
    NPL npl_3 = npl;
    SPHERAL_ASSERT_EQ(npl_2.size(), npl.size());
    SPHERAL_ASSERT_EQ(npl_3.size(), npl.size());
    NPLV npl2_view = npl_2.view();
    NPLV npl3_view = npl_3.view();
    NPLV npl_v = npl.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_NE(npl_2.data(), npl_v.data());
      SPHERAL_ASSERT_NE(npl_2.data(), npl_v.data());
    EXEC_IN_SPACE_END()
    RAJA::forall<TypeParam>(TRS_UINT(0, N),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
        SPHERAL_ASSERT_EQ(npl2_view[i], npl_v[i]);
        SPHERAL_ASSERT_EQ(npl3_view[i], npl_v[i]);
      });
  }
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, ConstructorFromContainer) {
  NPL npl = gpu_this->createContainer();
  NPLV npl_v = npl.view(gpu_this->callback());
  SPHERAL_ASSERT_EQ(npl_v.size(), N);
  SPHERAL_ASSERT_EQ(&(npl_v[0]), &(npl[0]));

  RAJA::forall<TypeParam>(TRS_UINT(0, N),
    [=] SPHERAL_HOST_DEVICE(size_t i) {
      NPIT nit = npl_v[i];
      SPHERAL_ASSERT_EQ(nit.i_node, i);
      SPHERAL_ASSERT_EQ(nit.i_list, i+1);
      SPHERAL_ASSERT_EQ(nit.j_node, 2*i);
      SPHERAL_ASSERT_EQ(nit.j_list, 2*i+1);
      SPHERAL_ASSERT_EQ(nit.f_couple, (double)i);
      npl_v[i].i_node *= 2;
      npl_v[i].i_list -= 1;
      npl_v[i].f_couple *= 2.;
    });

  npl_v.move(chai::CPU);

  for (size_t i = 0; i < N; ++i) {
    SPHERAL_ASSERT_EQ(npl[i].i_node, 2*i);
    SPHERAL_ASSERT_EQ(npl[i].i_list, i);
    SPHERAL_ASSERT_EQ(npl[i].f_couple, 2.*(double)i);
  }
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, Touch) {
  {
    NPL npl = gpu_this->createContainer();
    NPLV npl_v = npl.view(gpu_this->callback());

    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(npl_v.size(), N);
    EXEC_IN_SPACE_END()

    npl[0].i_list = 4; // Modify the data on the host
    npl_v.touch(chai::CPU); // Change the execution space for the MA
    npl_v = npl.view(gpu_this->callback()); // Create a new view

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
    NPL npl = gpu_this->createContainer();
    NPLV npl_v = npl.view(gpu_this->callback());

    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(npl_v.size(), N);
    EXEC_IN_SPACE_END()

    NPIT nit(4, 4, 4, 4, 4.);
    npl.push_back(nit);
    npl_v = npl.view(gpu_this->callback());

    RAJA::forall<TypeParam>(TRS_UINT(0, npl.size()),
      [=] SPHERAL_HOST_DEVICE(size_t i) {
        SPHERAL_ASSERT_EQ(npl_v.size(), N+1);
        if (i == N) {
          SPHERAL_ASSERT_EQ(npl_v[i].i_list, 4);
          npl_v[i].i_node = 6;
        }
      });
    npl_v.move(chai::CPU);
    SPHERAL_ASSERT_EQ(npl_v[N+1].i_node, npl[N+1].i_node);
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count = {2, 1, 0, 2, 0, 2};
  }
  COMP_COUNTERS(gpu_this->n_count, ref_count);
}

REGISTER_TYPED_TEST_SUITE_P(NPLViewTypedTest, DefaultConstructor, CopyAssign,
                            ConstructorFromContainer, Touch, Resize);

INSTANTIATE_TYPED_TEST_SUITE_P(NodePairListView, NPLViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
