// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/FieldList.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"
#include <Utilities/Logger.hh>

using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldViewDouble = Spheral::FieldView<DIM3, double>;
using FieldListDouble = Spheral::FieldList<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

// Default Testing Size.
static constexpr int N = 100;

class FieldListViewTest : public ::testing::Test {
public:
  GPUCounters fl_count;
  GPUCounters f_count;

  NodeList_t createNodeList(std::string name = "NL", size_t count = N) {
    return NodeList_t(name, count, 0);
  }

  // Increment variables for each action and space
  auto fl_callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
               chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? fl_count.DToHCopies++ : fl_count.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? fl_count.HNumAlloc++ : fl_count.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? fl_count.HNumFree++ : fl_count.DNumFree++;
    // DEBUG_LOG << "fl_callback: " << action << " " << space << " : " << chai::ACTION_MOVE << " " << chai::ACTION_ALLOC << " " << chai::ACTION_FREE << " : " << chai::CPU << " " << chai::GPU;
    // fl_count.print();
    };
  }

  auto f_callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
               chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? f_count.DToHCopies++ : f_count.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? f_count.HNumAlloc++ : f_count.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? f_count.HNumFree++ : f_count.DNumFree++;
    };
  }
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(FieldListViewTypedTest);
template <typename T> class FieldListViewTypedTest : public FieldListViewTest {};

//------------------------------------------------------------------------------
// Capture on device
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, BasicCapture) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, val);
    FieldListDouble field_list;
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());

    field_list.appendField(field1);
    field_list.appendField(field2);

    const size_t numFields = field_list.size();
    auto fl_v = field_list.view();
    fl_v.setCallback(gpu_this->fl_callback());

    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].numElements(), N);
      });
  }

  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = {1, 0, 1, 1, 2, 1};
    f_ref_count  = {2, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

// TODO: Add test for having multiple FL contain the same field

//------------------------------------------------------------------------------
// multi-scope
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, MultiScopeAndTouch) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, val);
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());

    FieldListDouble field_list;
    field_list.appendField(field1);
    field_list.appendField(field2);

    const size_t numFields = field_list.size() ;

    { // Scope 1
    auto fl_v = field_list.view(gpu_this->fl_callback());

    DEBUG_LOG << "Start Kernel 1";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].numElements(), N);
      });
    DEBUG_LOG << "Stop Kernel 1";

    fl_v.touch(chai::CPU, true);
    } // Scope 1

    { // Scope 2
    auto fl_v = field_list.view(gpu_this->fl_callback());

    DEBUG_LOG << "Start Kernel 2";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].numElements(), N);
      });
    DEBUG_LOG << "Stop Kernel 2";
    } // Scope 2
  }

  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = { 2, 0, 2, 2, 3, 2 };
    f_ref_count  = { 4, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

//------------------------------------------------------------------------------
// multi-scope no touch
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, MultiScopeNoTouch) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, val);
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());

    FieldListDouble field_list;
    field_list.appendField(field1);
    field_list.appendField(field2);

    const size_t numFields = field_list.size();

    { // Scope 1
    auto fl_v = field_list.view(gpu_this->fl_callback());

    DEBUG_LOG << "Start Kernel 1";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].numElements(), N);
      });
    DEBUG_LOG << "Stop Kernel 1";
    } // Scope 1

    { // Scope 2
    auto fl_v = field_list.view(gpu_this->fl_callback());

    DEBUG_LOG << "Start Kernel 2";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].numElements(), N);
      });
    DEBUG_LOG << "Stop Kernel 2";
    } // Scope 2
  }

  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = { 2, 0, 2, 2, 3, 2 };
    f_ref_count  = { 2, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

//------------------------------------------------------------------------------
// Test calling FieldListView move function
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, MoveTest) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, 2.*val);
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());
    
    FieldListDouble field_list;
    field_list.appendField(field1);
    field_list.appendField(field2);
    {
      auto fl_v = field_list.view(gpu_this->fl_callback());

      RAJA::forall<TypeParam>(TRS_UINT(0, N),
        [=] SPHERAL_HOST_DEVICE (size_t i) {
          fl_v(0, i) += (double)i;
          fl_v(1, i) += (double)i;
        });
      fl_v.move(chai::CPU);

      for (size_t i = 0; i < N; ++i) {
        SPHERAL_ASSERT_EQ(field_list(0,i), val+(double)i);
        SPHERAL_ASSERT_EQ(field_list(1,i), 2.*val+(double)i);
      }
    }
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = { 1, 1, 1, 1, 2, 1 };
    f_ref_count  = { 2, 2, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

//------------------------------------------------------------------------------
// Test implicit move of Field and FieldList from host RAJA launch
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, HostRajaTest) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, 2.*val);
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());

    FieldListDouble field_list;
    field_list.appendField(field1);
    field_list.appendField(field2);
    {
      auto fl_v = field_list.view(gpu_this->fl_callback());

      RAJA::forall<TypeParam>(TRS_UINT(0, N),
        [=] SPHERAL_HOST_DEVICE (size_t i) {
          fl_v(0, i) += (double)i;
          fl_v(1, i) += (double)i;
        });

      RAJA::forall<RAJA::seq_exec>(TRS_UINT(0, N),
        [=] SPHERAL_HOST_DEVICE (size_t i) {
          SPHERAL_ASSERT_EQ(fl_v(0, i), val+(double)i);
          SPHERAL_ASSERT_EQ(fl_v(1, i), 2.*val+(double)i);
        });
    }
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = { 1, 1, 1, 1, 1, 1 };
    f_ref_count  = { 2, 2, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

//------------------------------------------------------------------------------
// Test implicit copying after host side modification
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, HostResize) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, 2.*val);
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());

    FieldListDouble field_list;
    field_list.appendField(field1);
    {
      auto fl_v = field_list.view(gpu_this->fl_callback());

      EXEC_IN_SPACE_BEGIN(TypeParam)
        SPHERAL_ASSERT_EQ(fl_v.size(), 1);
      EXEC_IN_SPACE_END()

      field_list.appendField(field2);
      fl_v = field_list.view(gpu_this->fl_callback());

      EXEC_IN_SPACE_BEGIN(TypeParam)
        SPHERAL_ASSERT_EQ(fl_v.size(), 2);
      EXEC_IN_SPACE_END()
    }
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = {2, 0, 2, 2, 2, 2};
    f_ref_count  = {2, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

//------------------------------------------------------------------------------
// Test multiple FieldLists holding the same Field
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(FieldListViewTypedTest, MultiFieldList) {

  const double val = 4.;
  NodeList_t nl1 = gpu_this->createNodeList("NL1");
  NodeList_t nl2 = gpu_this->createNodeList("NL2");
  {
    FieldDouble field1("TestField1", nl1, val);
    FieldDouble field2("TestField2", nl2, 2.*val);
    field1.setCallback(gpu_this->f_callback());
    field2.setCallback(gpu_this->f_callback());

    FieldListDouble field_list, field_list2;
    field_list.appendField(field1);
    field_list.appendField(field2);
    field_list2.appendField(field2);
    {
      auto fl_v = field_list.view(gpu_this->fl_callback());
      auto fl_v2 = field_list2.view(gpu_this->fl_callback());
      RAJA::forall<RAJA::seq_exec>(TRS_UINT(0, N),
        [=] (size_t i) {
          fl_v(1, i) += (double)i;
        });

      EXEC_IN_SPACE_BEGIN(TypeParam)
        SPHERAL_ASSERT_EQ(fl_v[1].data(), fl_v2[0].data());
      EXEC_IN_SPACE_END()

      RAJA::forall<TypeParam>(TRS_UINT(0, N),
        [=] SPHERAL_HOST_DEVICE (size_t i) {
          SPHERAL_ASSERT_EQ(fl_v(1, i), fl_v2(0, i));
          fl_v(1, i) *= 2.;
        });

      RAJA::forall<RAJA::seq_exec>(TRS_UINT(0, N),
        [=] SPHERAL_HOST_DEVICE (size_t i) {
          SPHERAL_ASSERT_EQ(fl_v(1, i), 2.*(2.*val+(double)i));
          SPHERAL_ASSERT_EQ(fl_v(1, i), fl_v2(0, i));
        });
    }
  }
  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = {2, 2, 2, 2, 2, 2};
    f_ref_count  = {2, 2, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 0, 0, 0, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

REGISTER_TYPED_TEST_SUITE_P(FieldListViewTypedTest, BasicCapture, MultiScopeAndTouch,
                            MultiScopeNoTouch, MoveTest, HostRajaTest, HostResize,
                            MultiFieldList);

INSTANTIATE_TYPED_TEST_SUITE_P(FieldListView, FieldListViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
