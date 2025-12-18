#include "field-test-types.hh"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/FieldView.hh"
#include <test-basic-datatypes.hh>

/**
 * This is a type test for using common Spheral datatypes with FieldView in a
 * compile time switchable context. This ensures that Spheral datatypes can be
 * correctly moved to and from the GPU.
 */

using DIM3 = Spheral::Dim<3>;
using NodeList_t = Spheral::NodeList<DIM3>;

class FieldViewDataTypeTest : public ::testing::Test {
public:
  NodeList_t test_node_list = NodeList_t("DataNodeList", 10, 0);
  int someInt = 5;

  void SetUp() override {}
};

// Setting up G Test for Field
TYPED_TEST_SUITE_P(FieldViewDataTypeTypedTest);
template <typename T>
class FieldViewDataTypeTypedTest : public FieldViewDataTypeTest {};

/**
 * This test constructs a host side field object. The field returns a FieldView 
 * to be used within a RAJA context. The FieldView handles data migration between
 * the host and the device.
 *
 * In the working loop we perform some mutation on the underlying data. We then
 * ensure the data is correctly copied back to the host.
 *
 * This test will compile for CPU and GPU execution when available.
 */
GPU_TYPED_TEST_P(FieldViewDataTypeTypedTest, CopyAndMutate) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using DATATYPE = typename camp::at<TypeParam, camp::num<1>>::type;
  {
    using namespace Spheral;
    using FieldType = Field<DIM3, DATATYPE>;
    using FieldViewType = FieldView<DIM3, DATATYPE>;

    std::string field_name = "Field::OffloadCopyAndMutate";
    FieldType field(field_name, gpu_this->test_node_list, Spheral::TestUtils::initT(DATATYPE{}));

    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 10);

    FieldViewType field_v = field.view();
    SPHERAL_ASSERT_EQ(field.size(), 10);

    // clang-format off
    RAJA::forall<WORK_EXEC_POLICY>(
        TRS_UINT(0, field.size()),
        [=] SPHERAL_HOST_DEVICE(int i) {
          SPHERAL_ASSERT_EQ(field_v[i], TestUtils::initT(DATATYPE{}));
          TestUtils::mutateT(field_v[i]);
        }
    );

    RAJA::forall<LOOP_EXEC_POLICY>(
        TRS_UINT(0, field.size()),
        [=] (int i) {
          SPHERAL_ASSERT_EQ(field_v[i], TestUtils::resultT(DATATYPE{}));
        }
    );
    // clang-format on

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

REGISTER_TYPED_TEST_SUITE_P(FieldViewDataTypeTypedTest, CopyAndMutate);

INSTANTIATE_TYPED_TEST_SUITE_P(Field, FieldViewDataTypeTypedTest,
                               FIELD_TEST_TYPES, );
