#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

/**
 * These are unit tests for Spehral::Field with a basic double datatype.
 * Spheral::Field is a host only data structure.
 */

using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

class FieldTest : public ::testing::Test {
public:
  NodeList_t test_node_list = NodeList_t("DataNodeList", 10, 0);
  int someInt = 5;

  void SetUp() override {}
};

/**
 * Basic Host CTor test for the Field Ctor that only takes a name.
 */
TEST_F(FieldTest, NameCtor) {
  {
    std::string field_name = "Field::NameCtor";
    FieldDouble field(field_name);
    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 0);
    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * Basic Host CTor test for the Field Ctor that takes a name and a nodelist.
 * - Uses the FieldTest Nodelist constructed for the testing suite above.
 */
TEST_F(FieldTest, NameNodeListCtor) {
  {
    std::string field_name = "Field::NodeListCtor";
    FieldDouble field(field_name, this->test_node_list);
    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 10);
    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * CTor test for the Field Ctor that takes a name, nodelist and
 * initial value.
 * - Uses the FieldTest Nodelist constructed for the testing suite above.
 */
TEST_F(FieldTest, NameNodeListValCtor) {
  // using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::NodeListValCtor";
    FieldDouble field(field_name, this->test_node_list, 4);
    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 10);

    for (size_t i = 0; i < field.size(); ++i) {
      SPHERAL_ASSERT_EQ(field[i], 4);
    }

    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * Copy CTor test for the Field.
 * - Test w/ double and GeomPolygon.
 */
TEST_F(FieldTest, CopyCtor) {
  // using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::CopyCtor";
    FieldDouble field(field_name, this->test_node_list, 4);

    FieldDouble copy_field(field);

    SPHERAL_ASSERT_EQ(copy_field.name(), field_name);
    SPHERAL_ASSERT_EQ(copy_field.size(), 10);

    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 7);
    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * Assignment operator to a FieldBase pointer.
 */
TEST_F(FieldTest, AssignmentFieldBase) {
  {
    std::string field_name = "Field::AssignmentFieldBase";
    FieldDouble field(field_name, this->test_node_list, 4);

    FieldBase *base = &field;

    SPHERAL_ASSERT_EQ(base->name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), base->size());
    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * Assignment operator of a Field to another Field.
 */
TEST_F(FieldTest, AssignmentField) {
  {
    std::string field_name = "Field::AssignmentField";
    FieldDouble field(field_name, this->test_node_list, 4);

    FieldDouble copy_field("SomeOtherField");
    copy_field = field;

    SPHERAL_ASSERT_EQ(copy_field.size(), 10);

    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);

    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 7);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * Assignment operator of a Field to by a std::vector container.
 */
TEST_F(FieldTest, AssignmentContainerType) {
  {
    std::string field_name = "Field::AssignmentContainer";
    FieldDouble field(field_name, this->test_node_list, 4);

    using ContainerType = std::vector<double>;
    ContainerType data(10);

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 10),
                                   [&](int i) { data[i] = i; });

    field = data;

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 10),
                                   [&](int i) { field.at(i) *= 2; });

    RAJA::forall<LOOP_EXEC_POLICY>(
        TRS_UINT(0, 10), [&](int i) { SPHERAL_ASSERT_EQ(field.at(i), i * 2); });

    SPHERAL_ASSERT_NE(&field[0], &data[0]);

    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

/**
 * Assignment operator of a Field by a single value. Setting all values of
 * the field.
 */
TEST_F(FieldTest, AssignmentDataType) {
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
  {
    std::string field_name = "Field::AssignmentDataType";
    FieldDouble field(field_name, this->test_node_list, 4);

    SPHERAL_ASSERT_EQ(field.size(), 10);

    RAJA::forall<LOOP_EXEC_POLICY>(
        TRS_UINT(0, 10), [=](int i) { SPHERAL_ASSERT_EQ(field.at(i), 4); });

    field = double(3);

    RAJA::forall<LOOP_EXEC_POLICY>(
        TRS_UINT(0, 10), [=](int i) { SPHERAL_ASSERT_EQ(field.at(i), 3); });
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}
