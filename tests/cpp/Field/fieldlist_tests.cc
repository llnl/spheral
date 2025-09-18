#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

#include <random>
#include <algorithm>
#include <limits>

// #ifdef USE_MPI
// //------------------------------------------------------------------------------
// // Create our own Test framework main to initialize/finalize MPI
// //------------------------------------------------------------------------------
// #include <mpi.h>
// class MPIEnvironment : public testing::Environment {
// public:
//   void SetUp() override {
//     // MPI_Init is called in main, so no need to call it here.
//   }
//   void TearDown() override {
//     MPI_Finalize();
//   }
// };

// int main(int argc, char* argv[]) {
//   printf("******** MY MAIN ************\n");

//   // Add the MPI environment to handle MPI_Finalize
//   // ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
//   ::testing::InitGoogleTest(&argc, argv);

//   MPI_Init(&argc, &argv);
//   auto comm = MPI_COMM_WORLD;
//   Spheral::Communicator::communicator(comm);
//   int rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   printf("Straight rank call : %d\n", rank);
//   printf("Process::getRank() : %d\n", Spheral::Process::getRank());

//   int result = RUN_ALL_TESTS();

//   MPI_Finalize();
//   // MPI_Finalize is handled by MPIEnvironment::TearDown(),
//   // so no explicit call here.

//   printf("******** MY MAIN ************\n");
//   return result;
// }
// #endif

//------------------------------------------------------------------------------
// These are unit tests for Spheral::FieldList with a basic double datatype.
// Spheral::FieldList is a host only data structure.
//------------------------------------------------------------------------------
using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldListDouble = Spheral::FieldList<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

class FieldListTest : public ::testing::Test {
public:
  NodeList_t nodes1 = NodeList_t("nodes1", 10, 5);
  NodeList_t nodes2 = NodeList_t("nodes2", 20, 10);
  NodeList_t nodes3 = NodeList_t("nodes3", 40, 20);

  void SetUp() override {}
};

//------------------------------------------------------------------------------
// Helper to check that two FieldLists are equal
//------------------------------------------------------------------------------
void checkEquivalence(const FieldListDouble& fl1,
                      const FieldListDouble& fl2) {
  const auto nf = fl1.numFields();
  const auto st = fl1.storageType();
  SPHERAL_ASSERT_EQ(fl2.numFields(), nf);
  SPHERAL_ASSERT_EQ(fl1.storageType(), fl2.storageType());
  
  // Depending on storage type, the Fields of 2 may be references to those of 1
  if (st == Spheral::FieldStorageType::ReferenceFields) {
    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, nf),
                                   [&](size_t k) { SPHERAL_ASSERT_TRUE(fl2.haveField(*fl1[k])); });
  } else {
    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, nf),
                                   [&](size_t k) { SPHERAL_ASSERT_FALSE(fl2.haveField(*fl1[k])); });
  }

  // Check Fields are equal
  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, nf),
                                 [&](size_t k) { SPHERAL_ASSERT_EQ(*fl1[k], *fl2[k]); });
}

//------------------------------------------------------------------------------
// Fill the Field values with random numbers
//------------------------------------------------------------------------------
void fillRandom(FieldListDouble& fl,
                double minVal = -1e4,
                double maxVal =  1e4) {
  std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> realDistribution(minVal, maxVal);
  const auto nf = fl.numFields();
  for (auto k = 0u; k < nf; ++k) {
    const auto n = fl[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) fl(k,i) = realDistribution(generator);
    // RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
    //                                [&](size_t i) { fl(k,i) = realDistribution(generator); });
  }
}

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
TEST_F(FieldListTest, DefaultCtor) {
  FieldListDouble fl;
  SPHERAL_ASSERT_EQ(fl.numFields(), 0u);
  SPHERAL_ASSERT_EQ(fl.storageType(), Spheral::FieldStorageType::ReferenceFields);
}

TEST_F(FieldListTest, StorageTypeCtor) {
  FieldListDouble fl1(Spheral::FieldStorageType::ReferenceFields);
  FieldListDouble fl2(Spheral::FieldStorageType::CopyFields);
  SPHERAL_ASSERT_EQ(fl1.storageType(), Spheral::FieldStorageType::ReferenceFields);
  SPHERAL_ASSERT_EQ(fl2.storageType(), Spheral::FieldStorageType::CopyFields);
}

TEST_F(FieldListTest, CopyCtor) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("field1", nodes1, 0.0);
  fl1.appendNewField("field2", nodes2, 1.0);
  FieldListDouble fl2(fl1);
  checkEquivalence(fl1, fl2);
}

//------------------------------------------------------------------------------
// Assignment with a FieldList
//------------------------------------------------------------------------------
TEST_F(FieldListTest, AssignmentFL) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("field1", nodes1, 0.0);
  fl1.appendNewField("field2", nodes2, 1.0);

  // Assignment with a new FieldList on construction
  {
    FieldListDouble fl2 = fl1;
    checkEquivalence(fl1, fl2);
  }

  // Assignment of an existing FieldList (reference)
  {
    FieldListDouble fl2(Spheral::FieldStorageType::ReferenceFields);
    fl2 = fl1;
    checkEquivalence(fl1, fl2);
  }

  // Assignment of an existing FieldList (copy)
  {
    FieldListDouble fl2(Spheral::FieldStorageType::CopyFields);
    fl2 = fl1;
    checkEquivalence(fl1, fl2);
  }

}

//------------------------------------------------------------------------------
// Assignment with a Scalar
//------------------------------------------------------------------------------
TEST_F(FieldListTest, AssignmentScalar) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("field1", nodes1, 0.0);
  fl1.appendNewField("field2", nodes2, 1.0);
  fl1 = 14.5;
  for (auto* fptr: fl1) {
    SPHERAL_ASSERT_EQ(*fptr, 14.5);
  }
}

//------------------------------------------------------------------------------
// copyFields
//------------------------------------------------------------------------------
TEST_F(FieldListTest, CopyFields) {
  FieldListDouble fl(Spheral::FieldStorageType::ReferenceFields);
  fl.appendField(nodes1.mass());
  fl.appendField(nodes2.mass());
  SPHERAL_ASSERT_EQ(fl[0], &nodes1.mass());
  SPHERAL_ASSERT_EQ(fl[1], &nodes2.mass());
  SPHERAL_ASSERT_EQ(fl, 0.0);
  fl.copyFields();
  fl = 1.0;
  SPHERAL_ASSERT_NE(fl[0], &nodes1.mass());
  SPHERAL_ASSERT_NE(fl[1], &nodes2.mass());
  for (auto* fptr: fl) {
    SPHERAL_ASSERT_EQ(*fptr, 1.0);
  }
 }
 
//------------------------------------------------------------------------------
// haveField/haveNodeList
//------------------------------------------------------------------------------
TEST_F(FieldListTest, haveFieldAndNodeList) {
  FieldListDouble fl(Spheral::FieldStorageType::ReferenceFields);
  fl.appendField(nodes1.mass());
  SPHERAL_ASSERT_EQ(fl[0], &nodes1.mass());
  SPHERAL_ASSERT_TRUE(fl.haveField(nodes1.mass()));
  SPHERAL_ASSERT_TRUE(fl.haveNodeList(nodes1));
  SPHERAL_ASSERT_FALSE(fl.haveField(nodes2.mass()));
  SPHERAL_ASSERT_FALSE(fl.haveNodeList(nodes2));
}

//------------------------------------------------------------------------------
// assignFields
//------------------------------------------------------------------------------
TEST_F(FieldListTest, assignFields) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  FieldListDouble fl2(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff1", nodes1, 1.0);
  fl1.appendNewField("stuff1", nodes2, -3.0);
  fl2.appendNewField("stuff2", nodes1, 0.0);
  fl2.appendNewField("stuff2", nodes2, 0.0);
  fl2.assignFields(fl1);
  for (auto k = 0u; k < 2u; ++k) {
    SPHERAL_ASSERT_EQ(fl2[k]->name(), "stuff2");
    SPHERAL_ASSERT_EQ(*fl2[k], *fl1[k]);
  }
}

//------------------------------------------------------------------------------
// referenceFields
//------------------------------------------------------------------------------
TEST_F(FieldListTest, referenceFields) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  FieldListDouble fl2(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff1", nodes1, 1.0);
  fl1.appendNewField("stuff1", nodes2, -3.0);
  fl2.referenceFields(fl1);
  for (auto k = 0u; k < 2u; ++k) {
    SPHERAL_ASSERT_EQ(fl2[k], fl1[k]);
    SPHERAL_ASSERT_EQ(fl2[k]->name(), "stuff1");
  }
}

//------------------------------------------------------------------------------
// appendField/deleteField
//------------------------------------------------------------------------------
TEST_F(FieldListTest, appendAndDeleteField) {
  FieldListDouble fl1(Spheral::FieldStorageType::ReferenceFields);
  FieldDouble stuff("stuff", nodes1);
  fl1.appendField(nodes1.mass());
  fl1.appendField(nodes2.mass());
  SPHERAL_ASSERT_EQ(fl1.numFields(), 2u);
  SPHERAL_ASSERT_TRUE(fl1.haveField(nodes1.mass()));
  SPHERAL_ASSERT_TRUE(fl1.haveField(nodes2.mass()));
  SPHERAL_ASSERT_FALSE(fl1.haveField(stuff));

  fl1.deleteField(nodes2.mass());
  SPHERAL_ASSERT_EQ(fl1.numFields(), 1u);
  SPHERAL_ASSERT_TRUE(fl1.haveField(nodes1.mass()));
  SPHERAL_ASSERT_FALSE(fl1.haveField(nodes2.mass()));

  fl1.deleteField(nodes1.mass());
  SPHERAL_ASSERT_EQ(fl1.numFields(), 0u);
  SPHERAL_ASSERT_FALSE(fl1.haveField(nodes1.mass()));
  SPHERAL_ASSERT_FALSE(fl1.haveField(nodes2.mass()));
}

//------------------------------------------------------------------------------
// appendNewField
//------------------------------------------------------------------------------
TEST_F(FieldListTest, appendNewField) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff", nodes1, 2.0);
  fl1.appendNewField("stuff", nodes2, -3.0);
  SPHERAL_ASSERT_EQ(fl1.numFields(), 2u);
  SPHERAL_ASSERT_EQ(*fl1[0], 2.0);
  SPHERAL_ASSERT_EQ(*fl1[1], -3.0);

  FieldListDouble fl2(Spheral::FieldStorageType::ReferenceFields);
  size_t numFails = 0u, numSucceeds = 0u;
  try {
    fl2.appendNewField("stuff", nodes1, 8.0);
    ++numSucceeds;
  } catch (const Spheral::dbc::VERIFYError&) {
    ++numFails;
  }
  SPHERAL_ASSERT_EQ(numFails, 1u);
  SPHERAL_ASSERT_EQ(numSucceeds, 0u);
  SPHERAL_ASSERT_EQ(fl2.numFields(), 0u);
}

//------------------------------------------------------------------------------
// indexing
//------------------------------------------------------------------------------
TEST_F(FieldListTest, indexing) {
  FieldListDouble fl(Spheral::FieldStorageType::ReferenceFields);
  FieldDouble stuff1("stuff", nodes1), stuff2("stuff", nodes2);

  // We deliberately append these in backwards order to test Field
  // sorting/ordering in the FieldList
  fl.appendField(stuff2);
  fl.appendField(stuff1);
  SPHERAL_ASSERT_EQ(fl.numFields(), 2u);

  // Index for Field pointers
  SPHERAL_ASSERT_EQ(fl[0], &stuff1);
  SPHERAL_ASSERT_EQ(fl[1], &stuff2);
  SPHERAL_ASSERT_EQ(fl.at(0), &stuff1);
  SPHERAL_ASSERT_EQ(fl.at(1), &stuff2);
  
  // Double indexing for individual point values
  for (auto k = 0u; k < 2u; ++k) {
    const auto n = fl[k]->size();
    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                   [&](size_t i) { fl(k,i) = (k+1u)*n + i; });
    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                   [&](size_t i) { SPHERAL_ASSERT_EQ(fl(k,i), (k+1u)*n + i); });
  }
}

//------------------------------------------------------------------------------
// fieldForNodeList
//------------------------------------------------------------------------------
TEST_F(FieldListTest, fieldForNodeList) {
  FieldListDouble fl(Spheral::FieldStorageType::ReferenceFields);
  FieldDouble stuff1("stuff", nodes1), stuff2("stuff", nodes2);

  // We deliberately append these in backwards order to test Field
  // sorting/ordering in the FieldList
  fl.appendField(stuff2);
  fl.appendField(stuff1);
  SPHERAL_ASSERT_EQ(fl.numFields(), 2u);

  SPHERAL_ASSERT_EQ(*fl.fieldForNodeList(nodes1), &stuff1);
  SPHERAL_ASSERT_EQ(*fl.fieldForNodeList(nodes2), &stuff2);
}

//------------------------------------------------------------------------------
// Zero
//------------------------------------------------------------------------------
TEST_F(FieldListTest, zero) {
  FieldListDouble fl(Spheral::FieldStorageType::CopyFields);
  fl.appendNewField("stuff", nodes1, 1.0);
  fl.appendNewField("stuff", nodes2, 2.0);
  SPHERAL_ASSERT_EQ(*fl[0], 1.0);
  SPHERAL_ASSERT_EQ(*fl[1], 2.0);

  fl.Zero();
  SPHERAL_ASSERT_EQ(*fl[0], 0.0);
  SPHERAL_ASSERT_EQ(*fl[1], 0.0);
}  

//------------------------------------------------------------------------------
// Math operations
//------------------------------------------------------------------------------
TEST_F(FieldListTest, math) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields),
                  fl2(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff1", nodes1, 0.0);
  fl1.appendNewField("stuff1", nodes2, 0.0);
  fl2.appendNewField("stuff2", nodes1, 0.0);
  fl2.appendNewField("stuff2", nodes2, 0.0);
  fillRandom(fl1);
  fillRandom(fl2, 1.0);
  const auto nf = fl1.numFields();

  // Addition (FL + FL)
  {
    auto fl3 = fl1 + fl2;
    SPHERAL_ASSERT_EQ(fl3.numFields(), nf);
    for (auto k = 0u; k < 2u; ++k) {
      const auto n = fl1[k]->numInternalElements();
      RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                     [&](size_t i) { SPHERAL_ASSERT_EQ(fl3(k,i), fl1(k,i) + fl2(k,i)); });
    }
  }

  // Subtraction (FL - FL)
  {
    auto fl3 = fl1 - fl2;
    SPHERAL_ASSERT_EQ(fl3.numFields(), nf);
    for (auto k = 0u; k < 2u; ++k) {
      const auto n = fl1[k]->numInternalElements();
      RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                     [&](size_t i) { SPHERAL_ASSERT_EQ(fl3(k,i), fl1(k,i) - fl2(k,i)); });
    }
  }

  // Division (FL / FL)
  {
    auto fl3 = fl1 / fl2;
    SPHERAL_ASSERT_EQ(fl3.numFields(), nf);
    for (auto k = 0u; k < 2u; ++k) {
      const auto n = fl1[k]->numInternalElements();
      RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                     [&](size_t i) { SPHERAL_ASSERT_TRUE(Spheral::fuzzyEqual(fl3(k,i), fl1(k,i) / fl2(k,i), 1.0e-10)); });
    }
  }

  // We need a random number generator to pick some RHS numbers
  std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> realDistribution(0.001, 1e4);

  // Addition (FL + scalar)
  {
    const auto rhs = realDistribution(generator);
    auto fl3 = fl1 + rhs;
    SPHERAL_ASSERT_EQ(fl3.numFields(), nf);
    for (auto k = 0u; k < 2u; ++k) {
      const auto n = fl1[k]->numInternalElements();
      RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                     [&](size_t i) { SPHERAL_ASSERT_EQ(fl3(k,i), fl1(k,i) + rhs); });
    }
  }

  // Subtraction (FL - scalar)
  {
    const auto rhs = realDistribution(generator);
    auto fl3 = fl1 - rhs;
    SPHERAL_ASSERT_EQ(fl3.numFields(), nf);
    for (auto k = 0u; k < 2u; ++k) {
      const auto n = fl1[k]->numInternalElements();
      RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                     [&](size_t i) { SPHERAL_ASSERT_EQ(fl3(k,i), fl1(k,i) - rhs); });
    }
  }

  // Division (FL / scalar)
  {
    const auto rhs = realDistribution(generator);
    auto fl3 = fl1 / rhs;
    SPHERAL_ASSERT_EQ(fl3.numFields(), nf);
    for (auto k = 0u; k < 2u; ++k) {
      const auto n = fl1[k]->numInternalElements();
      RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, n),
                                     [&](size_t i) { SPHERAL_ASSERT_TRUE(Spheral::fuzzyEqual(fl3(k,i), fl1(k,i) / rhs, 1.0e-10)); });
    }
  }
}

//------------------------------------------------------------------------------
// sumElements, min, and max
//
// Note -- due to MPI goofiness with unit testing we only do the per-process
// versions of these.
//------------------------------------------------------------------------------
TEST_F(FieldListTest, SumMinMax) {
  FieldListDouble fl(Spheral::FieldStorageType::CopyFields);
  fl.appendNewField("stuff", nodes1, 0.0);
  fl.appendNewField("stuff", nodes2, 0.0);
  fillRandom(fl);
  const auto nf = fl.numFields();

  // Find the answers
  auto sumAnswer = 0.0;
  auto minAnswer = std::numeric_limits<double>::max();
  auto maxAnswer = std::numeric_limits<double>::lowest();
  for (auto k = 0u; k < nf; ++k) {
    const auto n = fl[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      sumAnswer += fl(k,i);
      minAnswer = std::min(minAnswer, fl(k,i));
      maxAnswer = std::max(maxAnswer, fl(k,i));
    }
  }
  SPHERAL_ASSERT_TRUE(Spheral::fuzzyEqual(fl.localSumElements(), sumAnswer, 1.0e-10));
  SPHERAL_ASSERT_EQ(fl.localMin(), minAnswer);
  SPHERAL_ASSERT_EQ(fl.localMax(), maxAnswer);
}

//------------------------------------------------------------------------------
// Comparison operators
//------------------------------------------------------------------------------
TEST_F(FieldListTest, comparisons) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields),
                  fl2(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff1", nodes1, 0.0);
  fl1.appendNewField("stuff1", nodes2, 0.0);
  fl2.appendNewField("stuff2", nodes1, 0.0);
  fl2.appendNewField("stuff2", nodes2, 0.0);
  auto minVal = -1e3, maxVal = 1e3;
  fillRandom(fl1, minVal, maxVal);
  fillRandom(fl2, minVal, maxVal);

  // FieldList : FieldList
  FieldListDouble fl3(fl1);
  SPHERAL_ASSERT_TRUE(fl1 == fl3);
  SPHERAL_ASSERT_TRUE(fl1 != fl2);

  // FieldList : Value
  fl3 = 1.0;
  auto delta = 0.1*(maxVal - minVal);
  SPHERAL_ASSERT_TRUE(fl3 == 1.0);
  SPHERAL_ASSERT_TRUE(fl1 >  (minVal - delta));
  SPHERAL_ASSERT_TRUE(fl1 <  (maxVal + delta));
  SPHERAL_ASSERT_TRUE(fl1 >= fl1.localMin(true));
  SPHERAL_ASSERT_TRUE(fl1 <= fl1.localMax(true));
  SPHERAL_ASSERT_FALSE(fl1 > fl1.localMin(true));
  SPHERAL_ASSERT_FALSE(fl1 < fl1.localMax(true));
}

//------------------------------------------------------------------------------
// nodeListPtrs
//------------------------------------------------------------------------------
TEST_F(FieldListTest, nodeListPtrs) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff1", nodes2, 0.0);  // Deliberately backwards to check ordering
  fl1.appendNewField("stuff1", nodes1, 0.0);

  auto nodeListPtrs = fl1.nodeListPtrs();
  SPHERAL_ASSERT_TRUE(nodeListPtrs.size() == 2u);
  SPHERAL_ASSERT_TRUE(nodeListPtrs[0] == &nodes1);
  SPHERAL_ASSERT_TRUE(nodeListPtrs[1] == &nodes2);
  SPHERAL_ASSERT_TRUE(std::find(nodeListPtrs.begin(), nodeListPtrs.end(), &nodes3) == nodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Flattening operations
//------------------------------------------------------------------------------
TEST_F(FieldListTest, flatten) {
  FieldListDouble fl1(Spheral::FieldStorageType::CopyFields);
  fl1.appendNewField("stuff1", nodes1, 0.0);
  fl1.appendNewField("stuff1", nodes2, 0.0);
  fillRandom(fl1);

  // internal
  {
    const auto vals = fl1.internalValues();
    SPHERAL_ASSERT_EQ(vals.size(), nodes1.numInternalNodes() + nodes2.numInternalNodes());
    const auto n1 = nodes1.numInternalNodes();
    const auto n = vals.size();
    RAJA::forall<LOOP_EXEC_POLICY>
      (TRS_UINT(0, n),
       [&](size_t i) {
         const size_t k = i < n1 ? 0 : 1;
         const size_t j = i - k*n1;
         // printf("%zu : (%zu, %zu)\n", i, k, j);
         SPHERAL_ASSERT_EQ(vals[i], fl1(k,j));
       });
  }

  // ghost
  {
    const auto vals = fl1.ghostValues();
    SPHERAL_ASSERT_EQ(vals.size(), nodes1.numGhostNodes() + nodes2.numGhostNodes());
    const auto n1 = nodes1.numGhostNodes();
    const auto n = vals.size();
    RAJA::forall<LOOP_EXEC_POLICY>
      (TRS_UINT(0, n),
       [&](size_t i) {
         const size_t k = i < n1 ? 0 : 1;
         const size_t j = i - k*n1;
         SPHERAL_ASSERT_EQ(vals[i], fl1(k,j));
       });
  }

  // all
  {
    const auto vals = fl1.allValues();
    SPHERAL_ASSERT_EQ(vals.size(), nodes1.numNodes() + nodes2.numNodes());
    const auto n1 = nodes1.numNodes();
    const auto n = vals.size();
    RAJA::forall<LOOP_EXEC_POLICY>
      (TRS_UINT(0, n),
       [&](size_t i) {
         const size_t k = i < n1 ? 0 : 1;
         const size_t j = i - k*n1;
         SPHERAL_ASSERT_EQ(vals[i], fl1(k,j));
       });
  }

}
