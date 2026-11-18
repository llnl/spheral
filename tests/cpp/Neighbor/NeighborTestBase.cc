//------------------------------------------------------------------------------
// C++ counterpart to tests/unit/Neighbor/NeighborTestBase.py.
//
// The ordering, names, full-mode setup values, and test bodies follow the Python
// implementation.  The default uses smaller populations; --full restores the
// original workload.
//------------------------------------------------------------------------------

#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "DataBase/DataBase.hh"
#include "Kernel/BSplineKernel.hh"
#include "Kernel/TableKernel.hh"
#include "Material/GammaLawGas.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/NodePairList.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "NodeList/FluidNodeList.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace Spheral;

bool fullTest = false;

// Match Python's process-global random module state.
std::mt19937 globalRandomGenerator(std::random_device{}());

//=============================================================================
// Apply a random rotation operation to the given symmetric tensor.
//=============================================================================
template<typename Dimension>
void
applyRandomRotation(typename Dimension::SymTensor& H,
                    std::mt19937& generator) {
  constexpr auto nDimensions = Dimension::nDim;
  std::uniform_real_distribution<double> angle(0.0, std::acos(-1.0));
  if constexpr (nDimensions == 2) {
    const auto theta = angle(generator);
    const typename Dimension::Tensor R(std::cos(theta), -std::sin(theta),
                                       std::sin(theta),  std::cos(theta));
    H.rotationalTransform(R);
  } else if constexpr (nDimensions == 3) {
    const auto theta1 = angle(generator);
    const auto theta2 = angle(generator);
    const auto theta3 = angle(generator);
    const typename Dimension::Tensor R1(std::cos(theta1), -std::sin(theta1), 0.0,
                                        std::sin(theta1),  std::cos(theta1), 0.0,
                                        0.0,               0.0,              1.0);
    const typename Dimension::Tensor R2(std::cos(theta2), 0.0, -std::sin(theta2),
                                        0.0,              1.0,  0.0,
                                        std::sin(theta2), 0.0,  std::cos(theta2));
    const typename Dimension::Tensor R3(1.0, 0.0,               0.0,
                                        0.0, std::cos(theta3), -std::sin(theta3),
                                        0.0, std::sin(theta3),  std::cos(theta3));
    const auto R = R1*R2*R3;
    H.rotationalTransform(R);
  }
}

//=============================================================================
// Base class for the tests, providing the generic setup methods.
//=============================================================================
template<typename Dimension>
class NeighborTestBase: public ::testing::Test {
public:
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using NodeList = Spheral::NodeList<Dimension>;
  using FluidNodeList = Spheral::FluidNodeList<Dimension>;
  using DataBase = Spheral::DataBase<Dimension>;
  using NodeID = typename Spheral::Neighbor<Dimension>::NodeID;
  using TableKernel = Spheral::TableKernel<Dimension>;
  using GammaLawGas = Spheral::GammaLawGas<Dimension>;
  using VolumeRange = std::array<std::array<double, 2>, Dimension::nDim>;

  static constexpr double kernelExtent = 2.0;

  struct RandomDistribution {
    std::vector<size_t> globalNodeIDs;
    std::vector<Vector> nodePositions;
    std::vector<SymTensor> H;
  };

  struct NodeSet {
    FluidNodeList nodes;

    NodeSet(const std::string& name,
            GammaLawGas& eos):
      nodes(name, eos, 0u, 0u,
            1.0e-20, 1.0e20, 0.1, 2.01, 500u,
            1.0e-10, 1.0e10) {
    }
  };

  NeighborTestBase():
    units(1.0, 1.0, 1.0),
    ncheck(0u),
    noverlapcheck(0u),
    generator(globalRandomGenerator) {
  }

  //---------------------------------------------------------------------------
  // Distribute nodes randomly amongst domains.
  //---------------------------------------------------------------------------
  RandomDistribution
  randomDistribute(const size_t n,
                   const VolumeRange& volrange,
                   const double nNodesPerh = 2.01) {
    constexpr auto nDimensions = Dimension::nDim;
    const auto adim = 1.0/double(nDimensions);
    auto volume = 1.0;
    for (auto i = 0; i < nDimensions; ++i) {
      REQUIRE(volrange[i][1] > volrange[i][0]);
      volume *= volrange[i][1] - volrange[i][0];
    }
    REQUIRE(volume > 0.0);
    const auto dx0 = std::pow(volume/double(n), adim);

    RandomDistribution result;
    result.globalNodeIDs.reserve(n);
    result.nodePositions.reserve(n);
    result.H.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      result.globalNodeIDs.push_back(i);
      Vector position;
      for (auto j = 0; j < nDimensions; ++j) {
        std::uniform_real_distribution<double> coordinate(volrange[j][0],
                                                          volrange[j][1]);
        position[j] = coordinate(generator);
      }
      result.nodePositions.push_back(position);

      Vector inverseDx;
      for (auto j = 0; j < nDimensions; ++j) {
        std::uniform_real_distribution<double> scale(0.5, 2.0);
        const auto dx = nNodesPerh*scale(generator)*dx0;
        REQUIRE(dx > 0.0);
        inverseDx[j] = 1.0/dx;
      }
      SymTensor Hi(inverseDx);
      applyRandomRotation<Dimension>(Hi, generator);
      result.H.push_back(Hi);
    }
    REQUIRE(result.globalNodeIDs.size() == n);
    REQUIRE(result.nodePositions.size() == n);
    REQUIRE(result.H.size() == n);
    return result;
  }

  //---------------------------------------------------------------------------
  // Setup routine common to all instances.
  //---------------------------------------------------------------------------
  void
  genericSetUp(const size_t n1,
               const size_t n2,
               const size_t n3,
               const VolumeRange& volrange1,
               const VolumeRange& volrange2,
               const VolumeRange& volrange3) {
    // Construct the NodeLists to be distributed.
    eos = std::make_unique<GammaLawGas>(2.0, 2.0, units,
                                       -1.0e200, 1.0e200,
                                       MaterialPressureMinType::PressureFloor,
                                       0.0);
    WT = std::make_unique<TableKernel>(BSplineKernel<Dimension>(), 100u);
    nodes1 = std::make_unique<NodeSet>("nodes 1", *eos);
    nodes2 = std::make_unique<NodeSet>("nodes 2", *eos);
    nodes3 = std::make_unique<NodeSet>("nodes 3", *eos);

    const std::array<NodeSet*, 3> nodeSets = {nodes1.get(), nodes2.get(), nodes3.get()};
    const std::array<size_t, 3> globalCounts = {n1, n2, n3};
    const std::array<const VolumeRange*, 3> volumeRanges = {
      &volrange1, &volrange2, &volrange3
    };
    for (auto nodeListID = 0u; nodeListID < nodeSets.size(); ++nodeListID) {
      auto& nodes = nodeSets[nodeListID]->nodes;
      const auto distribution = randomDistribute(globalCounts[nodeListID],
                                                  *volumeRanges[nodeListID]);
      const auto n = distribution.globalNodeIDs.size();
      nodes.numInternalNodes(n);
      for (size_t i = 0; i < n; ++i) {
        nodes.mass()[i] = 1.0;
        nodes.positions()[i] = distribution.nodePositions[i];
        nodes.Hfield()[i] = distribution.H[i];
      }
    }

    // Put the distributed NodeLists into a DataBase.
    dataBase = std::make_unique<DataBase>(kernelExtent);
    dataBase->appendNodeList(nodes1->nodes);
    dataBase->appendNodeList(nodes2->nodes);
    dataBase->appendNodeList(nodes3->nodes);

    // Update the DataBase-owned Neighbor info.
    dataBase->updateNeighbors();
  }

  //---------------------------------------------------------------------------
  // C++ equivalents of the SpheralTestUtilities neighbor helpers.
  //---------------------------------------------------------------------------
  std::vector<NodeID>
  findCandidateNodes(const Vector& position,
                     const SymTensor& H) const {
    std::vector<NodeID> master, coarse, result;
    const auto neighborVariant = dataBase->neighborView();
    visitNeighborView(neighborVariant, [&](const auto& neighbor) {
      const auto [ok, numMaster, numCoarse] = neighbor.setMasterList(position, H, master, coarse);
      VERIFY2(ok, "Insufficient storage while constructing the master/coarse neighbor lists.");
      master.resize(numMaster);
      coarse.resize(numCoarse);
      result.resize(numCoarse);
      result.resize(neighbor.setRefineNeighborList(position, H, numCoarse, coarse, result));
    });
    return result;
  }

  std::vector<size_t>
  findNeighborNodes(const Vector& position,
                    const SymTensor& H,
                    const NodeList& nodes) const {
    std::vector<size_t> result;
    for (auto i = 0u; i < nodes.numNodes(); ++i) {
      const auto displacement = nodes.positions()[i] - position;
      if (std::min((H*displacement).magnitude(),
                   (nodes.Hfield()[i]*displacement).magnitude()) <= kernelExtent) {
        result.push_back(i);
      }
    }
    return result;
  }

  std::vector<size_t>
  findGatherNeighborNodes(const Vector& position,
                          const SymTensor& H,
                          const NodeList& nodes) const {
    std::vector<size_t> result;
    for (auto i = 0u; i < nodes.numNodes(); ++i) {
      if ((H*(position - nodes.positions()[i])).magnitude() <= kernelExtent) {
        result.push_back(i);
      }
    }
    return result;
  }

  std::vector<size_t>
  findScatterNeighborNodes(const Vector& position,
                           const NodeList& nodes) const {
    std::vector<size_t> result;
    for (auto i = 0u; i < nodes.numNodes(); ++i) {
      if ((nodes.Hfield()[i]*(position - nodes.positions()[i])).magnitude() <= kernelExtent) {
        result.push_back(i);
      }
    }
    return result;
  }

  std::vector<std::vector<size_t>>
  findOverlapNeighbors(const Vector& position,
                       const SymTensor& H) const {
    const auto& nodeLists = dataBase->nodeListPtrs();
    std::vector<std::set<size_t>> uniqueResult(nodeLists.size());
    for (size_t nodeListID = 0; nodeListID < nodeLists.size(); ++nodeListID) {
      const auto neighbors = findNeighborNodes(position, H, *nodeLists[nodeListID]);
      uniqueResult[nodeListID].insert(neighbors.begin(), neighbors.end());
    }
    for (const auto* inodes: nodeLists) {
      const auto gather = findGatherNeighborNodes(position, H, *inodes);
      for (const auto i: gather) {
        for (size_t nodeListID = 0; nodeListID < nodeLists.size(); ++nodeListID) {
          const auto scatter = findScatterNeighborNodes(inodes->positions()[i],
                                                        *nodeLists[nodeListID]);
          uniqueResult[nodeListID].insert(scatter.begin(), scatter.end());
        }
      }
    }
    std::vector<std::vector<size_t>> result(nodeLists.size());
    for (size_t nodeListID = 0; nodeListID < nodeLists.size(); ++nodeListID) {
      result[nodeListID].assign(uniqueResult[nodeListID].begin(),
                                uniqueResult[nodeListID].end());
    }
    return result;
  }

  static bool
  checkNeighbors(const std::vector<size_t>& neighborList,
                 const std::vector<size_t>& answer) {
    return std::all_of(answer.begin(), answer.end(), [&](const auto value) {
      return std::find(neighborList.begin(), neighborList.end(), value) != neighborList.end();
    });
  }

  std::vector<size_t>
  randomSample(const size_t populationSize,
               const size_t sampleSize) {
    REQUIRE(sampleSize <= populationSize);
    if (sampleSize == 0u) return {};
    std::vector<size_t> result(populationSize);
    std::iota(result.begin(), result.end(), 0u);
    std::shuffle(result.begin(), result.end(), generator);
    result.resize(sampleSize);
    return result;
  }

  // In the default mode, explicitly cover both x boundaries of each NodeList
  // (including the interfaces between lists), then retain random interior nodes.
  // Full mode keeps the original sampling behavior.
  std::vector<size_t>
  sampleNodes(const NodeList& nodes, const size_t sampleSize) {
    if (fullTest) return randomSample(nodes.numInternalNodes() - 1u, sampleSize);
    auto result = randomSample(nodes.numInternalNodes(), sampleSize);
    if (result.empty()) return result;
    size_t lower = 0u, upper = 0u;
    for (size_t i = 1u; i < nodes.numInternalNodes(); ++i) {
      if (nodes.positions()[i][0] < nodes.positions()[lower][0]) lower = i;
      if (nodes.positions()[i][0] > nodes.positions()[upper][0]) upper = i;
    }
    const auto include = [&](const size_t node, const size_t slot) {
      const auto itr = std::find(result.begin(), result.end(), node);
      if (itr == result.end()) {
        result[slot] = node;
      } else {
        std::iter_swap(result.begin() + slot, itr);
      }
    };
    include(lower, 0u);
    if (result.size() > 1u and upper != lower) include(upper, 1u);
    return result;
  }

  // Cover every available within-list and cross-list pair category before
  // filling the remaining slots with random pairs in the default mode.
  std::vector<size_t>
  samplePairs(const std::vector<NodePairIdxType>& pairs) {
    const auto random = randomSample(pairs.size(), ncheck);
    if (fullTest) return random;
    std::vector<size_t> result;
    std::set<std::pair<size_t, size_t>> categories;
    for (size_t i = 0u; i < pairs.size() and result.size() < ncheck; ++i) {
      if (categories.emplace(std::min(pairs[i].i_list, pairs[i].j_list),
                             std::max(pairs[i].i_list, pairs[i].j_list)).second) {
        result.push_back(i);
      }
    }
    for (const auto i: random) {
      if (result.size() == ncheck) break;
      if (std::find(result.begin(), result.end(), i) == result.end()) result.push_back(i);
    }
    return result;
  }

  static std::vector<size_t>
  sorted(std::vector<size_t> values) {
    std::sort(values.begin(), values.end());
    return values;
  }

  //---------------------------------------------------------------------------
  // Exercise TreeNeighborView list construction in each configured execution
  // space.  Fixed-size ManagedArrays provide device-accessible output storage;
  // the View methods report the number of populated entries separately.
  //---------------------------------------------------------------------------
  template<typename ExecPolicy>
  void
  testTreeNeighborView() {
    const auto neighborVariant = dataBase->neighborView();
    const auto treeView = std::get<Spheral::TreeNeighborView<Dimension>>(neighborVariant);
    const auto& nodes = nodes1->nodes;
    ASSERT_GT(nodes.numInternalNodes(), 0u);

    const auto position = nodes.positions()[0u];
    const auto H = nodes.Hfield()[0u];
    size_t capacity = 0u;
    for (const auto* nodeList: dataBase->nodeListPtrs()) capacity += nodeList->numNodes();
    ASSERT_GT(capacity, 0u);

    chai::ManagedArray<NodeID> master(capacity);
    chai::ManagedArray<NodeID> coarse(capacity);
    chai::ManagedArray<NodeID> refine(capacity);
    RAJA::forall<ExecPolicy>(TRS_UINT(0u, 1u),
      [=] SPHERAL_HOST_DEVICE(size_t) mutable {
        const auto [ok, numMaster, numCoarse] =
          treeView.setMasterList(position, H, master, coarse);
        SPHERAL_ASSERT_TRUE(ok);
        SPHERAL_ASSERT_TRUE(numMaster <= capacity);
        SPHERAL_ASSERT_TRUE(numCoarse <= capacity);
        const auto numRefine =
          treeView.setRefineNeighborList(position, H, numCoarse, coarse, refine);
        SPHERAL_ASSERT_TRUE(numRefine <= capacity);

        auto foundExpectedNode = false;
        for (size_t k = 0u; k < numRefine; ++k) {
          if (refine[k].first == 0u and refine[k].second == 0u) foundExpectedNode = true;
        }
        SPHERAL_ASSERT_TRUE(foundExpectedNode);
      });

    master.free();
    coarse.free();
    refine.free();
  }

  //---------------------------------------------------------------------------
  // Test raw Neighbor object calls.
  //---------------------------------------------------------------------------
  void
  testNeighborCalls() {
    // Iterate over the NodeLists.
    for (const auto* nodes: dataBase->nodeListPtrs()) {
      const auto& position = nodes->positions();
      const auto& H = nodes->Hfield();

      // Select boundary and random nodes in default mode, random nodes in full mode.
      for (const auto nodeID: sampleNodes(*nodes, ncheck)) {
        const auto ri = position[nodeID];
        const auto Hi = H[nodeID];

        // Have the DataBase-owned Neighbor select neighbors for this node.
        const auto refineNeighbors = findCandidateNodes(ri, Hi);
        std::vector<size_t> neighborIDs;
        const auto& nodeLists = dataBase->nodeListPtrs();
        std::vector<size_t> offsets(nodeLists.size(), 0u);
        for (auto nodeListID = 1u; nodeListID < nodeLists.size(); ++nodeListID) {
          offsets[nodeListID] = offsets[nodeListID - 1u] + nodeLists[nodeListID - 1u]->numInternalNodes();
        }
        for (const auto [nodeListID, i]: refineNeighbors) {
          neighborIDs.push_back(offsets[nodeListID] + i);
        }

        // Now build the checks.
        std::vector<size_t> answerIDs;
        auto nodeListID = 0u;
        for (const auto* answerNodes: nodeLists) {
          for (const auto i: findNeighborNodes(ri, Hi, *answerNodes)) {
            answerIDs.push_back(i + offsets[nodeListID]);
          }
          ++nodeListID;
        }

        // Check the answer.
        EXPECT_TRUE(checkNeighbors(neighborIDs, answerIDs)) << "node " << nodeID;
      }
    }
  }

  //---------------------------------------------------------------------------
  // Test ConnectivityMap neighbors.
  //---------------------------------------------------------------------------
  void
  testConnectivityMapNeighbors() {
    dataBase->updateConnectivityMap(false, false, false);
    const auto& cm = dataBase->connectivityMap(false, false, false);
    const auto& nodeLists = dataBase->nodeListPtrs();

    // Iterate over the NodeLists.
    for (size_t iNL = 0; iNL < nodeLists.size(); ++iNL) {
      const auto* inodes = nodeLists[iNL];
      const auto& position = inodes->positions();
      const auto& H = inodes->Hfield();

      // Select boundary and random nodes in default mode, random nodes in full mode.
      for (const auto i: sampleNodes(*inodes, ncheck)) {
        const auto ri = position[i];
        const auto Hi = H[i];
        const auto& cmneighbors = cm.connectivityForNode(inodes, i);

        // Check ConnectivityMap vs. N^2 neighbor search.
        for (size_t jNL = 0; jNL < nodeLists.size(); ++jNL) {
          const auto cmcheck = sorted(cmneighbors[jNL]);
          auto answer = sorted(findNeighborNodes(ri, Hi, *nodeLists[jNL]));
          if (iNL == jNL) {
            const auto self = std::find(answer.begin(), answer.end(), i);
            ASSERT_NE(self, answer.end());
            answer.erase(self);
          }
          EXPECT_EQ(cmcheck, answer) << "node " << i;
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  // Test ConnectivityMap NodePairList.
  //---------------------------------------------------------------------------
  void
  testConnectivityMapNodePairs() {
    dataBase->updateConnectivityMap(false, false, false);
    const auto& cm = dataBase->connectivityMap(false, false, false);
    const auto& nodeLists = dataBase->nodeListPtrs();
    const auto numNodeLists = dataBase->numNodeLists();

    // Build the answer based on the node neighbors.
    std::vector<NodePairIdxType> answer;
    for (size_t iNL = 0; iNL < numNodeLists; ++iNL) {
      const auto* inodes = nodeLists[iNL];
      for (size_t i = 0; i < inodes->numInternalNodes(); ++i) {
        const auto& cmneighbors = cm.connectivityForNode(inodes, i);
        ASSERT_EQ(cmneighbors.size(), numNodeLists);
        for (size_t jNL = iNL; jNL < numNodeLists; ++jNL) {
          for (const auto j: cmneighbors[jNL]) {
            if (jNL > iNL or j > i) answer.emplace_back(i, iNL, j, jNL);
          }
        }
      }
    }

    std::vector<NodePairIdxType> pairs(cm.nodePairList().begin(),
                                       cm.nodePairList().end());
    ASSERT_EQ(pairs.size(), answer.size());
    std::sort(answer.begin(), answer.end());
    std::sort(pairs.begin(), pairs.end());
    EXPECT_EQ(answer, pairs);
  }

  //---------------------------------------------------------------------------
  // Test ConnectivityMap overlap neighbors.
  //---------------------------------------------------------------------------
  void
  testConnectivityMapOverlapNeighbors() {
    dataBase->updateConnectivityMap(false, true, false);
    const auto& cm = dataBase->connectivityMap(false, true, false);
    const auto& nodeLists = dataBase->nodeListPtrs();

    // Iterate over the NodeLists.
    for (size_t iNL = 0; iNL < nodeLists.size(); ++iNL) {
      const auto* inodes = nodeLists[iNL];

      // Select boundary nodes in default mode, random nodes in full mode.
      for (const auto i: sampleNodes(*inodes, noverlapcheck)) {
        const auto ri = inodes->positions()[i];
        const auto Hi = inodes->Hfield()[i];
        const auto& cmneighbors = cm.overlapConnectivityForNode(inodes, i);
        const auto fullanswer = findOverlapNeighbors(ri, Hi);

        // Check ConnectivityMap vs. N^2 neighbor search.
        for (size_t jNL = 0; jNL < nodeLists.size(); ++jNL) {
          const auto cmcheck = sorted(cmneighbors[jNL]);
          auto answer = sorted(fullanswer[jNL]);
          if (iNL == jNL) {
            const auto self = std::find(answer.begin(), answer.end(), i);
            ASSERT_NE(self, answer.end());
            answer.erase(self);
          }
          EXPECT_EQ(cmcheck, answer) << "node " << i;
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  // Test ConnectivityMap intersection methods.
  //---------------------------------------------------------------------------
  void
  testConnectivityComputeIntersection() {
    dataBase->updateConnectivityMap(false, true, true);
    const auto numNodeLists = dataBase->numNodeLists();
    const auto& cm = dataBase->connectivityMap(false, true, true);
    const auto& nodePairList = cm.nodePairList();
    std::vector<NodePairIdxType> pairs(nodePairList.begin(), nodePairList.end());

    // Compute the actual intersection of two points for each NodeList.
    const auto intersectionAnswer = [&](const NodePairIdxType& pair) {
      const auto& neighborsi = cm.connectivityForNode(pair.i_list, pair.i_node);
      const auto& neighborsj = cm.connectivityForNode(pair.j_list, pair.j_node);
      std::vector<std::set<size_t>> result(numNodeLists);
      for (size_t kNL = 0; kNL < numNodeLists; ++kNL) {
        const std::set<size_t> seti(neighborsi[kNL].begin(), neighborsi[kNL].end());
        const std::set<size_t> setj(neighborsj[kNL].begin(), neighborsj[kNL].end());
        std::set_intersection(seti.begin(), seti.end(),
                              setj.begin(), setj.end(),
                              std::inserter(result[kNL], result[kNL].begin()));
      }
      result[pair.i_list].insert(pair.i_node);
      result[pair.j_list].insert(pair.j_node);
      return result;
    };

    // Make vector<vector<size_t>> into a vector of sets for comparison.
    const auto convertToSets = [&](const auto& neighbors) {
      EXPECT_EQ(neighbors.size(), numNodeLists);
      std::vector<std::set<size_t>> result;
      result.reserve(neighbors.size());
      for (const auto& values: neighbors) result.emplace_back(values.begin(), values.end());
      return result;
    };

    // Include within-list and cross-list pairs in default mode.
    for (const auto pairIndex: samplePairs(pairs)) {
      const auto& pair = pairs[pairIndex];
      const auto answer = intersectionAnswer(pair);
      const auto intersect = convertToSets(
        cm.connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                            pair.j_list, pair.j_node));
      ASSERT_EQ(intersect.size(), numNodeLists);
      EXPECT_EQ(intersect, answer) << "pair " << pair;
    }
  }

protected:
  PhysicalConstants units;
  std::unique_ptr<GammaLawGas> eos;
  std::unique_ptr<TableKernel> WT;
  std::unique_ptr<NodeSet> nodes1, nodes2, nodes3;
  std::unique_ptr<DataBase> dataBase;
  size_t ncheck, noverlapcheck;
  std::mt19937& generator;
};

//=============================================================================
// Random node distribution -- 1-D.
//=============================================================================
class NeighborRandom1d: public NeighborTestBase<Dim<1>> {
public:
  void SetUp() override {
    ncheck = 10u;
    noverlapcheck = 2u;
    const auto n1 = fullTest ? 1000u : 40u;
    const auto n2 = fullTest ? 2500u : 100u;
    const auto n3 = fullTest ? 500u : 20u;
    const VolumeRange volrange1 = {{{-2.0, -1.0}}};
    const VolumeRange volrange2 = {{{-1.0,  0.5}}};
    const VolumeRange volrange3 = {{{ 0.5,  2.0}}};
    genericSetUp(n1, n2, n3, volrange1, volrange2, volrange3);
  }
};

TEST_F(NeighborRandom1d, testNeighborCalls) {
  this->testNeighborCalls();
}

TEST_F(NeighborRandom1d, testConnectivityMapNeighbors) {
  this->testConnectivityMapNeighbors();
}

TEST_F(NeighborRandom1d, testConnectivityMapNodePairs) {
  this->testConnectivityMapNodePairs();
}

TEST_F(NeighborRandom1d, testConnectivityMapOverlapNeighbors) {
  this->testConnectivityMapOverlapNeighbors();
}

TEST_F(NeighborRandom1d, testConnectivityComputeIntersection) {
  this->testConnectivityComputeIntersection();
}

//=============================================================================
// Random node distribution -- 2-D.
//=============================================================================
class NeighborRandom2d: public NeighborTestBase<Dim<2>> {
public:
  void SetUp() override {
    ncheck = 10u;
    noverlapcheck = 2u;
    const auto n1 = fullTest ? 10000u : 80u;
    const auto n2 = fullTest ? 25000u : 200u;
    const auto n3 = fullTest ? 5000u : 40u;
    const VolumeRange volrange1 = {{{-2.0, -1.0}, {0.0, 1.0}}};
    const VolumeRange volrange2 = {{{-1.0, -0.5}, {0.0, 1.0}}};
    const VolumeRange volrange3 = {{{-0.5,  0.0}, {0.0, 1.0}}};
    genericSetUp(n1, n2, n3, volrange1, volrange2, volrange3);
  }
};

TEST_F(NeighborRandom2d, testNeighborCalls) {
  this->testNeighborCalls();
}

TEST_F(NeighborRandom2d, testConnectivityMapNeighbors) {
  this->testConnectivityMapNeighbors();
}

TEST_F(NeighborRandom2d, testConnectivityMapNodePairs) {
  this->testConnectivityMapNodePairs();
}

TEST_F(NeighborRandom2d, testConnectivityMapOverlapNeighbors) {
  this->testConnectivityMapOverlapNeighbors();
}

TEST_F(NeighborRandom2d, testConnectivityComputeIntersection) {
  this->testConnectivityComputeIntersection();
}

//=============================================================================
// Random node distribution -- 3-D.
//=============================================================================
class NeighborRandom3d: public NeighborTestBase<Dim<3>> {
public:
  void SetUp() override {
    ncheck = 10u;
    noverlapcheck = fullTest ? 0u : 2u;
    const auto n1 = fullTest ? 1000u : 64u;
    const auto n2 = fullTest ? 2500u : 160u;
    const auto n3 = fullTest ? 1500u : 96u;
    const VolumeRange volrange1 = {{{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}};
    const VolumeRange volrange2 = {{{1.0, 1.5}, {0.0, 1.0}, {0.0, 1.0}}};
    const VolumeRange volrange3 = {{{1.5, 2.0}, {0.0, 1.0}, {0.0, 1.0}}};
    genericSetUp(n1, n2, n3, volrange1, volrange2, volrange3);
  }
};

// Run TreeNeighborView itself in every enabled execution space, including the
// device policy in CUDA/HIP builds.
TYPED_TEST_SUITE_P(TreeNeighborViewTypedTest);
template<typename T>
class TreeNeighborViewTypedTest: public NeighborRandom3d {};

GPU_TYPED_TEST_P(TreeNeighborViewTypedTest, DeviceNeighborQuery) {
  gpu_this->template testTreeNeighborView<TypeParam>();
}

REGISTER_TYPED_TEST_SUITE_P(TreeNeighborViewTypedTest,
                            DeviceNeighborQuery);

INSTANTIATE_TYPED_TEST_SUITE_P(TreeNeighborView,
                               TreeNeighborViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );

TEST_F(NeighborRandom3d, testNeighborCalls) {
  this->testNeighborCalls();
}

TEST_F(NeighborRandom3d, testConnectivityMapNeighbors) {
  this->testConnectivityMapNeighbors();
}

TEST_F(NeighborRandom3d, testConnectivityMapNodePairs) {
  this->testConnectivityMapNodePairs();
}

TEST_F(NeighborRandom3d, testConnectivityMapOverlapNeighbors) {
  this->testConnectivityMapOverlapNeighbors();
}

TEST_F(NeighborRandom3d, testConnectivityComputeIntersection) {
  this->testConnectivityComputeIntersection();
}

//=============================================================================
// Cylindrical node distribution -- 2-D.
//=============================================================================
class NeighborCylindrical2d: public NeighborTestBase<Dim<2>> {
public:
  void SetUp() override {
    ncheck = fullTest ? 50u : 20u;
    noverlapcheck = 2u;

    // Match GenerateNodeDistribution2d(..., distributionType="constantDTheta").
    eos = std::make_unique<GammaLawGas>(2.0, 2.0, units,
                                       -1.0e200, 1.0e200,
                                       MaterialPressureMinType::PressureFloor,
                                       0.0);
    WT = std::make_unique<TableKernel>(BSplineKernel<Dim<2>>(), 100u);
    nodes1 = std::make_unique<NodeSet>("cylindrical nodes 1", *eos);

    const auto nRadial = fullTest ? 100u : 12u;
    constexpr auto rmin = 0.0;
    constexpr auto rmax = 1.0;
    constexpr auto theta = 0.5*3.14159265358979323846;
    constexpr auto nNodesPerh = 2.01;
    constexpr auto rho = 1.0;
    const auto dr = (rmax - rmin)/double(nRadial);
    const auto h = 1.0/(nNodesPerh*dr);
    const SymTensor Hi(h, 0.0,
                       0.0, h);

    std::vector<Vector> positions;
    std::vector<double> masses;
    for (auto i = 0u; i < nRadial; ++i) {
      const auto rInner = rmin + i*dr;
      const auto rOuter = rmin + (i + 1u)*dr;
      const auto ri = rmin + (i + 0.5)*dr;
      const auto nominalNTheta = int(theta*ri/dr);
      const auto nTheta = std::max(1, nominalNTheta);
      const auto dTheta = theta/double(nTheta);
      const auto mRing = (rOuter*rOuter - rInner*rInner)*theta*0.5*rho;
      const auto mi = mRing/double(nTheta);
      for (auto j = 0; j < nTheta; ++j) {
        const auto thetai = (j + 0.5)*dTheta;
        positions.emplace_back(ri*std::cos(thetai), ri*std::sin(thetai));
        masses.push_back(mi);
      }
    }

    auto& nodes = nodes1->nodes;
    nodes.numInternalNodes(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
      nodes.mass()[i] = masses[i];
      nodes.positions()[i] = positions[i];
      nodes.Hfield()[i] = Hi;
    }

    dataBase = std::make_unique<DataBase>(kernelExtent);
    dataBase->appendNodeList(nodes);
    dataBase->updateNeighbors();
  }
};

TEST_F(NeighborCylindrical2d, testNeighborCalls) {
  this->testNeighborCalls();
}

TEST_F(NeighborCylindrical2d, testConnectivityMapNeighbors) {
  this->testConnectivityMapNeighbors();
}

TEST_F(NeighborCylindrical2d, testConnectivityMapNodePairs) {
  this->testConnectivityMapNodePairs();
}

TEST_F(NeighborCylindrical2d, testConnectivityMapOverlapNeighbors) {
  this->testConnectivityMapOverlapNeighbors();
}

TEST_F(NeighborCylindrical2d, testConnectivityComputeIntersection) {
  this->testConnectivityComputeIntersection();
}

} // anonymous namespace

int main(int argc, char* argv[]) {
  // Remove our option before GoogleTest parses its own flags.  Keep argv
  // compact and null-terminated so filters, repeats, and XML output still work.
  auto remaining = 1;
  for (auto i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--full") {
      fullTest = true;
    } else {
      argv[remaining++] = argv[i];
    }
  }
  argc = remaining;
  argv[argc] = nullptr;
  ::testing::InitGoogleTest(&argc, argv);
  std::cout << "NeighborTestBase mode: " << (fullTest ? "full" : "default")
            << " (use --full for the original workload)\n";
  return RUN_ALL_TESTS();
}
