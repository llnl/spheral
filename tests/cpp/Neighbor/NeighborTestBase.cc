//------------------------------------------------------------------------------
// C++ counterpart to tests/unit/Neighbor/NeighborTestBase.py.
//
// The ordering, names, setup values, and test bodies intentionally follow the
// Python implementation so the two versions can be compared side by side.
//------------------------------------------------------------------------------

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
#include <iterator>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {

using namespace Spheral;

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
  using TreeNeighbor = Spheral::TreeNeighbor<Dimension>;
  using DataBase = Spheral::DataBase<Dimension>;
  using TableKernel = Spheral::TableKernel<Dimension>;
  using GammaLawGas = Spheral::GammaLawGas<Dimension>;
  using VolumeRange = std::array<std::array<double, 2>, Dimension::nDim>;

  static constexpr double kernelExtent = 2.0;

  struct RandomDistribution {
    std::vector<int> globalNodeIDs;
    std::vector<Vector> nodePositions;
    std::vector<SymTensor> H;
  };

  struct NodeSet {
    FluidNodeList nodes;
    TreeNeighbor neighbor;

    NodeSet(const std::string& name,
            GammaLawGas& eos,
            const Vector& xmin,
            const Vector& xmax):
      nodes(name, eos, 0u, 0u,
            1.0e-20, 1.0e20, 0.1, 2.01, 500u,
            1.0e-10, 1.0e10),
      neighbor(nodes, NeighborSearchType::GatherScatter,
               kernelExtent, xmin, xmax) {
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
      result.globalNodeIDs.push_back(int(i));
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
    const auto xmin = Vector(-10.0);
    const auto xmax = Vector(10.0);
    nodes1 = std::make_unique<NodeSet>("nodes 1", *eos, xmin, xmax);
    nodes2 = std::make_unique<NodeSet>("nodes 2", *eos, xmin, xmax);
    nodes3 = std::make_unique<NodeSet>("nodes 3", *eos, xmin, xmax);

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
    dataBase = std::make_unique<DataBase>();
    dataBase->appendNodeList(nodes1->nodes);
    dataBase->appendNodeList(nodes2->nodes);
    dataBase->appendNodeList(nodes3->nodes);

    // Update all the Neighbor info.
    for (auto* nodeSet: nodeSets) nodeSet->nodes.neighbor().updateNodes();
  }

  //---------------------------------------------------------------------------
  // C++ equivalents of the SpheralTestUtilities neighbor helpers.
  //---------------------------------------------------------------------------
  std::vector<int>
  findNeighborNodes(const Vector& position,
                    const SymTensor& H,
                    const NodeList& nodes) const {
    std::vector<int> master, coarse, potentials;
    nodes.neighbor().setMasterList(position, H, master, coarse);
    nodes.neighbor().setRefineNeighborList(position, H, coarse, potentials);
    std::vector<int> result;
    for (const auto i: potentials) {
      const auto displacement = nodes.positions()[i] - position;
      if (std::min((H*displacement).magnitude(),
                   (nodes.Hfield()[i]*displacement).magnitude()) <= kernelExtent) {
        result.push_back(int(i));
      }
    }
    return result;
  }

  std::vector<int>
  findGatherNeighborNodes(const Vector& position,
                          const SymTensor& H,
                          const NodeList& nodes) const {
    std::vector<int> master, coarse, potentials;
    nodes.neighbor().setMasterList(position, H, master, coarse);
    nodes.neighbor().setRefineNeighborList(position, H, coarse, potentials);
    std::vector<int> result;
    for (const auto i: potentials) {
      if ((H*(position - nodes.positions()[i])).magnitude() <= kernelExtent) {
        result.push_back(int(i));
      }
    }
    return result;
  }

  std::vector<int>
  findScatterNeighborNodes(const Vector& position,
                           const NodeList& nodes) const {
    std::vector<int> master, coarse, potentials;
    nodes.neighbor().setMasterList(position, master, coarse);
    nodes.neighbor().setRefineNeighborList(position, coarse, potentials);
    std::vector<int> result;
    for (const auto i: potentials) {
      if ((nodes.Hfield()[i]*(position - nodes.positions()[i])).magnitude() <= kernelExtent) {
        result.push_back(int(i));
      }
    }
    return result;
  }

  std::vector<std::vector<int>>
  findOverlapNeighbors(const Vector& position,
                       const SymTensor& H) const {
    const auto& nodeLists = dataBase->nodeListPtrs();
    std::vector<std::set<int>> uniqueResult(nodeLists.size());
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
    std::vector<std::vector<int>> result(nodeLists.size());
    for (size_t nodeListID = 0; nodeListID < nodeLists.size(); ++nodeListID) {
      result[nodeListID].assign(uniqueResult[nodeListID].begin(),
                                uniqueResult[nodeListID].end());
    }
    return result;
  }

  static bool
  checkNeighbors(const std::vector<int>& neighborList,
                 const std::vector<int>& answer) {
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

  static std::vector<int>
  sorted(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    return values;
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

      // Randomly select nodes from each NodeList to explicitly test.
      for (const auto nodeID: randomSample(nodes->numInternalNodes() - 1u, ncheck)) {
        const auto ri = position[nodeID];
        const auto Hi = H[nodeID];

        // Have the neighbor objects select neighbors for this node.
        std::vector<std::vector<int>> masterLists, coarseNeighbors, refineNeighbors;
        dataBase->setMasterNodeLists(ri, Hi, masterLists, coarseNeighbors, false);
        dataBase->setRefineNodeLists(ri, Hi, coarseNeighbors, refineNeighbors);
        std::vector<int> neighborIDs;
        auto offset = 0u;
        const auto& nodeLists = dataBase->nodeListPtrs();
        for (size_t nodeListID = 0; nodeListID < nodeLists.size(); ++nodeListID) {
          for (const auto i: refineNeighbors[nodeListID]) {
            neighborIDs.push_back(i + int(offset));
          }
          offset += nodeLists[nodeListID]->numInternalNodes();
        }

        // Now build the checks.
        std::vector<int> answerIDs;
        offset = 0u;
        for (const auto* answerNodes: nodeLists) {
          for (const auto i: findNeighborNodes(ri, Hi, *answerNodes)) {
            answerIDs.push_back(i + int(offset));
          }
          offset += answerNodes->numInternalNodes();
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

      // Randomly select nodes from each NodeList to explicitly test.
      for (const auto i: randomSample(inodes->numInternalNodes() - 1u, ncheck)) {
        const auto ri = position[i];
        const auto Hi = H[i];
        const auto& cmneighbors = cm.connectivityForNode(inodes, int(i));

        // Check ConnectivityMap vs. N^2 neighbor search.
        for (size_t jNL = 0; jNL < nodeLists.size(); ++jNL) {
          const auto cmcheck = sorted(cmneighbors[jNL]);
          auto answer = sorted(findNeighborNodes(ri, Hi, *nodeLists[jNL]));
          if (iNL == jNL) {
            const auto self = std::find(answer.begin(), answer.end(), int(i));
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
        const auto& cmneighbors = cm.connectivityForNode(inodes, int(i));
        ASSERT_EQ(cmneighbors.size(), numNodeLists);
        for (size_t jNL = iNL; jNL < numNodeLists; ++jNL) {
          for (const auto j: cmneighbors[jNL]) {
            if (jNL > iNL or j > int(i)) answer.emplace_back(i, iNL, j, jNL);
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

      // Randomly select nodes from each NodeList to explicitly test.
      for (const auto i: randomSample(inodes->numInternalNodes() - 1u,
                                      noverlapcheck)) {
        const auto ri = inodes->positions()[i];
        const auto Hi = inodes->Hfield()[i];
        const auto& cmneighbors = cm.overlapConnectivityForNode(inodes, int(i));
        const auto fullanswer = findOverlapNeighbors(ri, Hi);

        // Check ConnectivityMap vs. N^2 neighbor search.
        for (size_t jNL = 0; jNL < nodeLists.size(); ++jNL) {
          const auto cmcheck = sorted(cmneighbors[jNL]);
          auto answer = sorted(fullanswer[jNL]);
          if (iNL == jNL) {
            const auto self = std::find(answer.begin(), answer.end(), int(i));
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
      std::vector<std::set<int>> result(numNodeLists);
      for (size_t kNL = 0; kNL < numNodeLists; ++kNL) {
        const std::set<int> seti(neighborsi[kNL].begin(), neighborsi[kNL].end());
        const std::set<int> setj(neighborsj[kNL].begin(), neighborsj[kNL].end());
        std::set_intersection(seti.begin(), seti.end(),
                              setj.begin(), setj.end(),
                              std::inserter(result[kNL], result[kNL].begin()));
      }
      result[pair.i_list].insert(int(pair.i_node));
      result[pair.j_list].insert(int(pair.j_node));
      return result;
    };

    // Make vector<vector<int>> into a vector of sets for comparison.
    const auto convertToSets = [&](const auto& neighbors) {
      EXPECT_EQ(neighbors.size(), numNodeLists);
      std::vector<std::set<int>> result;
      result.reserve(neighbors.size());
      for (const auto& values: neighbors) result.emplace_back(values.begin(), values.end());
      return result;
    };

    // Randomly select pairs to test.
    for (const auto pairIndex: randomSample(pairs.size(), ncheck)) {
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
    const auto n1 = 1000u;
    const auto n2 = 2500u;
    const auto n3 = 500u;
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
    const auto n1 = 10000u;
    const auto n2 = 25000u;
    const auto n3 = 5000u;
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
    noverlapcheck = 0u;
    const auto n1 = 1000u;
    const auto n2 = 2500u;
    const auto n3 = 1500u;
    const VolumeRange volrange1 = {{{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}};
    const VolumeRange volrange2 = {{{1.0, 1.5}, {0.0, 1.0}, {0.0, 1.0}}};
    const VolumeRange volrange3 = {{{1.5, 2.0}, {0.0, 1.0}, {0.0, 1.0}}};
    genericSetUp(n1, n2, n3, volrange1, volrange2, volrange3);
  }
};

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
    ncheck = 50u;
    noverlapcheck = 2u;

    // Match GenerateNodeDistribution2d(..., distributionType="constantDTheta").
    eos = std::make_unique<GammaLawGas>(2.0, 2.0, units,
                                       -1.0e200, 1.0e200,
                                       MaterialPressureMinType::PressureFloor,
                                       0.0);
    WT = std::make_unique<TableKernel>(BSplineKernel<Dim<2>>(), 100u);
    nodes1 = std::make_unique<NodeSet>("cylindrical nodes 1", *eos,
                                      Vector(-10.0), Vector(10.0));

    constexpr auto nRadial = 100u;
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

    dataBase = std::make_unique<DataBase>();
    dataBase->appendNodeList(nodes);
    nodes.neighbor().updateNodes();
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
