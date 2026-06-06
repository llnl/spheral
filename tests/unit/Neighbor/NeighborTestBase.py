from math import *
import unittest
import random

from Spheral import *

#===============================================================================
# Apply a random rotation operation to the given symmetric tensor.
#===============================================================================
def applyRandomRotation(t, nDimensions):
    if nDimensions == 2:
        theta = random.uniform(0.0, pi)
        R = Tensor2d(cos(theta), -sin(theta),
                     sin(theta), cos(theta))
        t.rotationalTransform(R)
    elif nDimensions == 3:
        theta1 = random.uniform(0.0, pi)
        theta2 = random.uniform(0.0, pi)
        theta3 = random.uniform(0.0, pi)
        R1 = Tensor3d(cos(theta1), -sin(theta1),  0.0,
                      sin(theta1),  cos(theta1),  0.0,
                      0.0,          0.0,          1.0)
        R2 = Tensor3d(cos(theta2),  0.0,        -sin(theta2),
                      0.0,          1.0,          0.0,
                      sin(theta2),  0.0,          cos(theta2))
        R3 = Tensor3d(1.0,          0.0,          0.0,
                      0.0,          cos(theta3), -sin(theta3),
                      0.0,          sin(theta3),  cos(theta3))
        R = R1*R2*R3
        t.rotationalTransform(R)
    return

#===============================================================================
# Base class for the tests, providing the generic set up methods.
#===============================================================================
class NeighborTestBase:

    #---------------------------------------------------------------------------
    # Distribute nodes randomly amongst domains.
    #---------------------------------------------------------------------------
    def randomDistribute(self,
                         n,                  # global number of nodes
                         volrange,              # total simulation volume
                         Vector,             # the vector type
                         SymTensor,          # tensor type
                         nNodesPerh = 2.01): # duh!
        assert Vector.nDimensions == SymTensor.nDimensions
        assert len(volrange) == Vector.nDimensions
        assert min([len(volrange[i]) == 2 for i in range(Vector.nDimensions)])

        adim = 1.0/SymTensor.nDimensions
        vol = 1.0
        for i in range(SymTensor.nDimensions):
            vol *= volrange[i][1] - volrange[i][0]
        assert vol > 0.0
        dx0 = (vol/n)**adim

        globalNodeIDs = []
        nodePositions = []
        H = []
        for i in range(n):
            globalNodeIDs.append(i)
            args = (random.uniform(volrange[j][0], volrange[j][1])
                    for j in range(Vector.nDimensions))
            nodePositions.append(Vector(*args))

            dx = [nNodesPerh * random.uniform(0.5, 2.0)*dx0
                  for j in range(SymTensor.nDimensions)]
            assert len(dx) == SymTensor.nDimensions
            assert min(dx) > 0.0
            args = [0.0] * SymTensor.nDimensions**2
            for j in range(SymTensor.nDimensions):
                args[j + j*SymTensor.nDimensions] = 1.0/dx[j]
            Hi = SymTensor(*tuple(args))
            applyRandomRotation(Hi, SymTensor.nDimensions)
            H.append(Hi)

        assert len(globalNodeIDs) == n
        assert len(nodePositions) == n
        assert len(H) == n
        return globalNodeIDs, nodePositions, H

    #---------------------------------------------------------------------------
    # setUp routine common to all instances.
    #---------------------------------------------------------------------------
    def genericSetUp(self, n1, n2, n3,
                     volrange1, volrange2, volrange3,
                     EOS,
                     makeFluidNodeList,
                     TableKernel,
                     BSplineKernel,
                     Vector,
                     SymTensor,
                     NeighborType,
                     DataBase,
                     Seeder):
        
        # Construct the NodeLists to be distributed
        self.eos = EOS(2.0, 2.0)
        self.WT = TableKernel(BSplineKernel(), 100)
        self.nodes1 = makeFluidNodeList("nodes 1", self.eos, NeighborType=NeighborType)
        self.nodes2 = makeFluidNodeList("nodes 2", self.eos, NeighborType=NeighborType)
        self.nodes3 = makeFluidNodeList("nodes 3", self.eos, NeighborType=NeighborType)
        for nodes, nGlobal, volrange in ((self.nodes1, n1, volrange1),
                                      (self.nodes2, n2, volrange2),
                                      (self.nodes3, n3, volrange3)):
            globalIDs, xyNodes, H = Seeder(nGlobal,
                                           volrange,
                                           Vector,
                                           SymTensor)
            n = len(globalIDs)
            nodes.numInternalNodes = n
            for i in range(n):
                nodes.mass()[i] = 1.0
                nodes.positions()[i] = xyNodes[i]
                nodes.Hfield()[i] = H[i]

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase()
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)

        # Update all the Neighbor info.
        #NeighborType.setBoundingBox()
        for n in (self.nodes1, self.nodes2, self.nodes3):
            n.neighbor().updateNodes()

        return

    #---------------------------------------------------------------------------
    # Test raw Neighbor object calls.
    #---------------------------------------------------------------------------
    def testNeighborCalls(self):
        from SpheralTestUtilities import findNeighborNodes, checkNeighbors
        import time

        # Iterate over the NodeLists.
        for nodes in self.dataBase.nodeLists:
            pos = nodes.positions()
            H = nodes.Hfield()

            # Randomly select nodes from each NodeList to explicitly test.
            for nodeID in random.sample(list(range(nodes.numInternalNodes - 1)), self.ncheck):
                ri = pos[nodeID]
                Hi = H[nodeID]

                # Have the neighbor objects select neighbors for this node.
                t0 = time.time()
                masterLists = vector_of_vector_of_int()
                coarseNeighbors = vector_of_vector_of_int()
                refineNeighbors = vector_of_vector_of_int()
                self.dataBase.setMasterNodeLists(ri, Hi, masterLists, coarseNeighbors, False)
                self.dataBase.setRefineNodeLists(ri, Hi, coarseNeighbors, refineNeighbors)
                neighborIDs = []
                offset = 0
                for inds, nds in enumerate(self.dataBase.nodeLists):
                    neighborIDs.extend([i + offset for i in refineNeighbors[inds]])
                    offset += nds.numInternalNodes
                t1 = time.time()

                # Now build the checks.
                answerIDs = []
                offset = 0
                for nds in self.dataBase.nodeLists:
                    answerIDs.extend([i + offset for i in findNeighborNodes(ri, Hi, self.kernelExtent, nds)])
                    offset += nds.numInternalNodes
                t2 = time.time()

                # Check the answer.
                test = checkNeighbors(neighborIDs, answerIDs)
                t3 = time.time()
                if not test:
                    neighborIDs.sort()
                    answerIDs.sort()
                    print('SPH Neighbor test FAILED')
                    print(' refine: ', neighborIDs)
                    print(' answer: ', answerIDs)
                    missing = [i for i in answerIDs if i not in neighborIDs]
                    print('missing: ', missing)
                    print('deltas: ', [((Hi*(pos[i] - ri)).x, (H[i]*(pos[i] - ri)).x) for i in missing])
                else:
                    print("Passed for node %i : %f %f %f" % (nodeID, t1 - t0, t2 - t1, t3 - t2))
                assert test

    #---------------------------------------------------------------------------
    # Test ConnectivityMap neighbors
    #---------------------------------------------------------------------------
    def testConnectivityMapNeighbors(self):
        from SpheralTestUtilities import findNeighborNodes, checkNeighbors
        import time

        self.dataBase.updateConnectivityMap(False, False, False)
        cm = self.dataBase.connectivityMap(False, False, False)

        # Iterate over the NodeLists.
        for iNL, inodes in enumerate(self.dataBase.nodeLists):
            pos = inodes.positions()
            H = inodes.Hfield()

            # Randomly select nodes from each NodeList to explicitly test.
            for i in random.sample(list(range(inodes.numInternalNodes - 1)), self.ncheck):
                ri = pos[i]
                Hi = H[i]
                cmneighbors = cm.connectivityForNode(inodes, i)

                # Check ConnectivityMap vs. N^2 neighbor search.
                for jNL, jnodes in enumerate(self.dataBase.nodeLists):
                    cmcheck = sorted(cmneighbors[jNL])
                    answer = sorted(findNeighborNodes(ri, Hi, self.kernelExtent, jnodes))
                    if iNL == jNL:
                        answer.remove(i)
                    if not cmcheck == answer:
                        print('SPH ConnectivityMap neighbor test FAILED for node %i' % i)
                        print('     CM: ', cmcheck)
                        print(' answer: ', answer)
                        missing = [i for i in answer if i not in cmcheck]
                        print('missing: ', missing)
                        print('deltas: ', [((Hi*(pos[i] - ri)).x, (H[i]*(pos[i] - ri)).x) for i in missing])
                        raise RuntimeError("Failed test")
                    else:
                        print("Passed for node %i" % i)

    #---------------------------------------------------------------------------
    # Test ConnectivityMap NodePairList
    #---------------------------------------------------------------------------
    def testConnectivityMapNodePairs(self):
        from SpheralTestUtilities import findNeighborNodes, checkNeighbors
        import time

        self.dataBase.updateConnectivityMap(False, False, False)
        cm = self.dataBase.connectivityMap(False, False, False)
        numNodeLists = self.dataBase.numNodeLists

        # Build the answer based on the node neighbors
        answer = []
        for iNL, inodes in enumerate(self.dataBase.nodeLists):
            for i in range(inodes.numInternalNodes):
                cmneighbors = cm.connectivityForNode(inodes, i)
                assert len(cmneighbors) == numNodeLists
                for jNL in range(iNL, numNodeLists):
                    for j in cmneighbors[jNL]:
                        if jNL > iNL or j > i:
                            answer.append(NodePairIdxType(i, iNL, j, jNL))

        pairs = list(cm.nodePairList)
        assert len(pairs) == len(answer)

        answer.sort()
        pairs.sort()

        # for i in xrange(len(pairs)):
        #     print "%s\t\t%s" % (pairs[i], answer[i])

        assert answer == pairs

        # for pair in answer:
        #     if not (pair in answer or (pair[2], pair[3], pair[0], pair[1]) in answer):
        #         raise RuntimeError, "Missing pair interaction for %s" % pair
        return

    #---------------------------------------------------------------------------
    # Test ConnectivityMap overlap neighbors
    #---------------------------------------------------------------------------
    def testConnectivityMapOverlapNeighbors(self):
        from SpheralTestUtilities import findNeighborNodes, findOverlapNeighbors, findOverlapRegion, checkNeighbors
        import time

        self.dataBase.updateConnectivityMap(False, True, False)
        cm = self.dataBase.connectivityMap(False, True, False)
        pos = self.dataBase.globalPosition
        H = self.dataBase.globalHfield

        # Local extra info method
        def etaNeighborStats(iNL, i, jNL, j, actual):
            ri = pos(iNL, i)
            Hi = H(iNL, i)
            rj = pos(jNL, j)
            Hj = H(jNL, j)
            for (kNL, nghbs) in enumerate(actual):
                for k in nghbs:
                    rk = pos(kNL, k)
                    print("    --> ", k, (Hi*(ri - rk)).magnitude(), (Hj*(rj - rk)).magnitude())

        def actualIntersection(iNL, i, jNL, j):
            actual = []
            ri = pos(iNL, i)
            Hi = H(iNL, i)
            rj = pos(jNL, j)
            Hj = H(jNL, j)
            for nodes in self.dataBase.nodeLists:
                actual.append(findOverlapRegion(ri, Hi, rj, Hj, self.kernelExtent, nodes))
            return actual

        # Iterate over the NodeLists.
        for iNL, inodes in enumerate(self.dataBase.nodeLists):

            # Randomly select nodes from each NodeList to explicitly test.
            for i in random.sample(list(range(inodes.numInternalNodes - 1)), self.noverlapcheck):
                print("Checking ", i)
                ri = pos(iNL, i)
                Hi = H(iNL, i)
                cmneighbors = cm.overlapConnectivityForNode(inodes, i)
                fullanswer = findOverlapNeighbors(ri, Hi, self.kernelExtent, self.dataBase)

                # Check ConnectivityMap vs. N^2 neighbor search.
                for jNL, jnodes in enumerate(self.dataBase.nodeLists):
                    cmcheck = sorted(cmneighbors[jNL])
                    answer = sorted(fullanswer[jNL])
                    if iNL == jNL:
                        answer.remove(i)
                    if not cmcheck == answer:
                        print('SPH ConnectivityMap overlap neighbor test FAILED for node %i' % i)
                        print('     CM: ', cmcheck)
                        print(' answer: ', answer)
                        missing = [j for j in answer if not j in cmcheck]
                        print('missing: ', missing)
                        print('intersections for missing:')
                        for j in missing:
                            rj = pos(jNL, j)
                            Hj = H(jNL, j)
                            print('   %i: ' % j, [list(x) for x in cm.connectivityIntersectionForNodes(iNL, i, jNL, j)])
                            actual = actualIntersection(iNL, i, jNL, j)
                            print('   %i: ' % j, actual)
                            etaNeighborStats(iNL, i, jNL, j, actual)
                            print("Finished")
                        extra = [j for j in cmcheck if not j in answer]
                        print('extra: ', extra)
                        print('intersections for extra:')
                        for j in extra:
                            print('   %i: ' % j, [list(x) for x in cm.connectivityIntersectionForNodes(iNL, i, jNL, j)])
                            actual = actualIntersection(iNL, i, jNL, j)
                            print('   %i: ' % j, actual)
                            etaNeighborStats(iNL, i, jNL, j, actual)
                            print("Finished")
                        print("REALLY Finished")
                        sys.stdin.flush()
                        raise RuntimeError("Failed test")
                print("    Passed for node %i" % i)

        return

    #---------------------------------------------------------------------------
    # Test ConnectivityMap intersection methods
    #---------------------------------------------------------------------------
    def testConnectivityComputeIntersection(self):
        from SpheralTestUtilities import findNeighborNodes, findOverlapNeighbors, findOverlapRegion, checkNeighbors
        import time

        self.dataBase.updateConnectivityMap(False, True, True)
        numNodeLists = self.dataBase.numNodeLists
        cm = self.dataBase.connectivityMap(False, True, True)
        pairs = list(cm.nodePairList)
        pos = self.dataBase.globalPosition
        H = self.dataBase.globalHfield

        # Compute the actual intersection of two points for a NodeList
        # Returns a list of Python sets (per NodeList)
        def intersectionAnswer(pair):
            neighborsi = cm.connectivityForNode(pair.i_list, pair.i_node)
            neighborsj = cm.connectivityForNode(pair.j_list, pair.j_node)
            result = [set(neighborsi[kNL]) & set(neighborsj[kNL]) for kNL in range(numNodeLists)]
            result[pair.i_list].add(pair.i_node)
            result[pair.j_list].add(pair.j_node)
            return result

        # Make our vector<vector<int>> into a list of sets for comparison with the answer
        def convertToSets(neighbors):
            assert len(neighbors) == numNodeLists
            return [set(x) for x in neighbors]

        # Randomly select pairs to test
        for pair in random.sample(pairs, self.ncheck):
            print("Checking ", pair)
            answer = intersectionAnswer(pair)
            intersect = convertToSets(answer)
            assert len(intersect) == numNodeLists
            print(answer)
            print(intersect)
            self.assertEqual(intersect, answer,
                            ("\nIntersection computed from Neighbor set intersections does not match: \nintersect: %s\nanswer: %s\n" % (intersect, answer)))
            print("    Passed for pair ", pair)

        return

#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
class NeighborRandom1d(unittest.TestCase, NeighborTestBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print("--------------------------------------------------------------------------------")
        print("1-D %s random test." % NeighborRandom1d._NeighborType.__name__)
        print("--------------------------------------------------------------------------------")

        self.ncheck = 10
        self.noverlapcheck = 2

        # Generic parameters for 1-D tests.
        n1 = 1000
        n2 = 2500
        n3 = 500

        volrange1 = ((-2.0, -1.0),)
        volrange2 = ((-1.0, 0.5),)
        volrange3 = ((0.5, 2.0),)

        searchType = GatherScatter
        numGridLevels = 20
        topGridCellSize = 100.0
        self.kernelExtent = 2.0

        self.genericSetUp(n1, n2, n3,
                          volrange1, volrange2, volrange3,
                          GammaLawGasMKS1d,
                          makeFluidNodeList1d,
                          TableKernel1d,
                          BSplineKernel1d,
                          Vector1d,
                          SymTensor1d,
                          NeighborRandom1d._NeighborType,
                          DataBase1d,
                          self.randomDistribute)

        return

    #---------------------------------------------------------------------------
    # Method called after test is completed.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.dataBase, self.nodes1, self.nodes2, self.nodes3
        return

#===============================================================================
# Radom node distribution -- 2-D.
#===============================================================================
class NeighborRandom2d(unittest.TestCase, NeighborTestBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print("--------------------------------------------------------------------------------")
        print("2-D %s random test." % NeighborRandom2d._NeighborType.__name__)
        print("--------------------------------------------------------------------------------")

        self.ncheck = 10
        self.noverlapcheck = 2

        # Generic parameters for 2-D tests.
        n1 = 10000
        n2 = 25000
        n3 = 5000

        volrange1 = ((-2.0, -1.0), (0.0, 1.0))
        volrange2 = ((-1.0, -0.5), (0.0, 1.0))
        volrange3 = ((-0.5, 0.0), (0.0, 1.0))

        self.kernelExtent = 2.0

        self.genericSetUp(n1, n2, n3,
                          volrange1, volrange2, volrange3,
                          GammaLawGasMKS2d,
                          makeFluidNodeList2d,
                          TableKernel2d,
                          BSplineKernel2d,
                          Vector2d,
                          SymTensor2d,
                          NeighborRandom2d._NeighborType,
                          DataBase2d,
                          self.randomDistribute)

        return

    #---------------------------------------------------------------------------
    # Method called after test is completed.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.dataBase, self.nodes1, self.nodes2, self.nodes3
        return

#===============================================================================
# Radom node distribution -- 3-D.
#===============================================================================
class NeighborRandom3d(unittest.TestCase, NeighborTestBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print("--------------------------------------------------------------------------------")
        print("3-D %s random test." % NeighborRandom3d._NeighborType.__name__)
        print("--------------------------------------------------------------------------------")

        self.ncheck = 10
        self.noverlapcheck = 0

        # Generic parameters for 3-D tests.
        n1 = 1000
        n2 = 2500
        n3 = 1500

        volrange1 = ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
        volrange2 = ((1.0, 1.5), (0.0, 1.0), (0.0, 1.0))
        volrange3 = ((1.5, 2.0), (0.0, 1.0), (0.0, 1.0))

        searchType = GatherScatter
        numGridLevels = 20
        topGridCellSize = 100.0
        origin = Vector3d(0.0, 0.0, 0.0)
        self.kernelExtent = 2.0

        self.genericSetUp(n1, n2, n3,
                          volrange1, volrange2, volrange3,
                          GammaLawGasMKS3d,
                          makeFluidNodeList3d,
                          TableKernel3d,
                          BSplineKernel3d,
                          Vector3d,
                          SymTensor3d,
                          NeighborRandom3d._NeighborType,
                          DataBase3d,
                          self.randomDistribute)

        return

    #---------------------------------------------------------------------------
    # Method called after test is completed.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.dataBase, self.nodes1, self.nodes2, self.nodes3
        return

#===============================================================================
# Cylindrical node distribution -- 2-D.
#===============================================================================
class NeighborCylindrical2d(unittest.TestCase, NeighborTestBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print("--------------------------------------------------------------------------------")
        print("2-D %s regular cylindrical test." % self._NeighborType)
        print("--------------------------------------------------------------------------------")

        self.ncheck = 50
        self.noverlapcheck = 2

        from GenerateNodeDistribution2d import GenerateNodeDistribution2d
        from DistributeNodes import distributeNodes2d
        self.eos = GammaLawGasMKS2d(2.0, 2.0)
        self.WT = TableKernel2d(BSplineKernel2d(), 100)
        self.nodes1 = makeFluidNodeList2d("cylindrical nodes 1", self.eos, NeighborType=NeighborCylindrical2d._NeighborType)
        self.kernelExtent = 2.0
        generator = GenerateNodeDistribution2d(nRadial = 100,
                                               nTheta = 100,
                                               rho = 1.0,
                                               distributionType = "constantDTheta",
                                               rmin = 0.0,
                                               rmax = 1.0,
                                               theta = 0.5*pi,
                                               nNodePerh = 2.01)
        distributeNodes2d((self.nodes1, generator))
        self.dataBase = DataBase2d()
        self.dataBase.appendNodeList(self.nodes1)

        return

    #---------------------------------------------------------------------------
    # Method called after test is completed.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.dataBase, self.nodes1
        return

