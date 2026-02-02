#ATS:test(SELF, np=8, label="PeanoHilbertOrderRedistribute unit tests")
from math import *
import unittest

from Spheral import *
from testDistributeByPosition2d import TestDistributeByPosition2d

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
nDomains = mpi.procs

#===============================================================================
# Main testing class for the 2-D tests.
#===============================================================================
class TestPeanoHilbertOrderRedistribute2d(TestDistributeByPosition2d):

    # The actual test itself!
    # Create a PeanoHilbertOrderRedistributeNodes object, have it redistribute the
    # nodes.
    def testIt(self):
        print("Testing PeanoHilbertOrderRedistributeNodes2d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

##         # Plot the initial conditions.
##         plot0 = plotNodePositions2d(self.dataBase,
##                                     colorNodeLists = False,
##                                     colorDomains = True,
##                                     title = "Initial conditions",
##                                     style = "linespoints",
##                                     persist = True)

        # Go ahead and redistribute those nodes!
        repartition = PeanoHilbertOrderRedistributeNodes2d(2.0)
        repartition.redistributeNodes(self.dataBase)

##         # Plot the final conditions.
##         plot1 = plotNodePositions2d(self.dataBase,
##                                     colorNodeLists = False,
##                                     colorDomains = True,
##                                     title = "After redistribution",
##                                     style = "linespoints",
##                                     persist = True)

        # Make sure that the numbers of nodes are correct.
        assert self.dataBase.numNodeLists == len(nNodesGlobal)
        i = 0
        for nodeList in self.dataBase.nodeLists:
            n = mpi.allreduce(nodeList.numInternalNodes, mpi.SUM)
            nGlobal = nNodesGlobal[i]
            if n != nGlobal:
                self.fail("Wrong number of nodes: %i != %i" % (n, nGlobal))
            i += 1

        # Do a consistency check of the new distribution.
        globalIDs = globalNodeIDsAll2d(self.dataBase)
        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase,
                                                                  globalIDs)
        if (not repartition.validDomainDecomposition(nodeDistribution,
                                                     self.dataBase)):
            self.fail("Invalid domain decomposition.")

##        import SpheralGnuPlotUtilities
##        self.p = SpheralGnuPlotUtilities.plotNodePositions2d(self.dataBase,
##                                                        colorNodeLists = 0,
##                                                        colorDomains = 1)

#===============================================================================
# Distribute nodes randomly amongst domains.
#===============================================================================
def randomDistribute3d(thisDomain,     # domain ID to be calculated
                       nDomains,       # total number of domains
                       nNodesGlobal,   # global number of nodes in this nodelist
                       xyzRangeTotal): # total simulation volume

    assert thisDomain >= 0 and thisDomain < nDomains
    assert nDomains > 0

    import random
    g = random.Random()
    globalNodeIDs = []
    nodePositions = []
    for globalNodeID in range(0, nNodesGlobal):
        mpi.barrier()
        domain0 = g.randint(0, nDomains - 1)
        domain = mpi.bcast(domain0)
        if domain == thisDomain:
            globalNodeIDs.append(globalNodeID)
            nodePositions.append(Vector3d(g.uniform(xyzRangeTotal[0][0],
                                                    xyzRangeTotal[1][0]),
                                          g.uniform(xyzRangeTotal[0][1],
                                                    xyzRangeTotal[1][1]),
                                          g.uniform(xyzRangeTotal[0][2],
                                                    xyzRangeTotal[1][2])))

    assert len(globalNodeIDs) == len(nodePositions)
    assert mpi.allreduce(len(globalNodeIDs), mpi.SUM) == nNodesGlobal
    return globalNodeIDs, nodePositions

#===============================================================================
# Calculate one over the smoothing scale for the given number of nodes and
# volume.
#===============================================================================
def determineH3d(nGlobal, xRange,
                 nNodesPerh = 2.01):
    assert nGlobal > 0
    
    vol = ((xRange[1][0] - xRange[0][0])*
           (xRange[1][1] - xRange[0][1])*
           (xRange[1][2] - xRange[0][2]))
    assert vol > 0.0
    dV = vol/nGlobal
    dx = dV**(1.0/3.0)
    hi = 1.0/(nNodesPerh*dx)
    Hi = SymTensor3d(hi, 0.0, 0.0,
                     0.0, hi, 0.0,
                     0.0, 0.0, hi)
    return Hi

#===============================================================================
# Main testing class for the 3-D tests.
#===============================================================================
class TestPeanoHilbertOrderRedistribute3d(unittest.TestCase):

    # Set up method called before test is run.
    def setUp(self):

        # Generic parameters for 2-D tests.
        n1 = 1000
        n2 = 2000
        n3 = 500

        range1 = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]
        range2 = [(1.0, 0.0, 0.0), (1.5, 1.0, 1.0)]
        range3 = [(1.5, 0.0, 0.0), (2.0, 1.0, 1.0)]

        # Construct the NodeLists to be distributed
        self.eos = GammaLawGasMKS3d(2.0, 2.0)
        self.WT = TableKernel3d(BSplineKernel3d())
        self.nodes1 = makeFluidNodeList3d("nodes1 3d", self.eos)
        self.nodes2 = makeFluidNodeList3d("nodes2 3d", self.eos)
        self.nodes3 = makeFluidNodeList3d("nodes3 3d", self.eos)
        for (nodes, nGlobal, globalRange) in ((self.nodes1, n1, range1),
                                              (self.nodes2, n2, range2),
                                              (self.nodes3, n3, range3)):
            globalIDs, xyzNodes = randomDistribute3d(domainID,
                                                     nDomains,
                                                     nGlobal,
                                                     globalRange)
            n = len(globalIDs)
            nodes.numInternalNodes = n
            Hi = determineH3d(nGlobal, globalRange)
            for i in range(n):
                nodes.mass()[i] = 1.0
                nodes.positions()[i] = xyzNodes[i]
                nodes.Hfield()[i] = Hi

            nodes.neighbor().updateNodes()

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase3d()
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)

        assert mpi.allreduce(self.nodes1.numInternalNodes, mpi.SUM) == n1
        assert mpi.allreduce(self.nodes2.numInternalNodes, mpi.SUM) == n2
        assert mpi.allreduce(self.nodes3.numInternalNodes, mpi.SUM) == n3

        return

    # Method called after test is completed.
    def tearDown(self):
        del self.nodes1, self.nodes2, self.nodes3
        return

    # The actual test itself!
    # Create a PeanoHilbertOrderRedistributeNodes object, have it redistribute the
    # nodes.
    def testIt(self):
        print("Testing PeanoHilbertOrderRedistributeNodes3d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

        # Go ahead and redistribute those nodes!
        repartition = PeanoHilbertOrderRedistributeNodes3d(2.0)
        repartition.redistributeNodes(self.dataBase)

        # Make sure that the numbers of nodes are correct.
        assert self.dataBase.numNodeLists == len(nNodesGlobal)
        i = 0
        for nodeList in self.dataBase.nodeLists:
            n = mpi.allreduce(nodeList.numInternalNodes, mpi.SUM)
            nGlobal = nNodesGlobal[i]
            if n != nGlobal:
                self.fail("Wrong number of nodes: %i != %i" % (n, nGlobal))
            i += 1

        # Do a consistency check of the new distribution.
        globalIDs = globalNodeIDsAll3d(self.dataBase)
        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase,
                                                                  globalIDs)
        if (not repartition.validDomainDecomposition(nodeDistribution,
                                                     self.dataBase)):
            self.fail("Invalid domain decomposition.")

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
