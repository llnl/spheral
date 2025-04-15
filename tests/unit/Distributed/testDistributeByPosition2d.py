#ATS:test(SELF, np=8, label="DistributeByPosition2d unit tests")
from math import *
import unittest

from Spheral import *

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
nDomains = mpi.procs

#===============================================================================
# Distribute nodes randomly amongst domains.
#===============================================================================
def randomDistribute(thisDomain,     # domain ID to be calculated
                     nDomains,       # total number of domains
                     nNodesGlobal,   # global number of nodes in this nodelist
                     xyRangeTotal):  # total simulation volume

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
            nodePositions.append(Vector2d(g.uniform(xyRangeTotal[0][0],
                                                    xyRangeTotal[1][0]),
                                          g.uniform(xyRangeTotal[0][1],
                                                    xyRangeTotal[1][1])))

    assert len(globalNodeIDs) == len(nodePositions)
    assert mpi.allreduce(len(globalNodeIDs), mpi.SUM) == nNodesGlobal
    return globalNodeIDs, nodePositions

#===============================================================================
# Calculate one over the smoothing scale for the given number of nodes and
# volume.
#===============================================================================
def determineH(nGlobal, xRange,
               nNodesPerh = 2.01):
    assert nGlobal > 0
    
    vol = (xRange[1][0] - xRange[0][0])*(xRange[1][1] - xRange[0][1])
    assert vol > 0.0
    dV = vol/nGlobal
    dx = sqrt(dV)
    hi = 1.0/(nNodesPerh*dx)
    Hi = SymTensor2d(hi, 0.0,
                     0.0, hi)
    return Hi

#===============================================================================
# Main testing class for the 2-D tests.
#===============================================================================
class TestDistributeByPosition2d(unittest.TestCase):

    # Set up method called before test is run.
    def setUp(self):

        # Generic parameters for 2-D tests.
        n1 = 1000
        n2 = 2500
        n3 = 500

        range1 = [(0.0, 0.0), (1.0, 1.0)]
        range2 = [(1.0, 0.0), (1.5, 1.0)]
        range3 = [(1.5, 0.0), (2.0, 1.0)]

        # Construct the NodeLists to be distributed
        self.eos = GammaLawGasMKS2d(2.0, 2.0)
        self.WT = TableKernel2d(BSplineKernel2d())
        self.nodes1 = makeFluidNodeList2d("nodes1 2d", self.eos)
        self.nodes2 = makeFluidNodeList2d("nodes2 2d", self.eos)
        self.nodes3 = makeFluidNodeList2d("nodes3 2d", self.eos)
        for (nodes, nGlobal, globalRange) in ((self.nodes1, n1, range1),
                                              (self.nodes2, n2, range2),
                                              (self.nodes3, n3, range3)):
            globalIDs, xyNodes = randomDistribute(domainID,
                                                 nDomains,
                                                 nGlobal,
                                                 globalRange)
            n = len(globalIDs)
            nodes.numInternalNodes = n
            Hi = determineH(nGlobal, globalRange)
            for i in range(n):
                nodes.mass()[i] = 1.0
                nodes.positions()[i] = xyNodes[i]
                nodes.Hfield()[i] = Hi

            nodes.neighbor().updateNodes()

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase2d()
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
    # Create a DistributeNodeByXPosition object, have it redistribute the
    # nodes, and check that each domains x distribution is contiguous and
    # in the right order.
    def testIt(self):
        print("Testing DistributeByXPosition2d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

        # Go ahead and redistribute those nodes!
        repartition = DistributeByXPosition2d()
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

        # Have each domain figure out it's min and max x.
        localxmin = 1e10
        localxmax = -1e10
        for nodeList in self.dataBase.nodeLists:
            if nodeList.numInternalNodes > 0:
                localxmin = min(localxmin, min([nodeList.positions()[i].x for i in range(nodeList.numNodes)]))
                localxmax = max(localxmax, max([nodeList.positions()[i].x for i in range(nodeList.numNodes)]))

        # Now make sure that our (xmin,xmax) range is greater than any processor
        # less than us, and less than any processor greater.
        localDomRange = [(domainID, localxmin, localxmax)]
        globalDomRange = mpi.allreduce(localDomRange, mpi.SUM)
        globalDomRange.sort()
        for i in range(1, len(globalDomRange)):
            if globalDomRange[i][1] < globalDomRange[i-1][1]:
                self.fail("(proc,xmin) not in order")
            if globalDomRange[i][2] < globalDomRange[i-1][2]:
                self.fail("(proc,xmax) not in order")

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
