#ATS:test(SELF, np=8, label="DistributeByPosition1d unit tests")
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
                     nDomains,     # total number of domains
                     nNodesGlobal, # global number of nodes in this nodelist
                     xRangeTotal): # total simulation volume

    assert thisDomain >= 0 and thisDomain < nDomains
    assert nDomains > 0

    import random
    g = random.Random()
    globalNodeIDs = []
    xNodePositions = []
    dxNodes = (xRangeTotal[1] - xRangeTotal[0])/nNodesGlobal
    for globalNodeID in range(0, nNodesGlobal):
        mpi.barrier()
        domain0 = g.randint(0, nDomains - 1)
        domain = mpi.bcast(domain0)
        if domain == thisDomain:
            globalNodeIDs.append(globalNodeID)
            xNodePositions.append(xRangeTotal[0] + (globalNodeID + 0.5)*dxNodes)

    assert len(globalNodeIDs) == len(xNodePositions)
    assert mpi.allreduce(len(globalNodeIDs), mpi.SUM) == nNodesGlobal
    return globalNodeIDs, xNodePositions

#===============================================================================
# Calculate one over the smoothing scale for the given number of nodes and
# volume.
#===============================================================================
def determineHField(nGlobal, nLocal, xRange,
                    nNodesPerh = 2.01):
    assert(xRange[1] - xRange[0] > 0.0)
    dx = (xRange[1] - xRange[0])/nGlobal
    H = [1.0/(nNodesPerh*dx)]*nLocal
    return H

#===============================================================================
# Main testing class for the 1-D tests.
#===============================================================================
class TestDistributeByPosition1d(unittest.TestCase):

    # Set up method called before test is run.
    def setUp(self):

        # Generic parameters for 1-D tests.
        nx1 = 100
        nx2 = 200
        nx3 = 50

        xRange1 = (0.0, 1.0)
        xRange2 = (1.0, 2.0)
        xRange3 = (2.0, 3.0)

        # Construct the NodeLists to be distributed
        self.eos = GammaLawGasMKS1d(2.0, 2.0)
        self.WT = TableKernel1d(BSplineKernel1d())
        self.nodes1 = makeFluidNodeList1d("nodes1", self.eos)
        self.nodes2 = makeFluidNodeList1d("nodes2", self.eos)
        self.nodes3 = makeFluidNodeList1d("nodes3", self.eos)
        for (nodes, nGlobal, globalRange) in ((self.nodes1, nx1, xRange1),
                                              (self.nodes2, nx2, xRange2),
                                              (self.nodes3, nx3, xRange3)):
            globalIDs, xNodes = randomDistribute(domainID,
                                                 nDomains,
                                                 nGlobal,
                                                 globalRange)
            n = len(globalIDs)
            nodes.numInternalNodes = n
            HNodes = determineHField(nGlobal, n, globalRange)
            for i in range(n):
                nodes.mass()[i] = 1.0
                nodes.positions()[i].x = xNodes[i]
                nodes.Hfield()[i].xx = HNodes[i]

            nodes.neighbor().updateNodes()

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase1d()
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)

        assert mpi.allreduce(self.nodes1.numInternalNodes, mpi.SUM) == nx1
        assert mpi.allreduce(self.nodes2.numInternalNodes, mpi.SUM) == nx2
        assert mpi.allreduce(self.nodes3.numInternalNodes, mpi.SUM) == nx3

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
        print("Testing DistributeByXPosition1d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

        # Go ahead and redistribute those nodes!
        repartition = DistributeByXPosition1d()
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
