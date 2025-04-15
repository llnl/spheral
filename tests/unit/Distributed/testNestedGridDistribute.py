from math import *
import unittest

from Spheral import *
from testDistributeByPosition1d import *
from testDistributeByPosition2d import *
from testParmetisDistribute import *

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
nDomains = mpi.procs

#===============================================================================
# Main testing class for the 1-D tests.
#===============================================================================
class TestNestedGridRedistribute1d(TestDistributeByPosition1d):

    # The actual test itself!
    def runTest(self):
        print("Testing NestedGridRedistribute1d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in [self.dataBase.nodeLists[i] for i in range(self.dataBase.numNodeLists)]:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes, mpi.SUM))
        print("Total num nodes: ", nNodesGlobal, sum(nNodesGlobal))

        # Go ahead and redistribute those nodes!
        repartition = NestedGridRedistributeNodes1d(2.0)
        repartition.redistributeNodes(self.dataBase)

        # Make sure that the numbers of nodes are correct.
        assert self.dataBase.numNodeLists == len(nNodesGlobal)
        for nodeList, nGlobal in zip([self.dataBase.nodeLists[i] for i in range(self.dataBase.numNodeLists)],
                                     nNodesGlobal):
            n = mpi.allreduce(nodeList.numInternalNodes, mpi.SUM)
            if n != nGlobal:
                self.fail("Wrong number of nodes: %i != %i" % (n, nGlobal))

        # Do a consistency check of the new distribution.
        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase)
        if (not repartition.validDomainDecomposition(nodeDistribution,
                                                     self.dataBase)):
            self.fail("Invalid domain decomposition.")

        # Have each domain figure out it's min and max x.  This is
        # just a report, since it's not necessarily true that each domain
        # is exclusive in x range.
        localxmin = 1e10
        localxmax = -1e10
        for nodeList in [self.dataBase.nodeLists[i] for i in range(self.dataBase.numNodeLists)]:
            if nodeList.numInternalNodes > 0:
                localxmin = min(localxmin, min([r.x for r in [nodeList.positions().internalValues()[i] for i in range(nodeList.numInternalNodes())]]))
                localxmax = max(localxmax, max([r.x for r in [nodeList.positions().internalValues()[i] for i in range(nodeList.numInternalNodes())]]))

        sys.stderr.write("Process %i in x range (%f, %f)\n" % (domainID, localxmin, localxmax))

#===============================================================================
# Main testing class for the 2-D tests.
#===============================================================================
class TestNestedGridRedistribute2d(TestParmetisRedistribute2d):

    # The actual test itself!
    # Create a NestedGridRedistributeNodes object, have it redistribute the
    # nodes.
    def testDistribute(self):
        print("Testing NestedGridRedistributeNodes2d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

        # Go ahead and redistribute those nodes!
        repartition = NestedGridRedistributeNodes2d(2.0)
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
        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase)
        if (not repartition.validDomainDecomposition(nodeDistribution,
                                                     self.dataBase)):
            self.fail("Invalid domain decomposition.")

##         import SpheralGnuPlotUtilities
##         self.p = SpheralGnuPlotUtilities.plotNodePositions2d(self.dataBase,
##                                                              colorNodeLists = 0,
##                                                              colorDomains = 1)

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
