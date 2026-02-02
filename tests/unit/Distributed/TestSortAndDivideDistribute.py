from math import *
import unittest

from Spheral import *
from testDistributeByPosition1d import TestDistributeByPosition1d
from testDistributeByPosition2d import TestDistributeByPosition2d 
from testParmetisDistribute import testParmetisRedistribute2d, testParmetisRedistribute3d

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
nDomains = mpi.procs

#===============================================================================
# Main testing class for the 1-D tests.
#===============================================================================
class TestSortAndDivideRedistribute1d(TestDistributeByPosition1d):

    # The actual test itself!
    def testIt(self):
        print("Testing SortAndDivideRedistribute1d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in [self.dataBase.nodeLists[i] for i in range(self.dataBase.numNodeLists)]:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes, mpi.SUM))
        print("Total num nodes: ", nNodesGlobal, sum(nNodesGlobal))

        # Go ahead and redistribute those nodes!
        repartition = SortAndDivideRedistributeNodes1d(2.0)
        repartition.redistributeNodes(self.dataBase)

        # Make sure that the numbers of nodes are correct.
        assert self.dataBase.numNodeLists == len(nNodesGlobal)
        for nodeList, nGlobal in zip([self.dataBase.nodeLists[i] for i in range(self.dataBase.numNodeLists)],
                                     nNodesGlobal):
            n = mpi.allreduce(nodeList.numInternalNodes, mpi.SUM)
            if n != nGlobal:
                self.fail("Wrong number of nodes: %i != %i" % (n, nGlobal))

##        # Do a consistency check of the new distribution.
##        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase)
##        if (not repartition.validDomainDecomposition(nodeDistribution,
##                                                     self.dataBase)):
##            self.fail("Invalid domain decomposition.")

        # Have each domain figure out it's min and max x.  This is
        # just a report, since it's not necessarily true that each domain
        # is exclusive in x range.
        localxmin = 1e10
        localxmax = -1e10
        for nodeList in [self.dataBase.nodeLists[i] for i in range(self.dataBase.numNodeLists)]:
            if nodeList.numInternalNodes > 0:
                localxmin = min(localxmin, min([r.x for r in nodeList.positions().internalValues()] + [1e60]))
                localxmax = max(localxmax, max([r.x for r in nodeList.positions().internalValues()] + [-1e60]))

        import sys
        sys.stderr.write("Process %i in x range (%f, %f)\n" % (domainID, localxmin, localxmax))

        # Build a diagnostic plot to help show how the domains are distributed.
        domain = ScalarFieldList1d()
        domain.copyFields()
        for nodes in self.dataBase.nodeLists:
            f = ScalarField1d(nodes.name(), nodes, mpi.rank)
            domain.appendField(f)
        self.p = plotFieldList(domain,
                               plotStyle = "points",
                               colorNodeLists = False,
                               colorDomains = True)

#===============================================================================
# Main testing class for the 2-D tests.
#===============================================================================
class TestSortAndDivideRedistribute2d(TestParmetisRedistribute2d):

    # The actual test itself!
    # Create a SortAndDivideRedistributeNodes object, have it redistribute the
    # nodes.
    def testIt(self):
        print("Testing SortAndDivideRedistributeNodes2d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

        # Go ahead and redistribute those nodes!
        repartition = SortAndDivideRedistributeNodes2d(2.0)
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

##        # Do a consistency check of the new distribution.
##        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase)
##        if (not repartition.validDomainDecomposition(nodeDistribution,
##                                                     self.dataBase)):
##            self.fail("Invalid domain decomposition.")

        self.p = plotNodePositions2d(self.dataBase,
                                     colorNodeLists = 0,
                                     colorDomains = 1)

#===============================================================================
# Main testing class for the 3-D tests.
#===============================================================================
class TestSortAndDivideRedistribute3d(TestParmetisRedistribute3d):

    # The actual test itself!
    # Create a SortAndDivideRedistributeNodes object, have it redistribute the
    # nodes.
    def testIt(self):
        print("Testing SortAndDivideRedistributeNodes3d on domain %i of %i domains" % \
              (domainID, nDomains))

        # Record how many nodes we're starting with.
        nNodesGlobal = []
        for nodeList in self.dataBase.nodeLists:
            nNodesGlobal.append(mpi.allreduce(nodeList.numInternalNodes,
                                              mpi.SUM))

        # Go ahead and redistribute those nodes!
        repartition = SortAndDivideRedistributeNodes3d(2.0)
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

##        # Do a consistency check of the new distribution.
##        nodeDistribution = repartition.currentDomainDecomposition(self.dataBase)
##        if (not repartition.validDomainDecomposition(nodeDistribution,
##                                                     self.dataBase)):
##            self.fail("Invalid domain decomposition.")
