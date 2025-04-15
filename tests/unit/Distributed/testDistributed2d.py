#ATS:test(SELF, np=8, label="TreeDistributedBoundary (2-D) unit tests")
import unittest
import random
import mpi

from Spheral import *
from DistributeNodes import *
from GenerateNodeDistribution2d import *
from SpheralTestUtilities import *

from generateGlobalIDs import *

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
numDomains = mpi.procs

#===============================================================================
# Main testing class for the 2-D tests.
#===============================================================================
class TestDistributedBoundary2d:

    #===========================================================================
    # Set up method called before test is run.
    #===========================================================================
    def genericSetup(self):

        # Generic parameters for 1-D tests.
        nx1, ny1 = 50, 25
        nx2, ny2 = 25, 25
        nx3, ny3 = 40, 40

        xmin1, xmax1 = ((0.0, 0.0), (1.0, 0.5))
        xmin2, xmax2 = ((0.0, 0.5), (0.5, 1.0))
        xmin3, xmax3 = ((0.5, 0.5), (1.0, 1.0))

        generator1 = GenerateNodeDistribution2d(nx1, ny1, 1.0, 'lattice',
                                                xmin = xmin1,
                                                xmax = xmax1)
        generator2 = GenerateNodeDistribution2d(nx2, ny2, 1.0, 'lattice',
                                                xmin = xmin2,
                                                xmax = xmax2)
        generator3 = GenerateNodeDistribution2d(nx3, ny3, 1.0, 'lattice',
                                                xmin = xmin3,
                                                xmax = xmax3)

        # Interpolation kernel.
        self.WT = TableKernel2d(BSplineKernel2d(), 100)
        self.kernelExtent = self.WT.kernelExtent

        # Construct the NodeLists to be distributed
        self.eos = GammaLawGasMKS2d(2.0, 2.0)
        self.nodes1 = makeFluidNodeList2d("nodes1", self.eos, NeighborType = TreeNeighbor2d)
        self.nodes2 = makeFluidNodeList2d("nodes2", self.eos, NeighborType = TreeNeighbor2d)
        self.nodes3 = makeFluidNodeList2d("nodes3", self.eos, NeighborType = TreeNeighbor2d)

        # Distribute the nodes.
        distributeNodes2d((self.nodes1, generator1),
                          (self.nodes2, generator2),
                          (self.nodes3, generator3))

        # Assign global IDs.
        self.globalIDField1, self.globalIDField2, self.globalIDField3 = generateGlobalIDs((self.nodes1,
                                                                                           self.nodes2,
                                                                                           self.nodes3),
                                                                                          globalNodeIDs2d,
                                                                                          numGlobalNodes2d)

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase2d()
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)

        return

    #===========================================================================
    # The actual test itself!
    # Try setting up the distributed boundary and see if it correctly builds
    # the communicated ghost nodes.
    #===========================================================================
    def testIt(self):
        print("Testing TreeDistributedBoundary2d on domain %i of %i domains" % \
              (domainID, numDomains))

        # Set the ghost nodes for each domain distributed NodeList.
        self.domainbc.setAllGhostNodes(self.dataBase)
        self.domainbc.finalizeGhostBoundary()
        for nodes in self.dataBase.nodeLists:
            nodes.neighbor().updateNodes()

        # Exchange the global node ID fields.
        self.domainbc.applyGhostBoundary(self.globalIDField1)
        self.domainbc.applyGhostBoundary(self.globalIDField2)
        self.domainbc.applyGhostBoundary(self.globalIDField3)
        self.domainbc.finalizeGhostBoundary()

        # Iterate over each domain.
        for testProc in range(mpi.procs):

            # Test each NodeList.
            for (nodes, globalIDField) in ((self.nodes1, self.globalIDField1),
                                           (self.nodes2, self.globalIDField2),
                                           (self.nodes3, self.globalIDField3)):

                # Tell everyone how many nodes we'll be testing, and iterate
                # over them
                n = mpi.bcast(nodes.numInternalNodes, testProc)
                for i in random.sample(list(range(n)), min(100, n)):

                    # Broadcast the position and H from the testing processor.
                    rilocal = Vector2d()
                    Hilocal = SymTensor2d()
                    if mpi.rank == testProc:
                        rilocal = nodes.positions()[i]
                        Hilocal = nodes.Hfield()[i]
                    ri = mpi.bcast(rilocal, testProc)
                    Hi = mpi.bcast(Hilocal, testProc)

                    # Get the global answer set for this node.
                    answer = mpi.allreduce([self.globalIDField1[j] for j in findNeighborNodes(ri, Hi, self.kernelExtent, self.nodes1)] +
                                           [self.globalIDField2[j] for j in findNeighborNodes(ri, Hi, self.kernelExtent, self.nodes2)] +
                                           [self.globalIDField3[j] for j in findNeighborNodes(ri, Hi, self.kernelExtent, self.nodes3)],
                                           mpi.SUM)

                    # Have the testing processor build it's own version.
                    if mpi.rank == testProc:
                        masterLists = vector_of_vector_of_int()
                        coarseNeighbors = vector_of_vector_of_int()
                        refineNeighbors = vector_of_vector_of_int()
                        self.dataBase.setMasterNodeLists(ri, Hi, masterLists, coarseNeighbors, False)
                        self.dataBase.setRefineNodeLists(ri, Hi, coarseNeighbors, refineNeighbors)
                        assert len(refineNeighbors) == 3
                        refine = []
                        for k, globalIDs in enumerate([self.globalIDField1,
                                                       self.globalIDField2,
                                                       self.globalIDField3]):
                            refine.extend([globalIDs[j] for j in refineNeighbors[k]])

                        # Check the answer.
                        test = checkNeighbors(refine, answer)
                        if not test:
                            sys.stderr.write("FAILED for node %i\n" % i)
                        else:
                            sys.stderr.write("PASSED for node %i of %i\n" % (i, n))
                        assert test

#===============================================================================
# TreeDistributedBoundary
#===============================================================================
class TestTreeDistributedBoundary2d(unittest.TestCase,
                                          TestDistributedBoundary2d):

    # Set up method called before test is run.
    def setUp(self):
        self.domainbc = TreeDistributedBoundary2d.instance()
        self.genericSetup()

    def tearDown(self):
        del self.nodes1, self.nodes2, self.nodes3
        while gc.collect():
            pass
        return

#===============================================================================
# BoundingVolumeDistributedBoundary
#===============================================================================
class TestBoundingVolumeDistributedBoundary2d(unittest.TestCase,
                                              TestDistributedBoundary2d):

    # Set up method called before test is run.
    def setUp(self):
        self.domainbc = BoundingVolumeDistributedBoundary2d.instance()
        self.genericSetup()

    def tearDown(self):
        del self.nodes1, self.nodes2, self.nodes3
        while gc.collect():
            pass
        return

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
