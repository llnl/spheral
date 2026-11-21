#ATS:test(SELF, np=8, label="TreeDistributedBoundary (3-D) unit tests")
import unittest
import random
import mpi

from Spheral import *
from DistributeNodes import *
from GenerateNodeDistribution3d import *
from SpheralTestUtilities import *

from generateGlobalIDs import *

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
numDomains = mpi.procs

#===============================================================================
# Main testing class for the 3-D tests.
#===============================================================================
class TestDistributedBoundary3d:

    #===========================================================================
    # Set up method called before test is run.
    #===========================================================================
    def genericSetup(self):

        # Generic parameters for 1-D tests.
        nx1, ny1, nz1 = 10, 10, 10
        nx2, ny2, nz2 = 20, 20, 20
        nx3, ny3, nz3 = 15, 15, 15

        xmin1, xmax1 = ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        xmin2, xmax2 = ((1.0, 0.0, 0.0), (2.0, 1.0, 1.0))
        xmin3, xmax3 = ((2.0, 0.0, 0.0), (3.0, 1.0, 1.0))

        generator1 = GenerateNodeDistribution3d(nx1, ny1, nz1, 1.0, 'lattice',
                                                xmin = xmin1,
                                                xmax = xmax1)
        generator2 = GenerateNodeDistribution3d(nx2, ny2, nz2, 1.0, 'lattice',
                                                xmin = xmin2,
                                                xmax = xmax2)
        generator3 = GenerateNodeDistribution3d(nx3, ny3, nz3, 1.0, 'lattice',
                                                xmin = xmin3,
                                                xmax = xmax3)

        # Interpolation kernel.
        self.WT = TableKernel3d(BSplineKernel3d(), 100)
        self.kernelExtent = self.WT.kernelExtent

        # Construct the NodeLists to be distributed
        self.eos = GammaLawGasMKS3d(2.0, 2.0)
        self.nodes1 = makeFluidNodeList3d("nodes 1", self.eos)
        self.nodes2 = makeFluidNodeList3d("nodes 2", self.eos)
        self.nodes3 = makeFluidNodeList3d("nodes 3", self.eos)

        # Distribute the nodes.
        distributeNodes3d((self.nodes1, generator1),
                          (self.nodes2, generator2),
                          (self.nodes3, generator3))

        # Assign global IDs.
        self.globalIDField1, self.globalIDField2, self.globalIDField3 = generateGlobalIDs((self.nodes1,
                                                                                           self.nodes2,
                                                                                           self.nodes3),
                                                                                          globalNodeIDs3d,
                                                                                          numGlobalNodes3d)

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase3d(self.kernelExtent)
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)
        self.dataBase.updateNeighbors()

        return

    #===========================================================================
    # The actual test itself!
    # Try setting up the distributed boundary and see if it correctly builds
    # the communicated ghost nodes.
    #===========================================================================
    def testIt(self):
        print("Testing TreeDistributedBoundary3d on domain %i of %i domains" % \
              (domainID, numDomains))

        # Set the ghost nodes for each domain distributed NodeList.
        self.domainbc.setAllGhostNodes(self.dataBase)
        self.domainbc.finalizeGhostBoundary()
        self.dataBase.updateNeighbors()

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
                for i in random.sample(list(range(n)), min(10, n)):

                    # Broadcast the position and H from the testing processor.
                    rilocal = Vector3d()
                    Hilocal = SymTensor3d()
                    if mpi.rank == testProc:
                        rilocal = nodes.positions()[i]
                        Hilocal = nodes.Hfield()[i]
                    ri = mpi.bcast(rilocal, testProc)
                    Hi = mpi.bcast(Hilocal, testProc)

                    # Get the global answer set for this node.
                    answer = mpi.reduce([self.globalIDField1[j] for j in findNeighborNodes(ri, Hi, self.kernelExtent, self.nodes1)] +
                                        [self.globalIDField2[j] for j in findNeighborNodes(ri, Hi, self.kernelExtent, self.nodes2)] +
                                        [self.globalIDField3[j] for j in findNeighborNodes(ri, Hi, self.kernelExtent, self.nodes3)],
                                        mpi.SUM, testProc)

                    # Have the testing processor build it's own version.
                    if mpi.rank == testProc:
                        masterNodes = vector_of_pair_ULL_ULL()
                        coarseNeighbors = vector_of_pair_ULL_ULL()
                        refineNeighbors = vector_of_pair_ULL_ULL()
                        neighbor = self.dataBase.neighbor()
                        ok, numMaster, numCoarse = neighbor.setMasterList(ri, Hi,
                                                                         masterNodes,
                                                                         coarseNeighbors,
                                                                         False)
                        assert ok
                        numRefine = neighbor.setRefineNeighborList(ri, Hi,
                                                                  numCoarse,
                                                                  coarseNeighbors,
                                                                  refineNeighbors)
                        refine = []
                        globalIDFields = [self.globalIDField1,
                                          self.globalIDField2,
                                          self.globalIDField3]
                        refine.extend(globalIDFields[nodeListID][nodeID]
                                      for nodeListID, nodeID in list(refineNeighbors)[:numRefine])

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
class TestTreeDistributedBoundary3d(unittest.TestCase,
                                          TestDistributedBoundary3d):

    # Set up method called before test is run.
    def setUp(self):
        self.domainbc = TreeDistributedBoundary3d.instance()
        self.genericSetup()

    def tearDown(self):
        del self.nodes1, self.nodes2, self.nodes3
        while gc.collect():
            pass
        return

#===============================================================================
# BoundingVolumeDistributedBoundary
#===============================================================================
class TestBoundingVolumeDistributedBoundary3d(unittest.TestCase,
                                              TestDistributedBoundary3d):

    # Set up method called before test is run.
    def setUp(self):
        self.domainbc = BoundingVolumeDistributedBoundary3d.instance()
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
