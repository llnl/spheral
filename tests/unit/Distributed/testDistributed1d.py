#ATS:test(SELF, np=8, label="TreeDistributedBoundary (1-D) unit tests")
import unittest
import mpi

from Spheral import *
from DistributeNodes import distributeNodesInRange1d
from SpheralTestUtilities import *

from generateGlobalIDs import *

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
domainID = mpi.rank
nDomains = mpi.procs

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
class TestDistributedBoundary1d:

    #===========================================================================
    # Set up method called before test is run.
    #===========================================================================
    def genericSetup(self):

        # Generic parameters for 1-D tests.
        nx1 = 100
        nx2 = 200
        nx3 = 50

        xRange1 = (0.0, 1.0)
        xRange2 = (0.0, 1.0)
        xRange3 = (0.5, 1.0)

        simulationVolume = (0.0, 1.0)

        neighborSearchType = GatherScatter
        numGridLevels = 20
        topGridCellSize = 0.25
        origin = Vector1d(0.0)

        # Interpolation kernel.
        self.WT = TableKernel1d(BSplineKernel1d(), 100)
        self.kernelExtent = self.WT.kernelExtent

        # Construct the NodeLists to be distributed
        self.eos = GammaLawGasMKS1d(2.0, 2.0)
        self.nodes1 = makeFluidNodeList1d("nodes1", self.eos, NeighborType = TreeNeighbor1d)
        self.nodes2 = makeFluidNodeList1d("nodes2", self.eos, NeighborType = TreeNeighbor1d)
        self.nodes3 = makeFluidNodeList1d("nodes3", self.eos, NeighborType = TreeNeighbor1d)

        # Distribute the nodes.
        distributeNodesInRange1d([(self.nodes1, nx1, 1.0, xRange1),
                                  (self.nodes2, nx2, 1.0, xRange2),
                                  (self.nodes3, nx3, 1.0, xRange3)])

        # Assign global IDs.
        self.globalIDField1, self.globalIDField2, self.globalIDField3 = generateGlobalIDs((self.nodes1,
                                                                                           self.nodes2,
                                                                                           self.nodes3),
                                                                                          globalNodeIDs1d,
                                                                                          numGlobalNodes1d)

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase1d()
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)
        print("Finished genericSetup")

        return

    #===========================================================================
    # The actual test itself!
    # Try setting up the distributed boundary and see if it correctly builds
    # the communicated ghost nodes.
    #===========================================================================
    def testIt(self):
        print("Testing TreeDistributedBoundary1d on domain %i of %i domains" % \
              (domainID, nDomains))

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
                for i in range(n):

                    # Broadcast the position and H from the testing processor.
                    rilocal = Vector1d()
                    Hilocal = SymTensor1d()
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
                            refine.sort()
                            answer.sort()
                            sys.stderr.write(" refine: %s\n" % str(refine))
                            sys.stderr.write(" answer: %s\n" % str(answer))
                            sys.stderr.write("  extra: %s\n" % str([x for x in refine if x not in answer]))
                            sys.stderr.write("missing: %s\n" % str([x for x in answer if x not in refine]))
                        else:
                            sys.stderr.write("PASSED for node %i\n" % i)
                        assert test

#===============================================================================
# TreeDistributedBoundary
#===============================================================================
class TestTreeDistributedBoundary1d(unittest.TestCase,
                                          TestDistributedBoundary1d):

    # Set up method called before test is run.
    def setUp(self):
        self.domainbc = TreeDistributedBoundary1d.instance()
        self.genericSetup()

    def tearDown(self):
        del self.nodes1, self.nodes2, self.nodes3
        while gc.collect():
            pass
        return

#===============================================================================
# BoundingVolumeDistributedBoundary
#===============================================================================
class TestBoundingVolumeDistributedBoundary1d(unittest.TestCase,
                                              TestDistributedBoundary1d):

    # Set up method called before test is run.
    def setUp(self):
        self.domainbc = BoundingVolumeDistributedBoundary1d.instance()
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
