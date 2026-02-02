#ATS:test(SELF,        label="LineMesh serial unit tests")
#ATS:test(SELF, np=2,  label="LineMesh parallel (2 proc) tests")
#ATS:test(SELF, np=4,  label="LineMesh parallel (4 proc) tests")
#ATS:test(SELF, np=10, label="LineMesh parallel (10 proc) tests")

from math import *
import unittest

from Spheral1d import *
from generateMesh import *
from SpheralTestUtilities import fuzzyEqual, testParallelConsistency

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
rank = mpi.rank
numDomains = mpi.procs

#===============================================================================
# Compute the min & max scales.
#===============================================================================
def meshScales(xnodes, xmin, xmax):
    nx = len(xnodes)
    xsort = list(xnodes)
    xsort.sort()
    dx = ([(0.5*(xsort[i + 1] + xsort[i]) - 0.5*(xsort[i] + xsort[i - 1])) for i in range(1, nx-1)] + 
          [0.5*(xsort[0] + xsort[1]) - xmin] +
          [xmax - 0.5*(xsort[-2] + xsort[-1])])
    dxmin = mpi.bcast(0.5*min(dx))
    dxmax = mpi.bcast(0.5*max(dx))
    return dxmin, dxmax

#===============================================================================
# Create a global random number generator.
#===============================================================================
import random
random.seed(589290234)

#===============================================================================
# Some boundary conditions.
#===============================================================================
x0, x1 = 0.0, 1.0
xmin = Vector(x0)
xmax = Vector(x1)
xbc0 = ReflectingBoundary(Plane(Vector(x0), Vector(1.0)))
bclist = [xbc0]
if numDomains > 1:
    bclist.append(BoundingVolumeDistributedBoundary.instance())

nx = 100
assert nx % numDomains == 0

#===============================================================================
# Test class for tests to apply to all meshes.
#===============================================================================
class LineMeshGenericTests:

    #---------------------------------------------------------------------------
    # Test numbers of elements.
    #---------------------------------------------------------------------------
    def testLineMeshNums0(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax,
                                      generateVoid = False,
                                      generateParallelConnectivity = False)
        assert mesh.numNodes == self.nodes.numInternalNodes + 1
        assert mesh.numEdges == self.nodes.numInternalNodes + 1
        assert mesh.numFaces == self.nodes.numInternalNodes + 1
        assert mesh.numZones == self.nodes.numInternalNodes

    #---------------------------------------------------------------------------
    # Test element IDs.
    #---------------------------------------------------------------------------
    def testLineMeshElementIDs(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax)
        sys.stderr.write("%i %i %i\n" % (self.nodes.numInternalNodes, self.nodes.numGhostNodes, mesh.numZones))
        for i in range(self.nodes.numInternalNodes):
            node = mesh.node(i)
            assert node.ID == i
        for i in range(mesh.numEdges):
            edge = mesh.edge(i)
            assert edge.ID == i
        for i in range(mesh.numFaces):
            face = mesh.face(i)
            assert face.ID == i
        for i in range(mesh.numZones):
            zone = mesh.zone(i)
            assert zone.ID == i
        return

    #---------------------------------------------------------------------------
    # Test that the zones in the mesh correspond to the correct seed nodes.
    #---------------------------------------------------------------------------
    def testLineMeshZoneOrder(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax)
        assert mesh.numZones >= self.nodes.numInternalNodes
        pos = self.nodes.positions()
        for i in range(self.nodes.numInternalNodes):
            zone = mesh.zone(i)
            zonehull = zone.convexHull
            self.assertTrue(zonehull.contains(pos[i]),
                            "Zone does not contain generator %s %s %s %s" %
                            (pos[i], zone.position, 
                             mesh.node(zone.nodeIDs[0]).position,
                             mesh.node(zone.nodeIDs[1]).position))
        return

    #---------------------------------------------------------------------------
    # Test that the nodes are sorted in x.
    #---------------------------------------------------------------------------
    def testLineMeshNodeOrder(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax)
        for i in range(mesh.numNodes - 1):
            assert mesh.node(i).position.x < mesh.node(i + 1).position.x
        return

    #---------------------------------------------------------------------------
    # Test the minimum scale.
    #---------------------------------------------------------------------------
    def testLineMeshMinimumScale(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax)
        self.assertTrue(mesh.minimumScale <= self.dxmin, 
                        "Scales don't match:  %g %g" % (mesh.minimumScale, self.dxmin))
        return

    #---------------------------------------------------------------------------
    # Test the parallel domain info.
    #---------------------------------------------------------------------------
    def testLineMeshParallel(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax,
                                      generateParallelConnectivity = True,
                                      generateVoid = False,
                                      removeBoundaryZones = False)

        msg = testParallelConsistency(mesh, xmin, xmax)
        self.assertTrue(msg == "ok", msg)

        # neighborDomains = [int(x) for x in mesh.neighborDomains]
        # sharedNodes = []
        # sharedFaces = []
        # for ll in mesh.sharedNodes:
        #     sharedNodes.append([int(x) for x in ll])
        # for ll in mesh.sharedFaces:
        #     sharedFaces.append([int(x) for x in ll])
        # assert len(neighborDomains) == len(sharedNodes)
        # assert len(neighborDomains) == len(sharedFaces)

        # # for irank in xrange(numDomains):
        # #     if rank == irank:
        # #         sys.stderr.write("node positions : %s\n" % [mesh.node(i).position.x for i in xrange(mesh.numNodes)])
        # #         sys.stderr.write("zone volumes   : %s\n" % [mesh.zone(i).volume() for i in xrange(mesh.numZones)])
        # #         sys.stderr.write("mesh.neighborDomains : %s\n" % neighborDomains)
        # #         for i in xrange(len(neighborDomains)):
        # #             sys.stderr.write("    mesh.sharedNodes[%i] : %s\n" % (neighborDomains[i], sharedNodes[i]))
        # #         for i in xrange(len(neighborDomains)):
        # #             sys.stderr.write("    mesh.sharedFaces[%i] : %s\n" % (neighborDomains[i], sharedFaces[i]))
        # #     mpi.barrier()

        # # Check the correct domains are talking to each other.
        # neighborDomainsAnswer = [i for i in xrange(max(0, rank - 1), min(numDomains, rank + 2)) if i != rank]
        # ok = mpi.allreduce((neighborDomains == neighborDomainsAnswer), mpi.MIN)
        # self.assertTrue(ok, "Strange neighbor domains for %i : %s ?= %s" % (rank, str(neighborDomains), str(neighborDomainsAnswer)))
      
        # # Check that the communicated mesh nodes are consistent.
        # boxInv = Vector(1.0/(xmax.x - xmin.x))
        # for sendProc in xrange(numDomains):
        #     numChecks = mpi.bcast(len(neighborDomains), root=sendProc)
        #     assert mpi.allreduce(numChecks, mpi.MIN) == mpi.allreduce(numChecks, mpi.MAX)
        #     for k in xrange(numChecks):
        #         if rank == sendProc:
        #             ksafe = k
        #         else:
        #             ksafe = 0
        #         recvProc = mpi.bcast(neighborDomains[ksafe], root=sendProc)
        #         recvHashes = mpi.bcast([hashPosition(mesh.node(i).position, xmin, xmax, boxInv) for i in sharedNodes[ksafe]], root=sendProc)
        #         ok = True
        #         msg = ""
        #         if rank == recvProc:
        #             assert sendProc in neighborDomains
        #             kk = neighborDomains.index(sendProc)
        #             assert kk < len(sharedNodes)
        #             ok = ([hashPosition(mesh.node(i).position, xmin, xmax, boxInv) for i in sharedNodes[kk]] == recvHashes)
        #             msg = "Shared node indicies don't match %i %i : %s != %s" % (rank, sendProc, 
        #                                                                          str([hashPosition(mesh.node(i).position, xmin, xmax, boxInv) for i in sharedNodes[kk]]),
        #                                                                          recvHashes)
        #         self.assertTrue(mpi.allreduce(ok, mpi.MIN), msg)

        # return

    #---------------------------------------------------------------------------
    # Test each mesh point hashes uniquely.
    #---------------------------------------------------------------------------
    def testLineMeshHash(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax)
        pos = [mesh.zone(i).position for i in range(mesh.numZones)] + [mesh.node(i).position for i in range(mesh.numNodes)]
        xpos = [x.x for x in pos]
        xpos.sort()
        boxInv = Vector(1.0/(xmax.x - xmin.x))
        hashes = [hashPosition(x, xmin, xmax, boxInv) for x in pos]
        hashes.sort()
        for x in hashes:
            if hashes.count(x) != 1:
                k = hashes.index(x)
                self.fail("%i: Non unique hash:  %s %i %i \n%g %g\n%s\n%s" % (mpi.rank, x, hashes.count(x), k, xpos[k], xpos[k + 1] - xpos[k], str(hashes), str(xpos)))
        return

    #---------------------------------------------------------------------------
    # Test the zones of the nodes.
    #---------------------------------------------------------------------------
    def testLineNodeZones(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax)
        zoneList = [(mesh.zone(i).position.x, i) for i in range(mesh.numZones)]
        zoneList.sort()
        for inode in range(1, mesh.numNodes - 1):
            zoneIDs = mesh.node(inode).zoneIDs
            assert len(zoneIDs) <= 2
            assert ((inode == 0             and len(zoneIDs) == 1 and zoneIDs[0] == zoneList[0][1]) or
                    (inode == mesh.numZones and len(zoneIDs) == 1 and zoneIDs[0] == zoneList[-1][1]) or
                    (len(zoneIDs) == 2 and
                     zoneIDs[0] == zoneList[inode - 1][1] and
                     zoneIDs[1] == zoneList[inode][1]))
        return

    #---------------------------------------------------------------------------
    # Test the bounding surface.
    #---------------------------------------------------------------------------
    def testBoundingSurface(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax,
                                      generateParallelConnectivity = True)
        bs = mesh.boundingSurface()
        assert fuzzyEqual(bs.xmin.x, x0, 1.0e-10)
        assert fuzzyEqual(bs.xmax.x, x1, 1.0e-10)

        # Check that all the generators are contained.
        pos = self.nodes.positions()
        for i in range(self.nodes.numInternalNodes):
            self.assertTrue(bs.contains(pos[i]),
                            "Failed containment for generator %i @ %s" % (i, pos[i]))

        # Check that all mesh nodes are contained.
        for i in range(mesh.numNodes):
            self.assertTrue(bs.contains(mesh.node(i).position),
                            "Failed containment for mesh node %i @ %s" % (i, mesh.node(i).position))

        return

#===============================================================================
# Create a uniformly spaced nodes/mesh.
#===============================================================================
class UniformLineMeshTests(unittest.TestCase, LineMeshGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        nxperdomain = nx // numDomains
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes", eos,
                                       numInternal = nxperdomain,
                                       nPerh = 2.01)
        pos = self.nodes.positions()
        H = self.nodes.Hfield()

        # Generate initial positions, and split them up between domains appropriately.
        dxavg = (x1 - x0)/nx
        xnodes = [x0 + (i + 0.5)*dxavg for i in range(nx)] # [random.uniform(x0, x1) for i in xrange(nx)]
        xnodes.sort()
        self.dxmin, self.dxmax = meshScales(xnodes, x0, x1)
        for proc in range(numDomains):
            xnodes[proc*nxperdomain : (proc + 1)*nxperdomain] = mpi.bcast(xnodes[proc*nxperdomain : (proc + 1)*nxperdomain])
        xnodes = xnodes[rank*nxperdomain : (rank + 1)*nxperdomain]
        assert len(xnodes) == nxperdomain
        assert mpi.allreduce(len(xnodes), mpi.SUM) == nx

        # We now have the positions for each domain appropriately divided, so shuffle
        # the local positions.
        random.shuffle(xnodes)

        # Now we can set the node conditions.
        for i in range(nxperdomain):
            pos[i] = Vector(xnodes[i])
            H[i] = SymTensor(1.0/(2.0*self.dxmax))
        self.nodes.neighbor().updateNodes()

        # Iterate the H tensors to somthing reasonable.
        db = DataBase()
        db.appendNodeList(self.nodes)
        for bc in bclist:
            bc.setAllGhostNodes(db)
        for bc in bclist:
            bc.finalizeGhostBoundary()
        db = DataBase()
        db.appendNodeList(self.nodes)
        WT = TableKernel(BSplineKernel(), 1000)
        smooth = SPHSmoothingScale(IdealH, WT)
        iterateIdealH(db,
                      vector_of_Physics([smooth]),
                      vector_of_Boundary(bclist))
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        return

    #---------------------------------------------------------------------------
    # Test the generation of the void node at the extreme end.
    #---------------------------------------------------------------------------
    def testUniformLineMeshVoidNode(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = 2.0*xmax)
        voidpos = void.positions()
        pos = self.nodes.positions()
        self.assertTrue(mpi.allreduce(void.numNodes, mpi.SUM) == 1, 
                        "Bad number of void nodes:  %i %s" % (mpi.allreduce(void.numNodes, mpi.SUM), str([x.x for x in void.positions().allValues()])))
        assert mpi.allreduce(mesh.numZones, mpi.SUM) == mpi.allreduce(self.nodes.numInternalNodes, mpi.SUM) + 1
        maxpos = mpi.allreduce(max([pos[i].x for i in range(self.nodes.numInternalNodes)]), mpi.MAX)
        if void.numNodes == 1:
            self.assertTrue(voidpos[0].x > maxpos, "%f %f" % (voidpos[0].x, maxpos))
            voidzone = mesh.zone(self.nodes.numInternalNodes)
            voidhull = voidzone.convexHull
            assert voidhull.contains(voidpos[0])
        return

#===============================================================================
# Create uniformly spaced nodes/mesh with a gap.
#===============================================================================
class UniformGapLineMeshTests(unittest.TestCase, LineMeshGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        nxperdomain = nx // numDomains
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes", eos,
                                       numInternal = nxperdomain,
                                       nPerh = 2.01)
        pos = self.nodes.positions()
        H = self.nodes.Hfield()

        # Generate initial positions, and split them up between domains appropriately.
        gap = 0.9
        dxavg = (x1 - x0 - gap)/nx
        xnodes = [x0 + (i + 0.5)*dxavg for i in range(nx//2)] + [0.5*(x0 + x1 + gap) + (i + 0.5)*dxavg for i in range(nx//2)]
        xnodes.sort()
        self.dxmin, self.dxmax = meshScales(xnodes, x0, x1)
        for proc in range(numDomains):
            xnodes[proc*nxperdomain : (proc + 1)*nxperdomain] = mpi.bcast(xnodes[proc*nxperdomain : (proc + 1)*nxperdomain])
        xnodes = xnodes[rank*nxperdomain : (rank + 1)*nxperdomain]
        assert len(xnodes) == nxperdomain
        assert mpi.allreduce(len(xnodes), mpi.SUM) == nx

        # We now have the positions for each domain appropriately divided, so shuffle
        # the local positions.
        #random.shuffle(xnodes)

        # Now we can set the node conditions.
        for i in range(nxperdomain):
            pos[i] = Vector(xnodes[i])
            H[i] = SymTensor(1.0/(2.0*dxavg))
        self.nodes.neighbor().updateNodes()

        # Iterate the H tensors to somthing reasonable.
        db = DataBase()
        db.appendNodeList(self.nodes)
        for bc in bclist:
            bc.setAllGhostNodes(db)
        for bc in bclist:
            bc.finalizeGhostBoundary()
        db = DataBase()
        db.appendNodeList(self.nodes)
        WT = TableKernel(BSplineKernel(), 1000)
        smooth = SPHSmoothingScale(IdealH, WT)
        iterateIdealH(db,
                      vector_of_Physics([smooth]),
                      vector_of_Boundary(bclist))
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        return

    #---------------------------------------------------------------------------
    # Test the generation of void nodes in the gap.
    #---------------------------------------------------------------------------
    def testUniformLineMeshVoidNodes(self):
        mesh, void = generateLineMesh([self.nodes],
                                      xmin = xmin,
                                      xmax = xmax,
                                      generateVoid = True)
        voidpos = void.positions()
        self.assertTrue(mpi.allreduce(void.numNodes, mpi.SUM) == 2, 
                        "Bad number of void nodes:  %i %s" % (mpi.allreduce(void.numNodes, mpi.SUM), str([x.x for x in void.positions().allValues()])))
        self.assertTrue(mpi.allreduce(mesh.numZones, mpi.SUM) == mpi.allreduce(self.nodes.numInternalNodes, mpi.SUM) + 2,
                        "Bad number of mesh zones:  %i %i %i" % (mpi.allreduce(mesh.numZones, mpi.SUM),
                                                                 mpi.allreduce(self.nodes.numInternalNodes, mpi.SUM),
                                                                 mpi.allreduce(void.numInternalNodes, mpi.SUM)))
        return

#===============================================================================
# Create a NodeList and mesh based on randomly distributed nodes.
#===============================================================================
class RandomLineMeshTests(unittest.TestCase, LineMeshGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        nxperdomain = nx // numDomains
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes", eos,
                                       numInternal = nxperdomain,
                                       nPerh = 2.01)
        pos = self.nodes.positions()
        H = self.nodes.Hfield()

        # Generate initial positions, and split them up between domains appropriately.
        dxavg = (x1 - x0)/nx
        xnodes = [random.uniform(x0, x1) for i in range(nx)]
        xnodes.sort()
        self.dxmin, self.dxmax = meshScales(xnodes, x0, x1)
        for proc in range(numDomains):
            xnodes[proc*nxperdomain : (proc + 1)*nxperdomain] = mpi.bcast(xnodes[proc*nxperdomain : (proc + 1)*nxperdomain])
        xnodes = xnodes[rank*nxperdomain : (rank + 1)*nxperdomain]
        assert len(xnodes) == nxperdomain
        assert mpi.allreduce(len(xnodes), mpi.SUM) == nx

        # We now have the positions for each domain appropriately divided, so shuffle
        # the local positions.
        random.shuffle(xnodes)

        # Now we can set the node conditions.
        for i in range(nxperdomain):
            pos[i] = Vector(xnodes[i])
            H[i] = SymTensor(1.0/(2.0*self.dxmax))
        self.nodes.neighbor().updateNodes()

        # Iterate the H tensors to somthing reasonable.
        db = DataBase()
        db.appendNodeList(self.nodes)
        for bc in bclist:
            bc.setAllGhostNodes(db)
        for bc in bclist:
            bc.finalizeGhostBoundary()
        db = DataBase()
        db.appendNodeList(self.nodes)
        WT = TableKernel(BSplineKernel(), 1000)
        smooth = SPHSmoothingScale(IdealH, WT)
        iterateIdealH(db,
                      vector_of_Physics([smooth]),
                      vector_of_Boundary(bclist))
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        return

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
