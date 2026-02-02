#ATS:test(SELF,        label="PolygonalMesh serial unit tests")
#ATS:test(SELF, np=4,  label="PolygonalMesh serial (2 proc) tests")
#ATS:test(SELF, np=9,  label="PolygonalMesh serial (4 proc) tests")
#ATS:test(SELF, np=16, label="PolygonalMesh serial (16 proc) tests")

from math import *
import unittest
import shutil, os
import random

from Spheral2d import *
from generateMesh import *
from siloMeshDump import siloMeshDump
from SpheralTestUtilities import fuzzyEqual

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
rank = mpi.rank
numDomains = mpi.procs
nxproc = int(sqrt(numDomains))
assert nxproc*nxproc == numDomains

#===============================================================================
# Create a global random number generator.
#===============================================================================
import random
random.seed(4599281940)

#===============================================================================
# Return a random string to help make test files unique.
#===============================================================================
def randomString():
    l = list(range(20))
    random.shuffle(l)
    result = ""
    for x in l:
        result += str(x)
    result = mpi.bcast(result, 0)
    return result

#===============================================================================
# Some boundary conditions.
#===============================================================================
x0, x1 = 0.0, 1.0
y0, y1 = 0.0, 1.0

nx, ny = 16*nxproc, 16*nxproc
n = nx*ny
nperdomain = n / numDomains
nxcell = KeyTraits.maxKey1d/4
nycell = nxcell
assert nx < nxcell
ncell = nxcell*nycell
dxcell = (x1 - x0)/nxcell
dycell = (y1 - y0)/nycell
dxhash = (x1 - x0)/(KeyTraits.maxKey1d - KeyTraits.two)

xmin = Vector(x0, y0)
xmax = Vector(x1, y1)

xbc0 = ReflectingBoundary(Plane(Vector(x0, y0), Vector(1.0, 0.0)))
ybc0 = ReflectingBoundary(Plane(Vector(x0, y0), Vector(0.0, 1.0)))
bclist = [xbc0, ybc0]
if numDomains > 1:
    bclist.append(BoundingVolumeDistributedBoundary.instance())
bclist = []

#===============================================================================
# Iterate the H's to something reasonable.
#===============================================================================
def iterateThoseHs(nodes):
    db = DataBase()
    db.appendNodeList(nodes)
    for bc in bclist:
        bc.setAllGhostNodes(db)
        nodes.neighbor().updateNodes()
    for bc in bclist:
        bc.finalizeGhostBoundary()
    nodes.neighbor().updateNodes()
    vecbound = vector_of_Boundary()
    for bc in bclist:
        vecbound.append(bc)
    WT = TableKernel(BSplineKernel(), 1000)
    smooth = SPHSmoothingScale()
    iterateIdealH(db, vecbound, WT, smooth,
                  tolerance = 1.0e-4)
    return

#===============================================================================
# A counter to help in creating unique NodeList names.
#===============================================================================
itest = 0

#===============================================================================
# Test class for tests to apply to all meshes.
#===============================================================================
class PolygonalMeshSiloGenericTests:

    #---------------------------------------------------------------------------
    # Test writing the mesh to a silo file.
    #---------------------------------------------------------------------------
    def testPolygonalMeshWriteSilo(self):
        mesh, void = generatePolygonalMesh([self.nodes],
                                           xmin = xmin,
                                           xmax = xmax,
                                           generateVoid = False,
                                           generateParallelConnectivity = True)
        siloMeshDump("silo_testPolygonalMeshWriteSilo_%s" % self.testext,
                     mesh,
                     label = "Test dumping a polygonal mesh",
                     nodeLists = [self.nodes])
        return

#===============================================================================
# Create a uniformly spaced nodes/mesh.
#===============================================================================
class UniformPolygonalMeshTests(unittest.TestCase, PolygonalMeshSiloGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        global itest
        self.testext = "Uniform%s_%idomain" % (randomString(), mpi.procs)
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes %i" % itest, eos,
                                       numInternal = nperdomain,
                                       nPerh = 2.01,
                                       hmin = 1e-5,
                                       hmax = 0.3)
        itest += 1
        self.pos = self.nodes.positions()
        self.H = self.nodes.Hfield()

        # Generate positions and split them up between domains appropriately.
        dxproc = (x1 - x0)/nxproc
        dyproc = (y1 - y0)/nxproc
        ixproc = rank % nxproc
        iyproc = rank / nxproc
        xminproc = Vector(x0 + ixproc*dxproc, y0 + iyproc*dyproc)
        xmaxproc = Vector(x0 + (ixproc + 1)*dxproc, y0 + (iyproc + 1)*dyproc)

        dxavg = (x1 - x0)/nx
        dyavg = (y1 - y0)/ny
        self.dxmin = dxavg
        xynodes_all = [Vector(x0 + (i % nx + 0.5)*dxavg, y0 + (i / nx + 0.5)*dyavg) for i in range(n)]
        xynodes = [v for v in xynodes_all if testPointInBox(v, xminproc, xmaxproc)]
        assert len(xynodes) == nperdomain
        assert mpi.allreduce(len(xynodes), mpi.SUM) == n

        # We now have the positions for each domain appropriately divided, so shuffle
        # the local positions.
        random.shuffle(xynodes)

        # Now we can set the node conditions.
        for i in range(nperdomain):
            self.pos[i] = xynodes[i]
            self.H[i] = SymTensor(1.0/(2.0*dxavg), 0.0,
                                  0.0, 1.0/(2.0*dyavg))
        self.nodes.neighbor().updateNodes()

        # Fix up the H's.
        #iterateThoseHs(self.nodes)
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        if mpi.rank == 0:
            os.remove("silo_testPolygonalMeshWriteSilo_%s.silo" % self.testext)
            shutil.rmtree("silo_testPolygonalMeshWriteSilo_%s" % self.testext, ignore_errors=True)
        return

#===============================================================================
# Create randomly spaced set of nodes in the unit square.
#===============================================================================
class RandomPolygonalMeshTests(unittest.TestCase, PolygonalMeshSiloGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        global itest
        self.testext = "Random%s_%idomain" % (randomString(), mpi.procs)
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes %i" % itest, eos,
                                       numInternal = nperdomain,
                                       nPerh = 2.01,
                                       hmin = 1.0e-5,
                                       hmax = 0.3)
        itest += 1
        self.pos = self.nodes.positions()
        self.H = self.nodes.Hfield()

        # Figure out the domain bounding volumes.
        dxproc = (x1 - x0)/nxproc
        dyproc = (y1 - y0)/nxproc
        ixproc = rank % nxproc
        iyproc = rank / nxproc
        xminproc = Vector(x0 + ixproc*dxproc, y0 + iyproc*dyproc)
        xmaxproc = Vector(x0 + (ixproc + 1)*dxproc, y0 + (iyproc + 1)*dyproc)

        # Randomly seed the generators.  We choose from random cells in order
        # to keep nodes from getting too close together.
        xynodes_all = []
        occupiedCells = set()
        for k in range(n):
            i = random.randint(0, ncell)
            while i in occupiedCells:
                i = random.randint(0, ncell)
            ix = i % nxcell
            iy = i / nxcell
            xynodes_all.append(Vector((ix + 0.5)*dxcell, (iy + 0.5)*dycell))
            occupiedCells.add(i)
        assert len(occupiedCells) == n
        xynodes_all = mpi.bcast(xynodes_all)
        xynodes = [v for v in xynodes_all if testPointInBox(v, xminproc, xmaxproc)]
        dxavg = (x1 - x0)/nx
        dyavg = (y1 - y0)/ny
        self.dxmin = dxavg
        assert mpi.allreduce(len(xynodes), mpi.SUM) == n

        # Now we can set the node conditions.
        self.nodes.numInternalNodes = len(xynodes)
        for i in range(len(xynodes)):
            self.pos[i] = xynodes[i]
            self.H[i] = SymTensor(1.0/(2.0*dxavg), 0.0,
                                  0.0, 1.0/(2.0*dyavg))
        self.nodes.neighbor().updateNodes()

        # Fix up the H's.
        iterateThoseHs(self.nodes)
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        if mpi.rank == 0:
            os.remove("silo_testPolygonalMeshWriteSilo_%s.silo" % self.testext)
            shutil.rmtree("silo_testPolygonalMeshWriteSilo_%s" % self.testext, ignore_errors=True)
        return

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
