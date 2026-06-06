from Spheral2d import *
import mpi
import os, random
from generateMesh import *
from siloMeshDump import *
from math import *
from testSharedElements import *

random.seed(4589281204)

x0, y0 = 0.0, 0.0
x1, y1 = 1.0, 1.0
nx = 10*mpi.procs

nxproc = int(sqrt(mpi.procs))
assert nxproc*nxproc == mpi.procs

dxproc = (x1 - x0)/nxproc
dyproc = (y1 - y0)/nxproc
ixproc = mpi.rank % nxproc
iyproc = mpi.rank / nxproc
xminproc = Vector(x0 + ixproc*dxproc, y0 + iyproc*dyproc)
xmaxproc = Vector(x0 + (ixproc + 1)*dxproc, y0 + (iyproc + 1)*dyproc)

fname = "generators_domain_%i_of_%i.txt" % (mpi.rank, mpi.procs)
if os.path.exists(fname):
    print("Reading existing gens.")
    f = open(fname, "r")
    xynodes = []
    for line in f:
        stuff = line.split()
        assert len(stuff) == 2
        xynodes.append(Vector(float(stuff[0]), float(stuff[1])))

else:
    # Randomly seed the generators.  We choose from random cells in order
    # to keep nodes from getting too close together.
    nxcell = KeyTraits.maxKey1d/4
    nycell = nxcell
    assert nx < nxcell
    ncell = nxcell*nycell
    dxcell = (x1 - x0)/nxcell
    dycell = (y1 - y0)/nycell
    xynodes_all = []
    occupiedCells = set()
    for k in range(nx*nx):
        i = random.randint(0, ncell)
        while i in occupiedCells:
            i = random.randint(0, ncell)
        ix = i % nxcell
        iy = i / nxcell
        xynodes_all.append(Vector((ix + 0.5)*dxcell, (iy + 0.5)*dycell))
        occupiedCells.add(i)
    assert len(occupiedCells) == nx*nx
    xynodes_all = mpi.bcast(xynodes_all)
    xynodes = [v for v in xynodes_all if testPointInBox(v, xminproc, xmaxproc)]

    # Write out the generator positions.
    f = open(fname, "w")
    for x in xynodes:
        f.write("%g %g\n" % (x.x, x.y))
    f.close()

# Build the set of generators.
gens = vector_of_Vector()
for x in xynodes:
    gens.append(x)

mesh = PolygonalMesh(gens,
                     xmin = Vector(x0, y0),
                     xmax = Vector(x1, y1))

#siloMeshDump("random_polygonal_mesh_%idomains" % mpi.procs, mesh)

# Test the mesh.
assert testSharedNodes(mesh)

# Now do the same thing through our NodeList interface.
eos = GammaLawGasMKS(5.0/3.0, 1.0)
W = TableKernel(BSplineKernel(), 1000)

def createNodes(gens, dx):
    nodes = makeFluidNodeList("some_nodes", eos,
                              numInternal = len(gens))
    pos = nodes.positions()
    H = nodes.Hfield()
    mass = nodes.mass()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    H0 = 1.0/(dx*nodes.nodesPerSmoothingScale) * SymTensor.one
    for i in range(len(gens)):
        xi = gens[i]
        pos[i] = xi
        H[i] = H0
        mass[i] = 1.0
        rho[i] = 1.0 + xi.magnitude2()
        vel[i] = xi

    db = DataBase()
    db.appendNodeList(nodes)
    nodes._neighbor.updateNodes()
    db.updateConnectivityMap()

    iterateIdealH(db, vector_of_Boundary(), W, SPHSmoothingScale(),
                  tolerance = 1.0e-4)
    db.updateConnectivityMap()

    return nodes, db

dx = 1.0/nx
nodes, db = createNodes(gens, dx)
mesh, void = generatePolygonalMesh([nodes], [],
                                   generateVoid = True,
                                   removeBoundaryZones = False)
siloMeshDump("random_polygonal_mesh_nodes_%idomains" % mpi.procs, mesh,
             nodeLists = [nodes, void])
