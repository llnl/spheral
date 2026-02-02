#-------------------------------------------------------------------------------
# Simplified method to generate and distribute a single NodeList as a lattice
# on (n**nu) domains in nu dimensions.
#-------------------------------------------------------------------------------
from math import *

from NodeGeneratorBase import *

from Spheral import Vector2d, SymTensor2d
from Spheral import Vector3d, SymTensor3d

import mpi

#-------------------------------------------------------------------------------
# Compute the beginning node index given the domain index, size, and remainder.
#-------------------------------------------------------------------------------
def nodeindex(idomain, ndomain, rdomain):
    return ndomain*idomain + max(0, min(idomain, rdomain))

#-------------------------------------------------------------------------------
# 2-D (squares)
#-------------------------------------------------------------------------------
class GenerateSquareNodeDistribution(NodeGeneratorBase):

    def __init__(self,
                 nx,
                 ny,
                 rho,
                 xmin,
                 xmax,
                 nNodePerh = 2.01,
                 SPH = False):
        assert nx > 0
        assert ny > 0
        assert len(xmin) == 2
        assert len(xmax) == 2
        assert xmax[0] >= xmin[0]
        assert xmax[1] >= xmin[1]
        assert nNodePerh > 0.0

        # Remember the input.
        self.nx = nx
        self.ny = ny
        self.xmin = xmin
        self.xmax = xmax

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        # Compute the number of domains in each direction.
        lx = xmax[0] - xmin[0]
        ly = xmax[1] - xmin[1]
        nxdomains = int(sqrt(lx/ly * mpi.procs) + 0.1)
        nydomains = int(mpi.procs/nxdomains + 0.1)
        assert nxdomains * nydomains == mpi.procs
        
        # The number of nodes per domain.
        nxperdomain = nx // nxdomains
        nxremainder = nx % nxdomains
        nyperdomain = ny // nydomains
        nyremainder = ny % nydomains
        assert nxremainder < nxdomains
        assert nyremainder < nydomains

        # Compute the value for H.
        dx = lx/nx
        dy = ly/ny
        hx = nNodePerh * dx
        hy = nNodePerh * dy
        assert hx > 0.0 and hy > 0.0
        H0 = SymTensor2d(1.0/hx, 0.0,
                         0.0, 1.0/hy)
        if SPH:
            hxy = sqrt(hx*hy)
            H0 = SymTensor2d.one / hxy

        # The mass per node.
        m0 = lx*ly*rho / (nx*ny)
        assert m0 > 0.0

        # Compute our domain indicies.
        ixdomain = mpi.rank % nxdomains
        iydomain = mpi.rank // nxdomains
        ixmin = nodeindex(ixdomain, nxperdomain,  nxremainder)
        ixmax = nodeindex(ixdomain + 1, nxperdomain,  nxremainder)
        iymin = nodeindex(iydomain, nyperdomain,  nyremainder)
        iymax = nodeindex(iydomain + 1, nyperdomain,  nyremainder)
        assert ixmin < ixmax
        assert ixmin >= 0 and ixmax <= nx
        assert iymin < iymax
        assert iymin >= 0 and iymax <= ny

        # Now fill in the node values for this domain.
        self.x = []
        self.y = []
        self.m = []
        self.H = []
        for iy in range(iymin, iymax):
            for ix in range(ixmin, ixmax):
                self.x.append(xmin[0] + (ix + 0.5)*dx)
                self.y.append(xmin[1] + (iy + 0.5)*dy)
                self.m.append(m0)
                self.H.append(H0)
        assert mpi.allreduce(len(self.x), mpi.SUM) == nx*ny

        # Initialize the base class.
        NodeGeneratorBase.__init__(self, False)

        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector2d(self.x[i], self.y[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.x)
        return self.rho(self.localPosition(i))

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

#-------------------------------------------------------------------------------
# 3-D (cubes)
#-------------------------------------------------------------------------------
class GenerateCubicNodeDistribution(NodeGeneratorBase):

    def __init__(self,
                 nx,
                 ny,
                 nz,
                 rho,
                 xmin,
                 xmax,
                 nNodePerh = 2.01,
                 SPH = False):
        assert nx > 0
        assert ny > 0
        assert nz > 0
        assert len(xmin) == 3
        assert len(xmax) == 3
        assert xmax[0] >= xmin[0]
        assert xmax[1] >= xmin[1]
        assert xmax[2] >= xmin[2]
        assert nNodePerh > 0.0

        # Remember the input.
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.xmin = xmin
        self.xmax = xmax

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        # Compute the number of domains in each direction.
        lx = xmax[0] - xmin[0]
        ly = xmax[1] - xmin[1]
        lz = xmax[2] - xmin[2]
        nxdomains = int((lx*lx/(ly*lz) * mpi.procs)**(1.0/3.0) + 0.1)
        nydomains = int(ly/lx * nxdomains + 0.1)
        nzdomains = int(lz/lx * nxdomains + 0.1)
        assert nxdomains * nydomains * nzdomains == mpi.procs
        
        # The number of nodes per domain.
        nxperdomain = nx / nxdomains
        nxremainder = nx % nxdomains
        nyperdomain = ny / nydomains
        nyremainder = ny % nydomains
        nzperdomain = nz / nzdomains
        nzremainder = nz % nzdomains
        assert nxremainder < nxdomains
        assert nyremainder < nydomains
        assert nzremainder < nzdomains

        # Compute the value for H.
        dx = lx/nx
        dy = ly/ny
        dz = lz/nz
        hx = nNodePerh * dx
        hy = nNodePerh * dy
        hz = nNodePerh * dz
        assert hx > 0.0 and hy > 0.0 and hz > 0.0
        H0 = SymTensor3d(1.0/hx, 0.0, 0.0,
                         0.0, 1.0/hy, 0.0,
                         0.0, 0.0, 1.0/hz)
        if SPH:
            hxyz = (hx*hy*hz)**(1.0/3.0)
            H0 = SymTensor3d.one / hxyz

        # The mass per node.
        m0 = lx*ly*lz*rho / (nx*ny*nz)
        assert m0 > 0.0

        # Compute our domain indicies.
        ixdomain = mpi.rank % nxdomains
        iydomain = (mpi.rank / nxdomains) % nydomains
        izdomain = mpi.rank / (nxdomains*nydomains)
        ixmin = nodeindex(ixdomain, nxperdomain,  nxremainder)
        ixmax = nodeindex(ixdomain + 1, nxperdomain,  nxremainder)
        iymin = nodeindex(iydomain, nyperdomain,  nyremainder)
        iymax = nodeindex(iydomain + 1, nyperdomain,  nyremainder)
        izmin = nodeindex(izdomain, nzperdomain,  nzremainder)
        izmax = nodeindex(izdomain + 1, nzperdomain,  nzremainder)
        assert ixmin < ixmax
        assert ixmin >= 0 and ixmax <= nx
        assert iymin < iymax
        assert iymin >= 0 and iymax <= ny
        assert izmin < izmax
        assert izmin >= 0 and izmax <= nz


        # Now fill in the node values for this domain.
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        for iz in range(izmin, izmax):
            for iy in range(iymin, iymax):
                for ix in range(ixmin, ixmax):
                    self.x.append(xmin[0] + (ix + 0.5)*dx)
                    self.y.append(xmin[1] + (iy + 0.5)*dy)
                    self.z.append(xmin[2] + (iz + 0.5)*dz)
                    self.m.append(m0)
                    self.H.append(H0)
        assert mpi.allreduce(len(self.x), mpi.SUM) == nx*ny*nz

        # Initialize the base class.
        NodeGeneratorBase.__init__(self, False)

        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        assert len(self.x) == len(self.z)
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.x)
        return self.rho(self.localPosition(i))

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
