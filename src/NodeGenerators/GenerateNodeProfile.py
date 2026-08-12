from math import *
import numpy as np

from NodeGeneratorBase import *

from Spheral import (Vector1d, Tensor1d, SymTensor1d,
                     Vector2d, Tensor2d, SymTensor2d, rotationMatrix2d, testPointInBox2d,
                     Vector3d, Tensor3d, SymTensor3d, rotationMatrix3d, testPointInBox3d)

#-------------------------------------------------------------------------------
# Class to generate 1-D node positions for a fixed node mass to fit the given
# density profile in a range (xmin, xmax).
#-------------------------------------------------------------------------------
class GenerateNodeProfile1d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 nx,                 # number of points to generate
                 rho,                # density profile
                 xmin,
                 xmax,
                 nNodePerh = 2.01,
                 numbins = 10000):

        assert nx > 0
        assert xmin < xmax
        assert nNodePerh > 0.0

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rhofunc = ConstantRho(rho)

            # In the constant rho case, no need to kill ourselves figuring out complicated fits...
            dx = (xmax - xmin)/nx
            mi = dx*rho
            self.x = [xmin + (i+0.5)*dx for i in range(nx)]
            self.H = [SymTensor1d(1.0/(nNodePerh*dx)) for i in range(nx)]
            self.m = [mi]*nx
            self.rho = [rho]*nx

        else:
            self.rhofunc = rho

            # Build the evenly sampled cumulative mass as a function of position.
            ok = False
            while not ok:
                dx = (xmax - xmin)/numbins
                mcum = np.cumsum(np.array([0.0] + [0.5*dx*(self.rhofunc(xmin + i*dx) + self.rhofunc(xmin + (i + 1)*dx)) for i in range(numbins)]))

                # Find the target mass per node.
                mi = mcum[-1]/nx

                # Do we need to have a finer binning?
                if mcum[-1]/mi > 0.5*numbins:
                    numbins = int(2*mcum[-1]/mi)
                    print("Warning, boosting numbins to %i to increase mass resolution for interpolation" % numbins)
                else:
                    ok = True

            # Now go through and bisect for positions to get the mass per point we want.
            xi = xmin
            self.x = []
            self.rho = []
            mtarget = -0.5*mi
            while xi < xmax:
                mtarget += mi
                if mtarget <= mcum[-1]:
                    i = np.searchsorted(mcum, mtarget) - 1
                    assert mtarget >= mcum[i] and mtarget <= mcum[i+1]
                    xi = xmin + (i + (mtarget - mcum[i])/(mcum[i+1] - mcum[i]))*dx
                    assert (xi >= xmin + i*dx) and (xi <= xmin + (i+1)*dx)
                    self.x.append(xi)
                    self.rho.append(self.rhofunc(xi))
                else:
                    xi = xmax
            n = len(self.x)
            print("Generated %i 1D points." % n)
            self.m = [mi]*n

            # Figure out the H.
            self.H = []
            for i in range(n):
                if i == 0:
                    dxavg = self.x[i+1] - self.x[i]
                elif i == n-1:
                    dxavg = self.x[i] - self.x[i-1]
                else:
                    dxavg = 0.5*(self.x[i+1] - self.x[i-1])
                self.H.append(SymTensor1d(1.0/(nNodePerh*dxavg)))

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.m, self.rho, self.H)

        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        return Vector1d(self.x[i])

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
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

#-------------------------------------------------------------------------------
# Similarly generate a 1D profile in 2D along the x-direction.
#-------------------------------------------------------------------------------
class GeneratePlanarNodeProfile2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 nx,                 # target number of points in x
                 ny,                 # target number of points in y
                 rho,                # density profile, must be 1D function
                 xmin,               # (xmin, ymin) coordinates
                 xmax,               # (xmax, ymax) coordinates
                 nNodePerh = 2.01,
                 numbins = 10000,
                 SPH = False):

        assert nx > 0
        assert ny > 0
        assert xmin[0] < xmax[0]
        assert xmin[1] < xmax[1]
        assert nNodePerh > 0.0

        # First use the 1D generator to generate a 1D slice profile along x.
        gen1d = GenerateNodeProfile1d(nx = nx,
                                      rho = rho,
                                      xmin = xmin[0],
                                      xmax = xmax[0],
                                      nNodePerh = nNodePerh,
                                      numbins = numbins)

        # Stitch the 1D profiles back into serial data.
        gen1d.x = mpi.allreduce(gen1d.x, mpi.SUM)
        gen1d.m = mpi.allreduce(gen1d.m, mpi.SUM)
        gen1d.rho = mpi.allreduce(gen1d.rho, mpi.SUM)
        gen1d.H = mpi.allreduce(gen1d.H, mpi.SUM)
        n1d = len(gen1d.x)

        # Replicate the 1D slices into the full 2D data.
        self.x = []
        self.y = []
        self.m = []
        self.rho = []
        self.H = []
        dy = (xmax[1] - xmin[1])/ny
        hyinv = 1.0/(nNodePerh*dy)
        for iy in range(ny):
            self.x += gen1d.x
            self.y += [xmin[1] + (iy + 0.5)*dy]*n1d
            self.m += [mi*(xmax[1] - xmin[1])/ny for mi in gen1d.m]
            self.rho += gen1d.rho
            self.H += [SymTensor2d(H1d.xx, 0.0, 0.0, hyinv) for H1d in gen1d.H]

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.rho, self.H)

        # If we're forcing round H tensors, do it.
        if SPH:
            self.makeHround()

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
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

#-------------------------------------------------------------------------------
# Similarly generate a 1D profile in 3D along the x-direction.
#-------------------------------------------------------------------------------
class GeneratePlanarNodeProfile3d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 nx,                 # target number of points in x
                 ny,                 # target number of points in y
                 nz,                 # target number of points in z
                 rho,                # density profile, must be 1D function
                 xmin,               # (xmin, ymin, zmin) coordinates
                 xmax,               # (xmax, ymax, zmax) coordinates
                 nNodePerh = 2.01,
                 numbins = 10000,
                 SPH = False):

        assert nx > 0
        assert ny > 0
        assert nz > 0
        assert xmin[0] < xmax[0]
        assert xmin[1] < xmax[1]
        assert xmin[2] < xmax[2]
        assert nNodePerh > 0.0

        # First use the 1D generator to generate a 1D slice profile along x.
        gen1d = GenerateNodeProfile1d(nx = nx,
                                      rho = rho,
                                      xmin = xmin[0],
                                      xmax = xmax[0],
                                      nNodePerh = nNodePerh,
                                      numbins = numbins)

        # Stitch the 1D profiles back into serial data.
        gen1d.x = mpi.allreduce(gen1d.x, mpi.SUM)
        gen1d.m = mpi.allreduce(gen1d.m, mpi.SUM)
        gen1d.rho = mpi.allreduce(gen1d.rho, mpi.SUM)
        gen1d.H = mpi.allreduce(gen1d.H, mpi.SUM)
        n1d = len(gen1d.x)

        # Replicate the 1D slices into the full 3D data.
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.rho = []
        self.H = []
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz
        hyinv = 1.0/(nNodePerh*dy)
        hzinv = 1.0/(nNodePerh*dz)
        for iz in range(nz):
            for iy in range(ny):
                self.x += gen1d.x
                self.y += [xmin[1] + (iy + 0.5)*dy]*n1d
                self.z += [xmin[2] + (iz + 0.5)*dz]*n1d
                self.m += [mi*(xmax[1] - xmin[1])*(xmax[2] - xmin[2])/(ny*nz) for mi in gen1d.m]
                self.rho += gen1d.rho
                self.H += [SymTensor3d(H1d.xx, 0.0, 0.0,
                                       0.0, hyinv, 0.0,
                                       0.0, 0.0, hzinv) for H1d in gen1d.H]

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.rho, self.H)

        # If we're forcing round H tensors, do it.
        if SPH:
            self.makeHround()

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
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

