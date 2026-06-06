#-------------------------------------------------------------------------------
# InteriorGenerator
#
# Specialized Spheral++ node generators for creating initial conditions inside
# some geometric boundary.
#-------------------------------------------------------------------------------
from math import *
from NodeGeneratorBase import *
from Spheral import Vector2d, Tensor2d, SymTensor2d

import random

#-------------------------------------------------------------------------------
# 2D
#-------------------------------------------------------------------------------
class InteriorGenerator2d(NodeGeneratorBase):

    def __init__(self, 
                 boundary,         # Some object that has "xmin", "xmax", & "contains"
                 dx,               # Nominal linear resolution
                 rho,              # initial mass density: constant, list, or function
                 nNodePerh = 2.01, # desired nPerh for H tensors
                 jitter = 0.0,     # (fraction of dx) any randomness to initial positions
                 SPH = False,      # Should we force round H tensors?
                 ):
        assert dx > 0.0
        assert nNodePerh > 0.0
        self.nNodePerh = nNodePerh

        # Start by clipping a lattice.
        xmin, xmax = boundary.xmin, boundary.xmax
        nx = max(1, int((xmax.x - xmin.x)/dx + 0.5))
        ny = max(1, int((xmax.y - xmin.y)/dx + 0.5))
        self.x, self.y = [], []
        for iy in range(ny):
            for ix in range(nx):
                posi = Vector2d(xmin.x + (ix + 0.5 + jitter*random.uniform(0,1))*dx,
                                xmin.y + (iy + 0.5 + jitter*random.uniform(0,1))*dx)
                if boundary.contains(posi):
                    self.x.append(posi.x)
                    self.y.append(posi.y)
        n = len(self.x)
        assert len(self.y) == n

        # Density and mass.
        if type(rho) is float:
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        # Mass per node.
        self.m = [self.rho(Vector2d(self.x[i], self.y[i])) for i in range(n)]

        # Set H.
        h0 = nNodePerh * dx
        H0 = SymTensor2d(1.0/h0, 0.0,
                         0.0, 1.0/h0)
        self.H = [H0]*len(self.x)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)

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
        return self.rho(Vector2d(self.x[i], self.y[i]))

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

