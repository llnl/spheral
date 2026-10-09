from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double
from SpheralTestUtilities import *
from SpheralTestUtilities import multiSort

import mpi
procID = mpi.rank
nProcs = mpi.procs

#-------------------------------------------------------------------------------
# Class to generate 3-D node positions in equal mass sheets in x-y plane
#-------------------------------------------------------------------------------
class GenerateEqualMassSheets3d(NodeGeneratorBase):
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, nx, ny, densityProfileMethod,
                 xmin = Vector3d(0,0,0),
                 xmax = Vector3d(1,1,1),
                 nNodePerh = 2.01,
                 rhoMin=None):
        assert nx > 0
        assert ny > 0
        assert nNodePerh > 0

        self.nx         = nx
        self.ny         = ny
        self.nNodePerh  = nNodePerh
        self.xmin       = xmin
        self.xmax       = xmax

        # no reason to support a constant density method here, just use a regular lattice for that
        self.densityProfileMethod = densityProfileMethod

        if rhoMin is None:
            rhoMin = densityProfileMethod.rhoMin    # hopefully your density method supports this!

        print(xmax)
        print(xmin)
        dz = float((xmax[0] - xmin[0])*(xmax[1]-xmin[1]))/(float(nx)*float(ny))
        self.m0 = rhoMin*dz*dz

        #print "dz,m0 = %f,%f" % (dz,self.m0)

        xl = []
        yl = []
        zl = []
        Hl = []
            
        zz = xmin[2]
        zc = 0
        rho = densityProfileMethod(zz)
        dx = dy = dz = sqrt(self.m0/rho)
        while (zz<xmax[2]):
            #print dz
            zz = xmin[2] + (zc+0.5)*dz
            rho = densityProfileMethod(zz)
            dx = dy = dz = sqrt(self.m0/rho)
            hx = 1.0/(nNodePerh*dx)
            hy = 1.0/(nNodePerh*dy)
            hz = 1.0/(nNodePerh*dz)
            H0 = SymTensor3d(hx, 0.0, 0.0,
                             0.0, hy, 0.0,
                             0.0, 0.0, hz)
            nxz = int((xmax[0]-xmin[0])/dx)
            nyz = int((xmax[1]-xmin[1])/dy)
            print("z = %f" % zz)
            for i in range(nxz):
                for j in range(nyz):
                    xx = xmin[0] + (i+0.5)*dx
                    yy = xmin[1] + (j+0.5)*dy
                    xl.append(xx)
                    yl.append(yy)
                    zl.append(zz)
                    Hl.append(H0)
            zc = zc + 1
        
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        
        for i in range(len(xl)):
            if self.within(Vector3d(xl[i],yl[i],zl[i]),xmin,xmax):
                self.x.append(xl[i])
                self.y.append(yl[i])
                self.z.append(zl[i])
                self.m.append(self.m0)
                self.H.append(Hl[i])
        del xl
        del yl
        del zl
        del Hl
            
        # Initialize the base class.  If "serialInitialization" is True, this
        # is where the points are broken up between processors as well.
        serialInitialization = True
        NodeGeneratorBase.__init__(self, serialInitialization,
                                   self.x, self.y, self.z, self.m, self.H)
            
        return


    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y) == len(self.z)
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
        loc = Vector3d(0,0,0)
        loc = self.localPosition(i)
        return self.densityProfileMethod(loc[2])
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

    #---------------------------------------------------------------------------
    # Helper function to determine within bounds
    #---------------------------------------------------------------------------
    def within(self,pos,xmin,xmax):
        a = 1
        for i in range(3):
            a = a*((pos[i]<xmax[i]) and (pos[i]>xmin[i]))
        return a
