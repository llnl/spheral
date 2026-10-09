from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double
from SpheralTestUtilities import *
#from Spheral import PairScalarFunctor, newtonRaphsonFindRoot
from SpheralTestUtilities import multiSort

import mpi
rank = mpi.rank
procs = mpi.procs

#-------------------------------------------------------------------------------
# Class to generate 3-D node positions in a stretched lattice
#-------------------------------------------------------------------------------
class GenerateStretchedLattice3d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, nr, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01,
                 offset=None,
                 m0ForMassMatching=None):
        
        assert nr > 0
        assert rmin >= 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0
        assert offset is None or len(offset)==3
        
        if offset is None:
            self.offset = Vector3d(0,0,0)
        else:
            self.offset = Vector3d(offset[0],offset[1],offset[2])
        
        self.nr         = nr
        self.rmin       = rmin
        self.rmax       = rmax
        self.thetaMin   = thetaMin
        self.thetaMax   = thetaMax
        self.phiMin     = phiMin
        self.phiMax     = phiMax
        self.nNodePerh  = nNodePerh
        
        self.xmin       = Vector3d(-2.0*rmax,-2.0*rmax,-2.0*rmax)
        self.xmax       = Vector3d(2.0*rmax,2.0*rmax,2.0*rmax)
        
        # no reason to support a constant density method here, just use a regular lattice for that
        self.densityProfileMethod = densityProfileMethod
        
        # Determine how much total mass there is in the system.
        targetMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 thetaMin, thetaMax,
                                                 phiMin, phiMax)

        
        #targetMass = self.integrateTotalMass(self.densityProfileMethod,
        #                                         rmax)
        
        targetN         = 4.0/3.0*pi*(nr**3)
        self.m0         = targetMass/targetN
        self.vol        = 4.0/3.0*pi*(rmax**3)
        # what this means is this currently only supports creating a full sphere and then
        # cutting out the middle to rmin if rmin > 0
        
        if m0ForMassMatching is None:
            self.rho0       = targetMass/self.vol
        else:
            self.m0 = m0ForMassMatching
            self.rho0 = targetN*self.m0/self.vol
        
        print("Found total mass = {0:3.3e} with rho0 = {1:3.3e}".format(targetMass,self.rho0))
    
        # compute kappa first
        # k = 3/(self.rho0*rmax**3) * targetMass/(4.0*pi)
        # print "Found kappa={0:3.3f}. Was that what you expected?".format(k)
        
        nlat = nr
        
        # create the unstretched lattice
        self.xl, self.yl, self.zl, self.ml, self.Hl = \
            self.latticeDistribution(nlat,
                                     self.rho0,
                                     self.m0,
                                     self.xmin,    # (xmin, ymin, zmin)
                                     self.xmax,    # (xmax, ymax, zmax)
                                     self.rmax,
                                     self.nNodePerh)


        self.rl = []
        for i in range(len(self.xl)):
            self.rl.append(sqrt(self.xl[i]**2+self.yl[i]**2+self.zl[i]**2))
        
        print("Sorting unstretched lattice... %d elements" % len(self.rl))
        
        multiSort(self.rl,self.xl,self.yl,self.zl)
        
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        
        nx  = 2*nlat+1
        eta = (self.xmax[0] - self.xmin[0])/nx
        
        print("Stretching lattice...")

        dr  = eta * 0.01    # this will essentially be the error in the new dumb way
        r0p = 0
        rp  = 0
        rn  = 0
        npp = len(self.rl)/procs
        rankmin = (npp*rank) if (rank>0) else (1)
        rankmax = (npp*(rank+1)) if (rank<procs-1) else (len(self.rl))
        for i in range(rankmin,rankmax):
            #print "%d / %d" % (i,len(self.rl))
            r0 = self.rl[i]
            if (abs(r0-r0p)/r0>1e-10):
                sol     = r0**3*self.rho0/3.0
                iter    = int(10*rmax // dr)
                fn      = 0
                for j in range(iter+1):
                    rj  = dr*j
                    rjj = dr*(j+1)
                    fj  = rj**2 * densityProfileMethod(rj)
                    fjj = rjj**2 * densityProfileMethod(rjj)
                    fn  = fn + 0.5*(fj+fjj)*dr
                    if (fn>=sol):
                        rn = rj
                        break
            r0p = r0
            if (rn <= rmax and rn > rmin):
                self.x.append(self.xl[i] * rn/r0)
                self.y.append(self.yl[i] * rn/r0)
                self.z.append(self.zl[i] * rn/r0)
                self.m.append(self.ml[i])
                self.H.append(self.Hl[i])
    
        seededMass = sum(self.m)
        seededMass = mpi.allreduce(seededMass,mpi.SUM)
        
        mAdj = targetMass / seededMass
        for i in range(len(self.m)):
            self.m[i] = self.m[i] * mAdj
        
        '''
        rp  = 0
        r0p = 0
        rn  = 0
        for i in xrange(1,len(self.rl)+1):
            if (self.rl[i] <= rmax):
                r0 = self.rl[i]
                if (abs(r0-r0p)/r0>1e-10):
                    eta = r0-r0p
                    print "Found eta={0:3.3e} @ r0,rp = {1:3.3e},{2:3.3e}".format(eta,r0,rp)
                    stretch = self.nrFunc(k,self.rho0,r0,eta,rp,self.densityProfileMethod,eta/100.0)
                    # want to do a check here to determine the range based on slope of function
                    lrmax = rmax
                    if self.dfunc(r0,eta,rp,eta/100.0,self.densityProfileMethod) > 0:
                        lrmax = r0
                    rn = newtonRaphsonFindRoot(stretch,0,lrmax,1.0e-15, 1.0e-15)
                rp  = rn
                r0p = r0
        
        
        for ix in xrange(nr):
            idx = (2*nr+1)*nr*(2*nr+2)
            i   = ix + nr + 1 + idx
            mi  = i - 2*ix
            j   = idx + (2*nr+1)*ix
            mj  = idx - (2*nr+1)*ix
            k   = idx + ix*(2*nr+1)**2
            mk  = idx - ix*(2*nr+1)**2
            print i
            r0 = sqrt(self.x[i]*self.x[i] + self.y[i]*self.y[i] + self.z[i]*self.z[i])
            rn = r0
            dr = eta/10.0
            
            myFunc = myClass()
            #rn = self.rootFind(rn,k,self.rho0,r0,eta,rp,dr,r0,self.densityProfileMethod,ns)
            #self.x[i] = self.x[i] * rn/r0
            #self.y[i] = self.y[i] * rn/r0
            #self.z[i] = self.z[i] * rn/r0
            #rp = rn
        
        '''
        
        # Initialize the base class.  If "serialInitialization" is True, this
        # is where the points are broken up between processors as well.
        serialInitialization = False
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
        loc = self.localPosition(i) - self.offset
        return self.densityProfileMethod(loc.magnitude())
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]


    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           rmin, rmax,
                           thetaMin, thetaMax,
                           phiMin, phiMax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nbins
        for i in range(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        return result

    #-------------------------------------------------------------------------------
    # Seed positions/masses on the unstretched lattice
    #-------------------------------------------------------------------------------
    def latticeDistribution(self, nr, rho0, m0,
                            xmin,
                            xmax,
                            rmax,
                            nNodePerh = 2.01):
        
        assert nr > 0
        assert rho0 > 0
        
        nx = 2*nr+1
        ny = 2*nr+1
        nz = 2*nr+1
        
        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz
        
        n = nx*ny*nz

        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        hz = 1.0/(nNodePerh*dz)
        H0 = SymTensor3d(hx, 0.0, 0.0,
                         0.0, hy, 0.0,
                         0.0, 0.0, hz)
                       
        x = []
        y = []
        z = []
        m = []
        H = []

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    xx = xmin[0] + (i + 0.5)*dx
                    yy = xmin[1] + (j + 0.5)*dy
                    zz = xmin[2] + (k + 0.5)*dz
                    x.append(xx)
                    y.append(yy)
                    z.append(zz)
                    m.append(m0)
                    H.append(H0)
       
        return x, y, z, m, H
