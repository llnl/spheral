from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double
from SpheralTestUtilities import *

import numpy as np

import mpi
procID = mpi.rank
nProcs = mpi.procs

#-------------------------------------------------------------------------------
# Class to generate 3-D node positions in a QVT fashion
# This is a direct port of Steen Hansen's QVT code. See Hansen et al. 2007
#-------------------------------------------------------------------------------
class QuaquaversalTiling3d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 n = 100,
                 xmin = 0.0,
                 xmax = 1.0,
                 rho = 1.0,
                 nNodePerh = 2.0,
                 offset=None,
                 rejecter=None,
                 maxLevel=6):

            self.maxLevel = maxLevel
            self.xmin = xmin
            self.xmax = xmax
            
            self.corners = []
            self.qua = []
            self.sph = []
            
            self.x = []
            self.y = []
            self.z = []
            self.m = []
            self.H = []
            
            self.mm = np.array([
                               [[1.0,0.0,0.0,0.0],
                                [0.5,0.5,0.0,0.0],
                                [0.5,0.0,0.5,0.0],
                                [0.5,0.0,0.0,0.5]],
                               [[0.0,1.0,0.0,0.0],
                                [0.5,0.5,0.0,0.0],
                                [0.5,0.5,-0.5,0.5],
                                [0.5,0.0,0.0,0.5]],
                               [[0.5, 0.5, -0.5, 0.5],
                                [0., 1., -0.5, 0.5],
                                [0., 1., 0., 0.],
                                [0., 0.5, 0.5, 0.]],
                               [[0.5, 0., 0.5, 0.],
                                [0., 0.5, 0.5, 0.],
                                [0., 0., 1., 0.],
                                [0., 0., 0.5, 0.5]],
                               [[1., 0., -0.5, 0.5],
                                [0.5, 0.5, -0.5, 0.5],
                                [0.5, 0., 0., 0.5],
                                [0.5, 0., -0.5, 1.]],
                               [[0., 1., -1., 1.],
                                [0.5, 0.5, -1., 1.],
                                [0.5, 0., -0.5, 1.],
                                [0.5, 0., 0., 0.5]],
                               [[0., 0., 0., 1.],
                                [0.25, 0.5, -0.75, 1.],
                                [0.5, 0., -0.5, 1.],
                                [0.5, 0., 0., 0.5]],
                               [[0., 0., 0.5, 0.5],
                                [0.25, 0.5, -0.25, 0.5],
                                [0., 1., -0.5, 0.5],
                                [0., 1., -1., 1.]]])
            self.fixpt = np.array([0.4285714, 0.2857143, 0.1428571])
            
            vec = np.array([0., - sqrt(3.), 0.,
                            0., 0., 0.,
                            0., 0., 1.,
                            -1., 0., 1.])
                            

            level = 0
            ii = 0

            n  = self.checkNorm(level,ii)
            dx = (xmax-xmin)
            vo = dx**3
            nd = n/vo
            m0 = rho/nd
            vi = 1.0/nd
            hi = 2.0*nNodePerh*pow(3.0/(4.0*pi)*vi,1.0/3.0)
            Hi = SymTensor3d(1.0/hi, 0.0, 0.0,
                             0.0, 1.0/hi, 0.0,
                             0.0, 0.0, 1.0/hi)
            
            
            vec = moveCenter(vec)
            #print vec
            
            vec = recurse(level,ii,vec)
            
            #print self.sph
            
            for i in xrange(len(self.sph)):
                self.x.append(self.scale(self.sph[i][0]))
                self.y.append(self.scale(self.sph[i][1]))
                self.z.append(self.scale(self.sph[i][2]))
                self.m.append(m0)
                self.H.append(Hi)
            

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
        loc = self.localPosition(i) - self.offset
        return self.densityProfileMethod(loc.magnitude())
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

    def checkNorm(self,ii):
        n = pow(8,self.maxLevel)
        print "This will produce %e points" % n

        for i in xrange(8):
            for j in xrange(4):
                b = 0
                for k in xrange(4):
                    b = b + self.mm[i][j][k]
                if (b!=1):
                    print "b = %f,%i,%i" %(b,i,j)
        return n

    def moveCenter(self,vec):
        for i in xrange(4):
            j = i*3
            print "%f %f %f" % (vec[j],vec[j+1],vec[j+2])

        for j in xrange(4):
            i = j*3
            vec[i] += (1.0-self.fixpt[0])
            vec[i+1] += sqrt(3.0)*self.fixpt[1]
            vec[i+2] += -self.fixpt[2]

        return vec

    def savePositions(self,vec):
        for i in xrange(4):
            self.corners.append([vec[i*3],vec[i*3+1],vec[i*3+2]])
        self.corners.append([vec[6],vec[7],vec[8]])
        self.corners.append([vec[0],vec[1],vec[2]])
        return

    def writeCenter(self,vec):
        cntr = np.zeros(3)
        for i in xrange(3):
            cntr[i] = vec[3+i] + (1.0-self.fixpt[0])*(vec[9+i]-vec[6+i])
            cntr[i] += self.fixpt[1]*(vec[i]-vec[3+i])
            cntr[i] += self.fixpt[2]*(vec[6+i]-vec[3+i])

        range = 0.14

        if ((abs(cntr[0]) < range) and (abs(cntr[1]) < range) and (abs(cntr[2]) < range)):
            self.qua.append([cntr[0],cntr[1],cntr[2]])

        radius = 0.0

        for i in xrange(3):
            radius += cntr[i]*cntr[i]
        radius = sqrt(radius)

        if (radius < range):
            self.sph.append([cntr[0],cntr[1],cntr[2]])
        return

    def recurse(self,level,ii,vec):
        vec2 = np.zeros(12)

        if (level == self.maxLevel):
            self.writeCenter(vec)
            return vec
        else:
            level += 1
            for ii in xrange(8):
                for i in xrange(12):
                    vec2[i] = 0
                for i in xrange(4):
                    for k in xrange(3):
                        for j in xrange(4):
                            vec2[k+3*i] += self.mm[ii][i][j] * vec[k+3*j]
                vec2 = recurse(level,ii,vec2)

            return vec2

    def scale(self,x):
        ymin = -0.14
        ymax = 0.14
        return (ymax-ymin)/(self.xmax-self.xmin)*(x-self.xmin) + ymin


