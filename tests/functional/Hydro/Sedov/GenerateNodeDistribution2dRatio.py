from math import *

from NodeGeneratorBase import *

from Spheral import Vector2d, Tensor2d, SymTensor2d, \
     rotationMatrix2d, testPointInBox2d

#-------------------------------------------------------------------------------
# Class to generate 2-D node positions.
#-------------------------------------------------------------------------------
class GenerateNodeDistribution2dRatio(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, nRadial, nTheta, rho,
                 xmin = None,
                 xmax = None,
                 rmin = None,
                 rmax = None,
                 nNodePerh = 2.01,
                 theta = pi/2.0,
                 azimuthalOffsetFraction = 0.0,
                 SPH = False,
                 rotation = 0.0,
                 offset = None,
                 xminreject = None,
                 xmaxreject = None,
                 rreject = None,
                 originreject = None,
                 reversereject = False,
                 relaxation = None,
                 rejecter = None,
                 xlmin = None,
                 xlmax = None):

        assert nRadial > 0
        assert nTheta > 0
        assert nNodePerh > 0.0
        assert offset is None or len(offset) == 2
        assert ((xminreject is None and xmaxreject is None) or
                (len(xminreject) == 2 and
                 len(xmaxreject) == 2))
        assert ((rreject is None and originreject is None) or
                (rreject and len(originreject) == 2))
        assert azimuthalOffsetFraction == 0.0 or theta == 2.0*pi

        self.nRadial = nRadial
        self.nTheta = nTheta
        self.xmin = xmin
        self.xmax = xmax
        self.rmin = rmin
        self.rmax = rmax
        self.nNodePerh = nNodePerh
        self.theta = theta
        self.azimuthalOffsetFraction = azimuthalOffsetFraction
        self.xlmin = xlmin
        self.xlmax = xlmax

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

 
        self.x, self.y, self.m, self.H = \
                self.latticeDistribution(self.nRadial, # nx
                                         self.nTheta,  # ny
                                         self.rho,
                                         self.xmin,
                                         self.xmax,
                                         self.rmin,
                                         self.rmax,
                                         self.nNodePerh,
                                         self.xlmin,
                                         self.xlmax)

        # If SPH has been specified, make sure the H tensors are round.
        if SPH:
            self.makeHround()

        # If requested, apply a rotation (in radians).
        if rotation:
            nhat = Vector2d(cos(rotation), -sin(rotation))
            R = rotationMatrix2d(nhat)
            for i in range(len(self.x)):
                v = R*Vector2d(self.x[i], self.y[i])
                self.x[i], self.y[i] = v.x, v.y
                self.H[i].rotationalTransform(R)

        # If requested, shift the nodes.
        if offset:
            for i in range(len(self.x)):
                self.x[i] += offset[0]
                self.y[i] += offset[1]

        # Did the user request to reject nodes based on a box?
        if xminreject:
            assert len(xminreject) == 2 and len(xmaxreject) == 2
            xminreject = Vector2d(*xminreject)
            xmaxreject = Vector2d(*xmaxreject)
            x = []
            y = []
            m = []
            H = []
            for xi, yi, mi, Hi in zip(self.x, self.y, self.m, self.H):
                t = testPointInBox2d(Vector2d(xi, yi), xminreject, xmaxreject)
                if ((t and (not reversereject)) or
                    ((not t) and reversereject)):
                    x.append(xi)
                    y.append(yi)
                    m.append(mi)
                    H.append(Hi)
            self.x = x
            self.y = y
            self.m = m
            self.H = H

        # Did the user request to reverse nodes based on a cylinder?
        if rreject:
            assert len(originreject) == 2
            originreject = Vector2d(*originreject)
            x = []
            y = []
            m = []
            H = []
            for xi, yi, mi, Hi in zip(self.x, self.y, self.m, self.H):
                t = (Vector2d(xi, yi) - originreject).magnitude() <= rreject
                if ((t and (not reversereject)) or
                    ((not t) and reversereject)):
                    x.append(xi)
                    y.append(yi)
                    m.append(mi)
                    H.append(Hi)
            self.x = x
            self.y = y
            self.m = m
            self.H = H

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.m, self.H = rejecter(self.x,
                                                      self.y,
                                                      self.m,
                                                      self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)

        # If requested, employ some relaxation on the NodeDistribution.
        if relaxation:
            relaxation(self)
        
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

    #---------------------------------------------------------------------------
    # Seed positions/masses on a lattice
    #---------------------------------------------------------------------------
    def latticeDistribution(self, nx, ny, rho,
                            xmin = (0.0, 0.0),
                            xmax = (1.0, 1.0),
                            rmin = None,
                            rmax = None,
                            nNodePerh = 2.01,
                            xlmin = 1e30,
                            xlmax = -1e30):

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny

        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        H0 = SymTensor2d(hx, 0.0, 0.0, hy)

        x = []
        y = []
        m = []
        H = []

        for j in range(ny):
            for i in range(nx):
                if i*dx < xlmin or i*dx > xlmax:
                    xx = xmin[0] + (i + 0.5)*dx
                    yy = xmin[1] + (j + 0.5)*dy
                    r = sqrt(xx*xx + yy*yy)
                    m0 = dx*dy*rho(Vector2d(xx, yy))
                    if ((r >= rmin or rmin is None) and
                        (r <= rmax or rmax is None)):
                        x.append(xx)
                        y.append(yy)
                        m.append(m0)
                        H.append(H0)
        dx = dx *0.5
        dy = dy *0.5
        xx = xlmin + 0.5*dx
        while (xx<=xlmax):
            for j in range(ny*2):
                yy = xmin[1] + (j + 0.5)*dy
                r = sqrt(xx*xx + yy*yy)
                m0 = dx*dy*rho(Vector2d(xx, yy))
                if ((r >= rmin or rmin is None) and
                        (r <= rmax or rmax is None)):
                        x.append(xx)
                        y.append(yy)
                        m.append(m0)
                        H.append(H0)
            xx += dx

        return x, y, m, H

    #------------------------------------------------------------------------------
    # Seed positions/masses on a cylindrical grid fading to a lattice at the edges
    #------------------------------------------------------------------------------
    def latticeCylindricalDistribution(self, nx, ny, rho,
                            xmin = (0.0, 0.0),
                            xmax = (1.0, 1.0),
                            rmin = None,
                            rmax = None,
                            nNodePerh = 2.01):
        k = 0
        np = 0
        deltar = rmax - rmin
        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        
        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        H0 = SymTensor2d(hx, 0.0, 0.0, hy)
        
        x = []
        y = []
        m = []
        H = []
        
        
        ml = []
        Hl = []
        xl = []
        yl = []
        
        xc = []
        yc = []
        mc = []
        Hc = []
        
        for j in range(ny):
            for i in range(nx):
                xx = xmin[0] + (i + 0.5)*dx
                yy = xmin[1] + (j + 0.5)*dy
                r = sqrt(xx*xx + yy*yy)
                m0 = dx*dy*rho(Vector2d(xx, yy))
                if (r>=rmin*0.8):
                    xl.append(xx)
                    yl.append(yy)
                    ml.append(m0)
                    Hl.append(H0)
                    k = k + 1
                if (r>=rmax):
                    x.append(xx)
                    y.append(yy)
                    m.append(m0)
                    H.append(H0)
                    np = np + 1
    
        # Start at the outermost radius, and work our way inward.
        theta = 2*3.14159
        ri = rmax+2.0*nNodePerh/nx
        
        #import random
        #random.seed()
        
        while ri > 0:
            
            # Get the nominal delta r, delta theta, number of nodes, and mass per
            # node at this radius.
            rhoi = rho(Vector2d(ri, 0.0))
            dr = sqrt(m0/rhoi)
            arclength = theta*ri
            arcmass = arclength*dr*rhoi
            nTheta = max(1, int(arcmass/m0))
            dTheta = theta/nTheta
            mi = arcmass/nTheta
            hi = nNodePerh*0.5*(dr + ri*dTheta)
            Hi = SymTensor2d(1.0/hi, 0.0,
                             0.0, 1.0/hi)

            # Now assign the nodes for this radius.
            for i in range(nTheta):
                thetai = (i + 0.5)*dTheta
                xc.append(ri*cos(thetai))
                yc.append(ri*sin(thetai))
                mc.append(mi)
                Hc.append(Hi)
                xi = ri*cos(thetai)
                yi = ri*sin(thetai)
                
                if(ri < rmin):
                    x.append(xi)
                    y.append(yi)
                    m.append(mi)
                    H.append(Hi)
                    np = np + 1
                elif(ri>=rmin):
                    #eps = random.random()
                    #func = ((ri-rmin)/deltar)**2
                    func = 1-ri/(rmin-rmax) - rmax/(rmax-rmin)
                    if(func>1.0):
                        func = 1.0
                    if(func<0.0):
                        func = 0.0
                    #if (eps <= func):
                    #x.append(ri*cos(thetai))
                    #y.append(ri*sin(thetai))
                    #m.append(mi)
                    #H.append(Hi)
                    #else:
                    minddr = nx
                    mink = 2*k
                    for j in range(k):
                        ddr = sqrt((xl[j]-xi)**2+(yl[j]-yi)**2)
                        if (ddr < minddr):
                            minddr = ddr
                            mink = j
                    xi = xi+(xl[mink]-xi)*func
                    yi = yi+(yl[mink]-yi)*func
                    
                    minddr = nx
                    for j in range(np):
                        ddr = sqrt((x[j]-xi)**2 + (y[j]-yi)**2)
                        if (ddr < minddr):
                            minddr = ddr
                    if(minddr > (1.0/hx)*0.5):
                        x.append(xi+(xl[mink]-xi)*func)
                        y.append(yi+(yl[mink]-yi)*func)
                        m.append(ml[mink])
                        H.append(Hl[mink])
        
     
            # Decrement to the next radial bin inward.
            ri = max(0.0, ri - dr)
    
        return x, y, m, H


    #---------------------------------------------------------------------------
    # Compute the area between offset circles in the given y range.
    #---------------------------------------------------------------------------
    def _areaBetweenOffsetCircles(self,
                                  r1, r2,
                                  c1, c2,
                                  y1, y2):
        assert r1 > 0.0
        assert r2 > 0.0
        dc = c2 - c1
        A1 = (y1*dc +
              0.5*(sqrt(r2*r2 - y1*y1) + r2*r2*asin(y1/r2)) -
              0.5*(sqrt(r1*r1 - y1*y1) + r1*r1*asin(y1/r1)))
        A2 = (y2*dc +
              0.5*(sqrt(r2*r2 - y2*y2) + r2*r2*asin(y2/r2)) -
              0.5*(sqrt(r1*r1 - y2*y2) + r1*r1*asin(y2/r1)))
        return A2 - A1
