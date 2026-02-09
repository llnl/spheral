from math import *

from NodeGeneratorBase import *

from Spheral import Vector2d, Tensor2d, SymTensor2d, rotationMatrix2d, testPointInBox2d, toroidalVolume
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Class to generate 2-D node positions.
#-------------------------------------------------------------------------------
class GenerateNodeDistribution2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, nRadial, nTheta, rho,
                 distributionType = "optimal",
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
                 rejecter = None):

        assert nRadial > 0
        assert nTheta > 0
        assert (((distributionType == "optimal" or
                  distributionType == "constantDTheta" or
                  distributionType == "constantNTheta" or
                  distributionType == "powerOf2NTheta") and
                 (rmin is not None and rmax is not None and
                  rmin < rmax and
                  theta is not None and theta > 0.0)) or
                ((distributionType == "lattice" or
                  distributionType == "xstaggeredLattice" or
                  distributionType == "rotatedLattice") and
                 xmin is not None and xmax is not None and
                 xmin[0] < xmax[0] and xmin[1] < xmax[1]) or
                (distributionType == "offsetCylindrical" and
                 xmin is not None and xmax is not None and
                 rmin is not None and rmax is not None) or
                (distributionType == "latticeCylindrical" and
                 xmin is not None and xmax is not None and
                 rmin is not None and rmax is not None))
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
        self.distributionType = distributionType
        self.xmin = xmin
        self.xmax = xmax
        self.rmin = rmin
        self.rmax = rmax
        self.nNodePerh = nNodePerh
        self.theta = theta
        self.azimuthalOffsetFraction = azimuthalOffsetFraction

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        # Generate the internal sets of positions, masses, and H
        if distributionType == "optimal":
            self.x, self.y, self.m, self.H = \
                    self.optimalCylindricalDistribution(self.nRadial,
                                                        self.nTheta,
                                                        self.rho,
                                                        self.xmin,
                                                        self.xmax,
                                                        self.rmin,
                                                        self.rmax,
                                                        self.nNodePerh,
                                                        self.theta,
                                                        self.azimuthalOffsetFraction)
        elif distributionType == "constantDTheta":
            self.x, self.y, self.m, self.H = \
                    self.constantDThetaCylindricalDistribution(self.nRadial,
                                                               self.rho,
                                                               self.xmin,
                                                               self.xmax,
                                                               self.rmin,
                                                               self.rmax,
                                                               self.nNodePerh,
                                                               self.theta,
                                                               self.azimuthalOffsetFraction)
        elif distributionType == "constantNTheta":
            self.x, self.y, self.m, self.H = \
                    self.constantNThetaCylindricalDistribution(self.nRadial,
                                                               self.nTheta,
                                                               self.rho,
                                                               self.xmin,
                                                               self.xmax,
                                                               self.rmin,
                                                               self.rmax,
                                                               self.nNodePerh,
                                                               self.theta,
                                                               self.azimuthalOffsetFraction)

        elif distributionType == "powerOf2NTheta":
            self.x, self.y, self.m, self.H = \
                    self.powerOf2NThetaCylindricalDistribution(self.nRadial,
                                                               self.rho,
                                                               self.xmin,
                                                               self.xmax,
                                                               self.rmin,
                                                               self.rmax,
                                                               self.nNodePerh,
                                                               self.theta,
                                                               self.azimuthalOffsetFraction)

        elif distributionType == "lattice":
            self.x, self.y, self.m, self.H = \
                    self.latticeDistribution(self.nRadial, # nx
                                             self.nTheta,  # ny
                                             self.rho,
                                             self.xmin,
                                             self.xmax,
                                             self.rmin,
                                             self.rmax,
                                             self.nNodePerh)

        elif distributionType == "xstaggeredLattice":
            self.x, self.y, self.m, self.H = \
                    self.xstaggeredLatticeDistribution(self.nRadial, # nx
                                                       self.nTheta,  # ny
                                                       self.rho,
                                                       self.xmin,
                                                       self.xmax,
                                                       self.rmin,
                                                       self.rmax,
                                                       self.nNodePerh)

        elif distributionType == "rotatedLattice":
            self.x, self.y, self.m, self.H = \
                    self.rotatedLatticeDistribution(self.nRadial, # nx
                                                    self.nTheta,  # ny
                                                    self.rho,
                                                    self.xmin,
                                                    self.xmax,
                                                    self.rmin,
                                                    self.rmax,
                                                    self.nNodePerh)

        elif distributionType == "offsetCylindrical":
            self.x, self.y, self.m, self.H = \
                    self.offsetCylindricalDistribution(self.nRadial,
                                                       self.nTheta,
                                                       self.rho,
                                                       self.xmin,
                                                       self.xmax,
                                                       self.rmin,
                                                       self.rmax,
                                                       self.nNodePerh)
        
        elif distributionType == "latticeCylindrical":
            self.x, self.y, self.m, self.H = \
                    self.latticeCylindricalDistribution(self.nRadial,
                                                        self.nTheta,
                                                        self.rho,
                                                        self.xmin,
                                                        self.xmax,
                                                        self.rmin,
                                                        self.rmax,
                                                        self.nNodePerh)

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
    # This version tries to combine the optimal configuration of constant
    # spacing for the innermost node distribution, and constant nodes per ring
    # for the outer regions.
    #---------------------------------------------------------------------------
    def optimalCylindricalDistribution(self, nRadial, nTheta, rho,
                                       xmin = None,
                                       xmax = None,
                                       rmin = 0.0,
                                       rmax = 1.0,
                                       nNodePerh = 2.01,
                                       theta = pi/2.0,
                                       azimuthalOffsetFraction = 0.0):

        # Determine the radius at which we want to shift from constant spacing
        # to constant nodes per radius.
        dr = (rmax - rmin)/nRadial
        dTheta = theta/nTheta

        r0 = rmin + 2**(log(nTheta*theta/(2.0*pi))/log(2))*dr/theta
        #r0 = rmin + dr/dTheta

        print("optimalCylindricalDistribution: Cutoff radius is ", r0)
        nRadial0 = max(0, min(nRadial, int((r0 - rmin)/dr)))
        nRadial1 = nRadial - nRadial0
        r0 = rmin + nRadial0*dr
        print("Shifted to ", r0)
        print(nRadial0, nRadial1)

        if nRadial0 and not nRadial1:
            # Only use constant spacing for the nodes.
            x, y, m, H = self.constantDThetaCylindricalDistribution(nRadial0, rho,
                                                                    xmin,
                                                                    xmax,
                                                                    rmin,
                                                                    rmax,
                                                                    nNodePerh,
                                                                    theta,
                                                                    azimuthalOffsetFraction)

        elif not nRadial0 and nRadial1:
            # Only use constant nodes per radial bin
            x, y, m, H = self.constantNThetaCylindricalDistribution(nRadial1,
                                                                    nTheta,
                                                                    rho,
                                                                    xmin,
                                                                    xmax,
                                                                    rmin,
                                                                    rmax,
                                                                    nNodePerh,
                                                                    theta,azimuthalOffsetFraction)

        else:
            # Combine the two schemes.
            x0, y0, m0, H0 = self.constantDThetaCylindricalDistribution(nRadial0, rho,
                                                                        xmin,
                                                                        xmax,
                                                                        rmin,
                                                                        r0,
                                                                        nNodePerh,
                                                                        theta,
                                                                        azimuthalOffsetFraction)
            x1, y1, m1, H1 = self.constantNThetaCylindricalDistribution(nRadial1, nTheta, rho,
                                                                        xmin,
                                                                        xmax,
                                                                        r0,
                                                                        rmax,
                                                                        nNodePerh,
                                                                        theta,
                                                                        azimuthalOffsetFraction)
            x, y, m, H = x0 + x1, y0 + y1, m0 + m1, H0 + H1

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry
    # This version tries to seed nodes in circular rings with constant spacing.
    #---------------------------------------------------------------------------
    def constantDThetaCylindricalDistribution(self, nRadial, rho,
                                              xmin = None,
                                              xmax = None,
                                              rmin = 0.0,
                                              rmax = 1.0,
                                              nNodePerh = 2.01,
                                              theta = pi/2.0,
                                              azimuthalOffsetFraction = 0.0):

        from Spheral import SymTensor2d

        dr = (rmax - rmin)/nRadial
        x = []
        y = []
        m = []
        H = []
        
        h = 1.0/(nNodePerh*dr)
        Hi = SymTensor2d(h, 0.0, 0.0, h)

        for i in range(0, nRadial):
            rInner = rmin + i*dr
            rOuter = rmin + (i + 1)*dr
            ri = rmin + (i + 0.5)*dr
            li = theta*ri
            #nominalNTheta = 2**int(log(li/dr)/log(2) + 0.5)*2.0*pi/theta
            nominalNTheta = int(li/dr)
            nTheta = max(1, int(nominalNTheta))
            dTheta = theta/nTheta
            mRing = (rOuter**2 - rInner**2) * theta/2.0 * rho(Vector2d(ri, 0.0))
            mi = mRing/nTheta
            for j in range(nTheta):
                thetai = (j + 0.5 + i*azimuthalOffsetFraction)*dTheta
                xi = ri*cos(thetai)
                yi = ri*sin(thetai)
                use = True
                if xmin:
                    if xi < xmin[0] or yi < xmin[1]:
                        use = False
                if xmax:
                    if xi > xmax[0] or yi > xmax[1]:
                        use = False
                if use:
                    x.append(xi)
                    y.append(yi)
                    m.append(mi)
                    H.append(Hi)

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry.
    # This version seeds a constant number of nodes per radial bin.
    #---------------------------------------------------------------------------
    def constantNThetaCylindricalDistribution(self, nRadial, nTheta, rho,
                                              xmin = None,
                                              xmax = None,
                                              rmin = 0.0,
                                              rmax = 1.0,
                                              nNodePerh = 2.01,
                                              theta = pi/2.0,
                                              azimuthalOffsetFraction = 0.0):

        from Spheral import Tensor2d
        from Spheral import SymTensor2d

        dr = (rmax - rmin)/nRadial
        dTheta = theta/nTheta
        hr = 1.0/(nNodePerh*dr)

        x = []
        y = []
        m = []
        H = []

        for i in range(0, nRadial):
            rInner = rmin + i*dr
            rOuter = rmin + (i + 1)*dr
            ri = rmin + (i + 0.5)*dr

            mRing = (rOuter**2 - rInner**2) * theta/2.0 * rho(Vector2d(ri, 0.0))
            mi = mRing/nTheta

            hTheta = 1.0/(nNodePerh*ri*dTheta)

            for j in range(nTheta):
                thetai = (j + 0.5 + i*azimuthalOffsetFraction)*dTheta
                xi = ri*cos(thetai)
                yi = ri*sin(thetai)
                use = True
                if xmin:
                    if xi < xmin[0] or yi < xmin[1]:
                        use = False
                if xmax:
                    if xi > xmax[0] or yi > xmax[1]:
                        use = False
                if use:
                    x.append(xi)
                    y.append(yi)
                    m.append(mi)

                    Hi = SymTensor2d(hr, 0.0, 0.0, hTheta)
                    rot = Tensor2d(cos(thetai), sin(thetai), -sin(thetai), cos(thetai))
                    Hii = Hi
                    Hii.rotationalTransform(rot)
                    H.append(Hii)

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry
    # This is a modified form of constantDTheta that restricts rings to jump by
    # powers of 2.
    #---------------------------------------------------------------------------
    def powerOf2NThetaCylindricalDistribution(self, nRadial, rho,
                                              xmin = None,
                                              xmax = None,
                                              rmin = 0.0,
                                              rmax = 1.0,
                                              nNodePerh = 2.01,
                                              theta = pi/2.0,
                                              azimuthalOffsetFraction = 0.0):

        from Spheral import SymTensor2d

        dr = (rmax - rmin)/nRadial
        x = []
        y = []
        m = []
        H = []
        
        h = 1.0/(nNodePerh*dr)
        Hi = SymTensor2d(h, 0.0, 0.0, h)

        # Figure out the innermost rings nTheta
        rInner = rmin
        rOuter = rmin + dr
        ri = rmin + 0.5*dr
        li = theta*ri
        baseNTheta = max(1, int(li/dr + 0.5))

        for i in range(0, nRadial):
            rInner = rmin + i*dr
            rOuter = rmin + (i + 1)*dr
            ri = rmin + (i + 0.5)*dr
            li = theta*ri
            nominalNTheta = max(1, int(li/dr))
            nTheta = max(1, baseNTheta * 2**max(0, int(log(float(nominalNTheta)/float(baseNTheta))/log(2.0))))
            dTheta = theta/nTheta
            mRing = (rOuter**2 - rInner**2) * theta/2.0 * rho(Vector2d(ri, 0.0))
            mi = mRing/nTheta
            for j in range(nTheta):
                thetai = (j + 0.5 + i*azimuthalOffsetFraction)*dTheta
                xi = ri*cos(thetai)
                yi = ri*sin(thetai)
                use = True
                if xmin:
                    if xi < xmin[0] or yi < xmin[1]:
                        use = False
                if xmax:
                    if xi > xmax[0] or yi > xmax[1]:
                        use = False
                if use:
                    x.append(xi)
                    y.append(yi)
                    m.append(mi)
                    H.append(Hi)

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions between circles of different radii and centers!
    #---------------------------------------------------------------------------
    def offsetCylindricalDistribution(self, nr, ny, rho,
                                      xmin, xmax,
                                      rmin, rmax,
                                      nNodePerh = 2.01):
        c1 = xmin[0]
        c2 = xmax[0]
        y1 = xmin[1]
        y2 = xmax[1]
        A = self._areaBetweenOffsetCircles(rmin, rmax, c1, c2, y1, y2)
        y12 = 0.5*(xmin[1] + xmax[1])
        dx12 = sqrt(rmax*rmax - y12*y12) - sqrt(rmin*rmin - y12*y12) + c2 - c1
        if ny == 0:
            dy = dx12/nr
            ny = int((y2 - y1)/dy + 0.5)
        dy = (y2 - y1)/ny
        hy = 1.0/(nNodePerh*dy)

        x = []
        y = []
        m = []
        H = []

        for iy in range(ny):
            yi = y1 + (iy + 0.5)*dy
            x1 = c1 + sqrt(rmin*rmin - yi*yi)
            x2 = c2 + sqrt(rmax*rmax - yi*yi)
            assert x2 > x1
            dx = (x2 - x1)/nr
            hx = 1.0/(nNodePerh*dx)
            Hi = SymTensor2d(hx, 0.0, 0.0, hy)
            volrow = self._areaBetweenOffsetCircles(rmin, rmax, c1, c2,
                                                    y1 + iy*dy,
                                                    y1 + (iy + 1)*dy)
            mi = volrow*rho(Vector2d(0.5*(x1 + x2), yi))/nr
            for ix in range(nr):
                xi = x1 + (ix + 0.5)*dx
                assert xi > x1 and xi < x2
                x.append(xi)
                y.append(yi)
                m.append(mi)
                H.append(Hi)

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions/masses on a lattice
    #---------------------------------------------------------------------------
    def latticeDistribution(self, nx, ny, rho,
                            xmin = (0.0, 0.0),
                            xmax = (1.0, 1.0),
                            rmin = None,
                            rmax = None,
                            nNodePerh = 2.01):

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
                xx = xmin[0] + (i + 0.5)*dx
                yy = xmin[1] + (j + 0.5)*dy
                r = sqrt(xx*xx + yy*yy)
                m0 = dx*dy*rho(Vector2d(xx, yy))
                if ((rmin is None or r >= rmin) and
                    (rmax is None or r <= rmax)):
                    x.append(xx)
                    y.append(yy)
                    m.append(m0)
                    H.append(H0)

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
    # Seed positions on a staggered (or offset) lattice
    #---------------------------------------------------------------------------
    def xstaggeredLatticeDistribution(self, nx, ny, rho,
                                      xmin = (0.0, 0.0),
                                      xmax = (1.0, 1.0),
                                      rmin = None,
                                      rmax = None,
                                      nNodePerh = 2.01):

        # Generate the basic distribution from the lattice option.
        x, y, m, H = self.latticeDistribution(nx, ny, rho,
                                              xmin, xmax,
                                              rmin, rmax,
                                              nNodePerh)

        # Now go through and stagger the positions in x.
        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        ddx = 0.25*dx

        assert len(x) == len(y)
        assert dy > 0.0
        for i in range(len(x)):
            x[i] += ddx * (-1.0)**int((y[i] - xmin[1])/dy + 0.1)

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions on a rotated lattice.
    #---------------------------------------------------------------------------
    def rotatedLatticeDistribution(self, nx, ny, rho,
                                   xmin = (0.0, 0.0),
                                   xmax = (1.0, 1.0),
                                   rmin = None,
                                   rmax = None,
                                   nNodePerh = 2.01):

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny

        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        H0 = SymTensor2d(hx, 0.0, 0.0, hy)

        x = []
        y = []
        m = []
        H = []

        for j in range(ny + 1):
            jmod = (j + 1) % 2
            nxrow = nx + jmod
            for i in range(nxrow):
                xx = xmin[0] + (i + 0.5 - 0.5*jmod)*dx
                yy = xmin[1] + j*dy
                r = sqrt(xx*xx + yy*yy)
                m0 = dx*dy*rho(Vector2d(xx, yy))
                if j == 0 or j == ny:
                    m0 /= 2.0
                if jmod == 1 and (i == 0 or i == nxrow - 1):
                    m0 /= 2.0
                if ((rmin is None or r >= rmin) and
                    (rmax is None or r <= rmax)):
                    x.append(xx)
                    y.append(yy)
                    m.append(m0)
                    H.append(H0)

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

#-------------------------------------------------------------------------------
# Specialized version that generates a variable radial stepping to try
# and match a given density profile with (nearly) constant mass nodes.  This
# only supports the equivalent of the constant Dtheta method.
#-------------------------------------------------------------------------------
class GenerateNodesMatchingProfile2d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, n, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi/2.0,
                 nNodePerh = 2.01,
                 offset = [0,0],
                 m0 = 0.0):
        
        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert nNodePerh > 0.0
        assert m0 >= 0.0
        
        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.thetaMin = thetaMin
        self.thetaMax = thetaMax
        self.nNodePerh = nNodePerh
        self.offset = offset
        
        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(densityProfileMethod) == type(1.0):
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.densityProfileMethod = densityProfileMethod
        
        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 thetaMin, thetaMax)
        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax))

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(self.n*self.n*pi)
        if (m0 > 0.0):
            self.m0 = m0
        assert self.m0 > 0.0
        print("Nominal mass per node of %g." % self.m0)

        # OK, we now know enough to generate the node positions.
        self.x, self.y, self.m, self.H = \
            self.constantDThetaCylindricalDistribution(self.densityProfileMethod,
                                                       self.m0,
                                                       rmin, rmax,
                                                       thetaMin, thetaMax,
                                                       nNodePerh)
        print("Generated a total of %i nodes." % len(self.x))

        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0
        massCorrection = self.totalMass/sumMass
        for i in range(len(self.m)):
            self.m[i] *= massCorrection
        print("Applied a mass correction of %f to ensure total mass is %f." % (
                                                                               massCorrection, self.totalMass))

        # Have the base class break up the serial node distribution
        # for parallel cases.
        
        for i in range(len(self.m)):
            self.x[i] = self.x[i] + self.offset[0]
            self.y[i] = self.y[i] + self.offset[1]
        
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
        return self.densityProfileMethod(sqrt(pow(self.localPosition(i)[0]-self.offset[0],2.0)+pow(self.localPosition(i)[1]-self.offset[1],2.0)))
    
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
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        theta = thetaMax - thetaMin
        h = (rmax - rmin)/nbins
        result = (rmin*densityProfileMethod(rmin) +
                  rmax*densityProfileMethod(rmax))
        for i in range(1, nbins):
            ri = rmin + i*h
            if i % 2 == 0:
                result += 4.0*ri*densityProfileMethod(ri)
            else:
                result += 2.0*ri*densityProfileMethod(ri)

        result *= theta*h/3.0
        return result
    
    #---------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry
    # This version tries to seed nodes in circular rings with constant spacing.
    #---------------------------------------------------------------------------
    def constantDThetaCylindricalDistribution(self, densityProfileMethod,
                                              m0,
                                              rmin, rmax,
                                              thetaMin, thetaMax,
                                              nNodePerh = 2.01):
        
        from Spheral import SymTensor2d
        
        # Return lists for positions, masses, and H's.
        x = []
        y = []
        m = []
        H = []
        
        # Start at the outermost radius, and work our way inward.
        theta = thetaMax - thetaMin
        ri = rmax
        while ri > rmin:

            
            
            # Get the nominal delta r, delta theta, number of nodes, and mass per
            # node at this radius.
            rhoi = densityProfileMethod(ri)
            dr = sqrt(m0/rhoi)
            rii = ri+dr/2.0
            rhoi = densityProfileMethod(rii)
            arclength = theta*rii
            arcmass = arclength*dr*rhoi
            nTheta = max(1, int(arcmass/m0))
            dTheta = theta/nTheta
            mi = arcmass/nTheta
            hi = nNodePerh*0.5*(dr + rii*dTheta)
            Hi = SymTensor2d(1.0/hi, 0.0,
                             0.0, 1.0/hi)
                             
            # Now assign the nodes for this radius.
            for i in range(nTheta):
                thetai = thetaMin + (i + 0.5)*dTheta
                x.append(rii*cos(thetai))
                y.append(rii*sin(thetai))
                m.append(mi)
                H.append(Hi)

            # Decrement to the next radial bin inward.
            ri = max(0.0, ri - dr)
        
        return x, y, m, H




#-------------------------------------------------------------------------------
# Specialized version that generates a variable mass distribution to try
# and match a given mass function with a constant density.  This
# only supports the equivalent of the constant Dtheta method.
#-------------------------------------------------------------------------------
class GenerateNodesMatchingMassProfile2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, n, densityProfileMethod,
                 massProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi/2.0,
                 nNodePerh = 2.01):

        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert nNodePerh > 0.0

        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.thetaMin = thetaMin
        self.thetaMax = thetaMax
        self.nNodePerh = nNodePerh

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        # CR: i'm currently going to assume this is used in constantRho form
        # every time. i may come back to this and make it work with both profiles
        # later.
        if type(densityProfileMethod) == type(1.0):
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.densityProfileMethod = densityProfileMethod
        
        self.massProfileMethod = massProfileMethod

        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 thetaMin, thetaMax)
        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g)" % \
              (self.totalMass, rmin, rmax, thetaMin, thetaMax))

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/self.n
        assert self.m0 > 0.0
        print("Nominal mass per node of %g." % self.m0)

        # OK, we now know enough to generate the node positions.
        self.x, self.y, self.m, self.H = \
                self.constantDThetaCylindricalDistribution(self.densityProfileMethod,
                                                           self.massProfileMethod,
                                                           rmin, rmax,
                                                           thetaMin, thetaMax,
                                                           nNodePerh)
        print("Generated a total of %i nodes." % len(self.x))

        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0
        massCorrection = self.totalMass/sumMass
        for i in range(len(self.m)):
            self.m[i] *= massCorrection
        print("Applied a mass correction of %f to ensure total mass is %f." % (
            massCorrection, self.totalMass))

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
        return self.densityProfileMethod(self.localPosition(i).magnitude())

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
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0

        theta = thetaMax - thetaMin
        h = (rmax - rmin)/nbins
        result = (rmin*densityProfileMethod(rmin) +
                  rmax*densityProfileMethod(rmax))
        for i in range(1, nbins):
            ri = rmin + i*h
            if i % 2 == 0:
                result += 4.0*ri*densityProfileMethod(ri)
            else:
                result += 2.0*ri*densityProfileMethod(ri)

        result *= theta*h/3.0
        return result

    #---------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry
    # This version tries to seed nodes in circular rings with constant spacing.
    #---------------------------------------------------------------------------
    def constantDThetaCylindricalDistribution(self, densityProfileMethod,
                                              massProfileMethod,
                                              rmin, rmax,
                                              thetaMin, thetaMax,
                                              nNodePerh = 2.01):

        from Spheral import SymTensor2d

        # Return lists for positions, masses, and H's.
        x = []
        y = []
        m = []
        H = []
        
        # Start at the outermost radius, and work our way inward.
        theta = thetaMax - thetaMin
        ri = rmax
        while ri > rmin:

            # Get the nominal delta r, delta theta, number of nodes, and mass per
            # node at this radius.
            rhoi = densityProfileMethod(ri)
            mpi = massProfileMethod(ri)
            dr = sqrt(mpi/rhoi)
            arclength = theta*ri
            arcmass = arclength*dr*rhoi
            nTheta = max(1, int(arcmass/mpi))
            dTheta = theta/nTheta
            mi = arcmass/nTheta
            hi = nNodePerh*0.5*(dr + ri*dTheta)
            Hi = SymTensor2d(1.0/hi, 0.0,
                             0.0, 1.0/hi)

            # Now assign the nodes for this radius.
            for i in range(nTheta):
                if nTheta > 1:
                    thetai = thetaMin + (i + 0.5)*dTheta
                    x.append(ri*cos(thetai))
                    y.append(ri*sin(thetai))
                else:
                    x.append(0)
                    y.append(0)
                m.append(mi)
                H.append(Hi)

            # Decrement to the next radial bin inward.
            ri = max(0.0, ri - dr)

        return x, y, m, H

#-------------------------------------------------------------------------------
# Specialized version to generate a variable y stepping to try
# and match a given density profile with (nearly) constant mass nodes.  This
# only supports the equivalent of the lattice or xstaggeredLattice methods.
#-------------------------------------------------------------------------------
class GenerateNodesMatchingYProfile2d(GenerateNodeDistribution2d):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, nx, ny, densityProfileMethod,
                 distributionType = "lattice",
                 xmin = None,
                 xmax = None,
                 nNodePerh = 2.01,
                 SPH = False):

        assert nx > 0
        assert ny > 0
        assert ((distributionType == "lattice" or
                 distributionType == "xstaggeredLattice") and
                xmin is not None and xmax is not None and
                xmin[0] < xmax[0] and xmin[1] < xmax[1])
        assert nNodePerh > 0.0

        self.nx = nx
        self.ny = ny
        self.distributionType = distributionType
        self.xmin = xmin
        self.xmax = xmax
        self.nNodePerh = nNodePerh

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(densityProfileMethod) == type(1.0):
            self.rho = ConstantRho(densityProfileMethod)
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.rho = densityProfileMethod
            self.densityProfileMethod = densityProfileMethod

        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 xmin, xmax)
        print("Total mass of %g in the range xmin = (%g, %g), xmax = (%g, %g)" % \
              (self.totalMass, xmin[0], xmin[1], xmax[0], xmax[1]))

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(nx*ny)
        assert self.m0 > 0.0
        print("Nominal mass per node of %g." % self.m0)

        # Generate the internal sets of positions, masses, and H
        if distributionType == "lattice":
            self.x, self.y, self.m, self.H = \
                    self.latticeDistribution(self.nx,
                                             self.ny,
                                             self.rho,
                                             self.xmin,
                                             self.xmax,
                                             self.nNodePerh)

        elif distributionType == "xstaggeredLattice":
            self.x, self.y, self.m, self.H = \
                    self.xstaggeredLatticeDistribution(self.nx,
                                                       self.ny,
                                                       self.rho,
                                                       self.xmin,
                                                       self.xmax,
                                                       self.nNodePerh)

        # If SPH has been specified, make sure the H tensors are round.
        if SPH:
            self.makeHround()

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)

        return

    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           xmin, xmax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0

        xwidth = xmax[0] - xmin[0]
        h = (xmax[1] - xmin[1])/nbins
        ymin = xmin[1]
        result = (xmin[1]*densityProfileMethod(Vector2d(*xmin)) +
                  xmax[1]*densityProfileMethod(Vector2d(*xmax)))
        for i in range(1, nbins):
            yi = ymin + i*h
            if i % 2 == 0:
                result += 4.0*yi*densityProfileMethod(Vector2d(0.0, yi))
            else:
                result += 2.0*yi*densityProfileMethod(Vector2d(0.0, yi))

        result *= xwidth*h/3.0
        return result

    #---------------------------------------------------------------------------
    # Seed positions/masses on a lattice
    #---------------------------------------------------------------------------
    def latticeDistribution(self, nx, ny, rho,
                            xmin = (0.0, 0.0),
                            xmax = (1.0, 1.0),
                            nNodePerh = 2.01):

        dx = (xmax[0] - xmin[0])/nx
        hx = 1.0/(nNodePerh*dx)

        x = []
        y = []
        m = []
        H = []

        # Start at ymax, and work our way down.
        yy = xmax[1]
        while yy > xmin[1]:
            rhoi = self.densityProfileMethod(Vector2d(0.0, yy))
            dy = sqrt(self.m0/rhoi)
            hy = 1.0/(nNodePerh*dy)
            H0 = SymTensor2d(hx, 0.0, 0.0, hy)
            for i in range(nx):
                xx = xmin[0] + (i + 0.5)*dx
                m0 = dx*dy*rho(Vector2d(xx, yy))
                x.append(xx)
                y.append(yy)
                m.append(m0)
                H.append(H0)
            # Decrement to the next y bin inward.
            yy = max(xmin[1], yy - dy)

        return x, y, m, H

    #---------------------------------------------------------------------------
    # Seed positions on a staggered (or offset) lattice
    #---------------------------------------------------------------------------
    def xstaggeredLatticeDistribution(self, nx, ny, rho,
                                      xmin = (0.0, 0.0),
                                      xmax = (1.0, 1.0),
                                      nNodePerh = 2.01):

        # Generate the basic distribution from the lattice option.
        x, y, m, H = self.latticeDistribution(nx, ny, rho,
                                              xmin, xmax,
                                              nNodePerh)

        # Now go through and stagger the positions in x.
        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        ddx = 0.25*dx

        assert len(x) == len(y)
        assert dy > 0.0
        for i in range(len(x)):
            x[i] += ddx * (-1.0)**int((y[i] - xmin[1])/dy + 0.1)

        return x, y, m, H


#-------------------------------------------------------------------------------
# Specialized 2-D NodeGenerator which initializes lattices in a slanted box.
#-------------------------------------------------------------------------------
class SlantedBoxNodeDistribution2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, 
                 n01,                 # Number aligned with x0->x1
                 n12,                 # Number aligned with x1->x2
                 rho,                 # Mass density
                 x0,                  # 4 points delineating the box 
                 x1,                  # in counter-clockwise order
                 x2,
                 x3,
                 nNodePerh = 2.01,
                 SPH = False):

        # Vectors defining the box.
        v01 = x1 - x0
        v12 = x2 - x1
        v23 = x3 - x2
        v30 = x0 - x3
        assert fuzzyEqual(v01.dot(v23), -(v01.magnitude()*v23.magnitude()), 1.0e-10)
        assert fuzzyEqual(v12.dot(v30), -(v12.magnitude()*v30.magnitude()), 1.0e-10)
        nhat01 = v01.unitVector()
        nhat12 = v12.unitVector()
        d01 = v01.magnitude()/n01
        d12 = v12.magnitude()/n12

        # Seed the points.
        self.x = []
        self.y = []
        for i01 in range(n01):
            for i12 in range(n12):
                r = x0 + (i01 + 0.5)*d01*nhat01 + (i12 + 0.5)*d12*nhat12
                self.x.append(r.x)
                self.y.append(r.y)

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        # Seed the masses.
        ntot = n01*n12
        voli = v01.cross(v12).z/ntot
        assert voli > 0.0
        mtot = sum([self.rho(Vector2d(xi, yi)) for xi, yi in zip(self.x, self.y)])*voli
        mi = mtot/ntot
        self.m = [mi for i in range(ntot)]

        # Seed the H tensors.
        Hi = SymTensor2d(1.0/(nNodePerh*d01), 
                         0.0, 
                         0.0,
                         1.0/(nNodePerh*d12*Vector2d(-nhat01.y, nhat01.x).dot(nhat12)))
        R = rotationMatrix2d(nhat01)
        Hi.rotationalTransform(R)
        if SPH:
            hi = sqrt(Hi.Determinant())
            Hi = SymTensor2d(hi, 0.0, 0.0, hi)
        self.H = [SymTensor2d(Hi) for i in range(ntot)]

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
        assert i >= 0 and i < len(self.x)
        return self.rho(self.localPosition(i))

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

#-------------------------------------------------------------------------------
# Factory function to convert any of the 2D generators to RZ.
#-------------------------------------------------------------------------------
def RZGenerator(generator):

    # Correct the mass.
    n = len(generator.m)
    for i in range(n):
        rhoi = generator.localMassDensity(i)
        assert rhoi > 0.0
        Ai = generator.localMass(i) / rhoi
        Ri = sqrt(Ai/pi)
        ri = generator.localPosition(i).y
        Vi = toroidalVolume(Ri, ri)
        generator.m[i] = generator.localPosition(i).y

    return generator
