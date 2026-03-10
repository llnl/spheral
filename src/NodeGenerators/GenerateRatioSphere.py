from math import *
import mpi

from NodeGeneratorBase import *

from Spheral import Vector2d, Tensor2d, SymTensor2d, CylindricalBoundary, rotationMatrix2d, Polygon
from Spheral import Vector3d, Tensor3d, SymTensor3d, CylindricalBoundary, rotationMatrix3d, Polyhedron
from Spheral import CylindricalBoundary, generateCylDistributionFromRZ
from Spheral import vector_of_int, vector_of_double, vector_of_vector_of_double, vector_of_SymTensor3d
from Spheral import polySecondMoment2d, polySecondMoment3d

#-------------------------------------------------------------------------------
# This version ratios from the center out in 2D.  Kind of a misnomer with the
# sphere thing...
#-------------------------------------------------------------------------------
class GenerateRatioSphere2d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 drStart, drRatio,
                 rho,
                 rmin,
                 rmax,
                 startFromCenter = True,
                 thetamin = 0.0,
                 thetamax = 0.5*pi,
                 ntheta = 1,
                 center = (0.0, 0.0),
                 distributionType = "constantDTheta",   # one of (constantDTheta, constantNTheta)
                 aspectRatio = 1.0,                     # only for constantDTheta
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None,
                 perturbFunc = None,
                 skipFinalInitialization = False):     # Only for use by 3D ratio generator (not users)

        nNodePerh = float(nNodePerh)  # Just to be sure...

        assert drStart > 0.0
        assert drRatio > 0.0
        assert nNodePerh > 0.0
        assert rmin >= 0.0
        assert rmax > rmin
        assert thetamax > thetamin
        assert distributionType.lower() in ("constantdtheta", "constantntheta")
        
        self.center = center

        # Did we get passed a function or a constant for the density?
        if type(rho) == type(1.0):
            def rhofunc(posi):
                return rho
        else:
            rhofunc = rho
        self.rhofunc = rhofunc

        # Do we have a perturbation function?
        if not perturbFunc:
            perturbFunc = lambda x: x

        self.x, self.y, self.m, self.H = [], [], [], []

        constantN = (distributionType.lower() == "constantntheta")
        Dtheta = thetamax - thetamin

        nthetamin = max(2, int(Dtheta/(0.5*pi) + 0.5)*2)

        # Decide the actual drStart we're going to use to arrive at an integer number of radial bins.
        if abs(drRatio - 1.0) > 1e-4:
            neff = max(1, int(log(1.0 - (rmax - rmin)*(1.0 - drRatio)/drStart)/log(drRatio) + 0.5))
            drStart = (rmax - rmin)*(1.0 - drRatio)/(1.0 - drRatio**neff)
        else:
            neff = max(1, int((rmax - rmin)/drStart + 0.5))
            drStart = (rmax - rmin)/neff
        print("Adjusting initial radial spacing to %g in order to create an integer radial number of bins %i." % (drStart, neff))

        # Step in radius (in or out) until we span the full radial range.
        dr = drStart
        for i in range(neff):
            if abs(drRatio - 1.0) > 1e-4:
                if startFromCenter:
                    r0 = min(rmax, rmin + drStart*(1.0 - drRatio**i)/(1.0 - drRatio))
                    r1 = min(rmax, rmin + drStart*(1.0 - drRatio**(i + 1))/(1.0 - drRatio))
                    r0hr = rmin + drStart*(1.0 - drRatio**max(0, i - nNodePerh))/(1.0 - drRatio)
                    r1hr = rmin + drStart*(1.0 - drRatio**(      i + nNodePerh))/(1.0 - drRatio)
                else:
                    r0 = max(rmin, rmax - drStart*(1.0 - drRatio**(i + 1))/(1.0 - drRatio))
                    r1 = max(rmin, rmax - drStart*(1.0 - drRatio**i)/(1.0 - drRatio))
                    r0hr = rmax - drStart*(1.0 - drRatio**(      i + nNodePerh))/(1.0 - drRatio)
                    r1hr = rmax - drStart*(1.0 - drRatio**max(0, i - nNodePerh))/(1.0 - drRatio)
            else:
                r0 = min(rmax, rmin + i*drStart)
                r1 = min(rmax, rmin + (i + 1)*drStart)
                r0hr = rmin + (i - nNodePerh)*drStart
                r1hr = rmin + (i + nNodePerh)*drStart

            dr = r1 - r0
            ri = 0.5*(r0 + r1)
            li = Dtheta*ri
            if constantN:
                ntheta = ntheta
            else:
                ntheta = max(nthetamin, int(li/dr*aspectRatio))
            dtheta = Dtheta/ntheta

            # Find the radial and azimuthal smoothing lengths we should use.  We have to be
            # careful for extrememely high aspect ratios that the points will overlap the expected
            # number of neighbors taking into account the curvature of the local point distribution.
            # This means hr might need to be larger than we would naively expect...
            #hdelta = 2.0*ri*(sin(0.5*nNodePerh*dtheta))**2
            r0hr -= 2.0*r1hr*(sin(0.5*nNodePerh*dtheta))**2
            r1hr += 2.0*r1hr*(sin(0.5*nNodePerh*dtheta))**2
            hr = max(r1hr - ri, ri - r0hr)
            ha = nNodePerh * ri*dtheta
            # box = Polygon([Vector2d(r0hr, -ha), Vector2d(r1hr, -ha),
            #                Vector2d(r1hr,  ha), Vector2d(r0hr,  ha)])
            # Hi = polySecondMoment2d(box, box.centroid).sqrt().Inverse()
            
            for j in range(ntheta):
                theta0 = thetamin + j*dtheta
                theta1 = thetamin + (j + 1)*dtheta
                pos0 = perturbFunc(Vector2d(r0*cos(theta0), r0*sin(theta0)))
                pos1 = perturbFunc(Vector2d(r1*cos(theta0), r1*sin(theta0)))
                pos2 = perturbFunc(Vector2d(r1*cos(theta1), r1*sin(theta1)))
                pos3 = perturbFunc(Vector2d(r0*cos(theta1), r0*sin(theta1)))
                areai = 0.5*((pos1 - pos0).cross(pos2 - pos0).z +
                             (pos2 - pos0).cross(pos3 - pos0).z)
                posi = 0.5*(r0 + r1)*Vector2d(cos(0.5*(theta0 + theta1)),
                                              sin(0.5*(theta0 + theta1)))
                mi = areai*self.rhofunc(posi)
                self.x.append(posi.x + center[0])
                self.y.append(posi.y + center[1])
                self.m.append(mi)
                if SPH:
                    hi = sqrt(hr*ha)
                    self.H.append(SymTensor2d(1.0/hi, 0.0, 0.0, 1.0/hi))
                else:
                    self.H.append(SymTensor2d(1.0/hr, 0.0, 0.0, 1.0/ha))
                    runit = posi.unitVector()
                    T = rotationMatrix2d(runit).Transpose()
                    self.H[-1].rotationalTransform(T)

        # # Do a numerical integral to get the expected total mass.
        # class integfunc(ScalarFunctor):
        #     def __init__(self, rho, Dtheta):
        #         ScalarFunctor.__init__(self)
        #         self.rho = rho
        #         self.Dtheta = Dtheta
        #         return
        #     def __call__(self, ri):
        #         return Dtheta*ri*self.rho(ri)
        # M1 = simpsonsIntegrationDouble(integfunc(rhofunc, Dtheta), rmin, rmax, 10000)

        # # Make sure the total mass is what we intend it to be, by applying
        # # a multiplier to the particle masses.
        # M0 = sum(self.m)
        # assert M0 > 0.0
        # massCorrection = M1/M0
        # for i in xrange(len(self.m)):
        #     self.m[i] *= massCorrection
        # print "Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, M1)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.m, self.H = rejecter(self.x,
                                                      self.y,
                                                      self.m,
                                                      self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        if not skipFinalInitialization:
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
        ri = sqrt((self.x[i] - self.center[0])**2 + (self.y[i] - self.center[1])**2)
        return self.rhofunc(ri)
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
    
#--------------------------------------------------------------------------------
# 3D version, actual sphere.
# Based on spinning the case above.
#--------------------------------------------------------------------------------
class GenerateRatioSphere3d(NodeGeneratorBase):

    def __init__(self, 
                 drCenter, drRatio,
                 rho,
                 rmin,
                 rmax,
                 startFromCenter = True,
                 thetamin = 0.0,
                 thetamax = 0.5*pi,
                 phi = pi,
                 ntheta = 1,
                 center = (0.0, 0.0, 0.0),
                 distributionType = "constantDTheta",   # one of (constantDTheta, constantNTheta)
                 aspectRatio = 1.0,                     # only for constantDTheta
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None):

        assert thetamax <= pi

        self.gen2d = GenerateRatioSphere2d(drStart = drCenter, 
                                           drRatio = drRatio, 
                                           rho = rho, 
                                           rmin = rmin, 
                                           rmax = rmax, 
                                           startFromCenter = startFromCenter, 
                                           thetamin = thetamin, 
                                           thetamax = thetamax, 
                                           ntheta = ntheta, 
                                           center = (0.0, 0.0), 
                                           distributionType = distributionType, 
                                           aspectRatio = aspectRatio,
                                           nNodePerh = nNodePerh, 
                                           SPH = SPH,
                                           skipFinalInitialization = True)

        # At this point every process generated the full 2D distribution redundantly.
        # Add a z coordinate array of the same size.
        self.x = self.gen2d.x
        self.y = self.gen2d.y
        self.m = self.gen2d.m
        self.H = self.gen2d.H
        n = len(self.x)
        self.z = [0.0]*n
        self.globalIDs = [0]*n

        # Convert the 2-D H tensors to 3-D, and correct the masses.
        for i in range(n):
            xi = self.x[i]
            yi = self.y[i]
            H2d = SymTensor2d(self.H[i])
            H2dinv = H2d.Inverse()

            hxy0 = 0.5*(H2dinv.Trace())
            dphi = CylindricalBoundary.angularSpacing(yi, hxy0, nNodePerh, 2.0)
            assert dphi > 0.0
            nsegment = max(1, int(phi/dphi + 0.5))
            dphi = phi/nsegment

            hz = dphi*yi*nNodePerh
            self.H[i] = SymTensor3d(H2d.xx, H2d.xy, 0.0,
                                    H2d.yx, H2d.yy, 0.0,
                                    0.0,    0.0,    1.0/hz)
            if SPH:
                h0 = self.H[i].Determinant()**(1.0/3.0)
                self.H[-1] = SymTensor3d(h0, 0.0, 0.0,
                                         0.0, h0, 0.0,
                                         0.0, 0.0, h0)

            # Convert the mass to the full hoop mass, which will then be used in
            # generateCylDistributionFromRZ to compute the actual nodal masses.
            mi = self.m[i]
            circ = 2.0*pi*yi
            mhoop = mi*circ
            self.m[i] = mhoop

        assert len(self.m) == n
        assert len(self.H) == n

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
        extras = []
        xvec = vector_of_double(self.x)
        yvec = vector_of_double(self.y)
        zvec = vector_of_double(self.z)
        mvec = vector_of_double(self.m)
        Hvec = vector_of_SymTensor3d(self.H)
        globalIDsvec = vector_of_int(self.globalIDs)
        extrasVec = vector_of_vector_of_double()
        for extra in extras:
            extrasVec.append(vector_of_double(extra))
        generateCylDistributionFromRZ(xvec, yvec, zvec, mvec, Hvec, globalIDsvec,
                                      extrasVec,
                                      nNodePerh, kernelExtent, phi,
                                      mpi.rank, mpi.procs)
        self.x = [x + center[0] for x in xvec]
        self.y = [x + center[1] for x in yvec]
        self.z = [z + center[2] for z in zvec]
        self.m = list(mvec)
        self.H = [SymTensor3d(x) for x in Hvec]
        self.globalIDs = list(globalIDsvec)
        for i in range(len(extras)):
            extras[i] = list(extrasVec[i])

        self.center = Vector2d(*center)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.z,
                                                              self.m,
                                                              self.H)
        # Initialize the base class.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        return self.gen2d.rhofunc((Vector2d(self.x[i], self.y[i]) - self.center).magnitude())

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]

