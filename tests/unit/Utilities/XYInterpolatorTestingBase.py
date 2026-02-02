from Spheral import *
from SpheralTestUtilities import *
from testCubicHermiteInterpolator import err

import unittest
import numpy as np
from numpy.polynomial.polynomial import polyval2d as npP2D
import numpy.random

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Create a global random number generator.
import random
random.seed(4599281940)

#===============================================================================
# A generator for creating a range of x values to test
#===============================================================================
def xygen(n, xmin, xmax, ymin, ymax):
    yield xmin, ymin
    yield xmax, ymin
    yield xmax, ymax,
    yield xmin, ymax
    count = 0
    while count < n:
        yield random.uniform(xmin, xmax), random.uniform(ymin, ymax)
        count += 1

#===============================================================================
# A polynomial functor
#===============================================================================
class PolynomialFunctor(ScalarScalarScalarFunctor):
    def __init__(self, order, cmin, cmax):
        ScalarScalarScalarFunctor.__init__(self)
        self.order = order
        self.c = np.random.rand(order + 1, order + 1)
        self.c *= cmax - cmin
        self.c += cmin
        return

    def __call__(self, x, y):
        return npP2D(x, y, self.c)

# Take the gradient of above
class GradPolynomialFunctor(ScalarScalarSymTensor2dFunctor):
    def __init__(self, Fpoly):
        ScalarScalarSymTensor2dFunctor.__init__(self)
        self.Fpoly = Fpoly
        self.c = Fpoly.c
        self.order1 = Fpoly.order + 1
        return

    def __call__(self, x, y):
        result = SymTensor2d()
        for i in range(1, self.order1):
            for j in range(self.order1):
                result.xx = result.xx + self.c[i,j] * i * x**(i-1) * y**j
        for i in range(self.order1):
            for j in range(1, self.order1):
                result.yy = result.yy + self.c[i,j] * j * x**i * y**(j - 1)
        for i in range(1, self.order1):
            for j in range(1, self.order1):
                result.xy = result.xy + self.c[i,j] * i * j * x**(i - 1) * y**(j - 1)
        return result

#===============================================================================
# XYInterpolatorTestingBase
#===============================================================================
class XYInterpolatorTestingBase:

    #===========================================================================
    # Analytic coordinate of a grid point
    #===========================================================================
    def coord_ans(self, ix, nx, xmin, xmax, xlog):
        assert ix <= nx
        if xlog:
            return xmin + dx*exp(ix - nx)
        else:
            return xmin + ix*dx

    #===========================================================================
    # Analytic normalized coordinate
    #===========================================================================
    def eta_ans(self, x, ix, nx, xmin, xmax, xlog):
        x0 = self.coord_ans(ix, nx, xmin, xmax, xlog)
        x1 = self.coord_ans(ix + 1, nx, xmin, xmax, xlog)
        return max(0.0, min(1.0, (x - x0)/(x1 - x0)))

    #===========================================================================
    # Analytic lowerBound
    #===========================================================================
    def lowerBound_ans(self,
                       x, y,
                       xmin, xmax, ymin, ymax,
                       nx, ny, xlog, ylog):
        nx1 = nx - 1
        ny1 = ny - 1
        if xlog:
            Bx = (xmax - xmin)/(1.0 - exp(-nx1))
            Ax = xmax - Bx
            ix = int(nx1 + log((x - Ax)/Bx))
        else:
            ix = int((x - xmin)/(xmax - xmin)*nx1)
        if ylog:
            By = (ymax - ymin)/(1.0 - exp(-ny1))
            Ay = ymax - By
            iy = int(ny1 + log((y - Ay)/By))
        else:
            iy = int((y - ymin)/(ymax - ymin)*ny1)
        ix = min(nx1 - 1, max(0, ix))
        iy = min(ny1 - 1, max(0, iy))
        i0 = (self.order+1)*(self.order+1)*(nx1*iy + ix)
        return ix, iy, i0

    #===========================================================================
    # Plotting fun
    #===========================================================================
    def plotem(self, x, y, F, Finterp):
        x = np.arange(self.xmin, 1.001*self.xmax, (self.xmax - self.xmin)/100)
        y = np.arange(self.ymin, 1.001*self.ymax, (self.ymax - self.ymin)/100)
        x, y = np.meshgrid(x, y)
        NX, NY = x.shape
        z0 = np.array([[F(x[j][i], y[j][i]) for j in range(NY)] for i in range(NX)])
        z1 = np.array([[Finterp(x[j][i], y[j][i]) for j in range(NY)] for i in range(NX)])
        zdiff = z1 - z0

        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111, projection='3d')
        surf0 = ax0.plot_surface(x, y, z0, cmap=cm.coolwarm,
                                 linewidth=0, antialiased=False)
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, projection='3d')
        surf1 = ax1.plot_surface(x, y, z1, cmap=cm.coolwarm,
                                 linewidth=0, antialiased=False)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        surf2 = ax2.plot_surface(x, y, zdiff, cmap=cm.coolwarm,
                                 linewidth=0, antialiased=False)
        plt.show()

    #===========================================================================
    # lowerBound
    #===========================================================================
    def test_lowerBound(self):
        for xlog in (False, True):
            for ylog in (False, True):
                for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
                    F = PolynomialFunctor(2, -10.0, 10.0)
                    Finterp = self.generateInterpolator(nx, ny, xlog, ylog, F)
                    for itest in range(self.ntests):
                        x = random.uniform(self.xmin, self.xmax)
                        y = random.uniform(self.ymin, self.ymax)
                        ix0, iy0, i0 = self.lowerBound_ans(x, y,
                                                           self.xmin, self.xmax,
                                                           self.ymin, self.ymax,
                                                           nx, ny,
                                                           xlog, ylog)
                        ix, iy, i = Finterp.lowerBound(x, y)
                        passing = (ix == ix0) and (iy == iy0) and (i == i0)
                        self.assertTrue(passing,
                                        "lower bound failure: query=(ix=%i, iy=%i, i0=%i), ans=((ix=%i, iy=%i, i0=%i) @ (%g, %g), (nx,ny)=(%i,%i), (xlog,ylog)=(%s,%s)" % (ix, iy, i, ix0, iy0, i0, x, y, nx, ny, xlog, ylog))

    #===========================================================================
    # Interpolate a linear function
    #===========================================================================
    def test_interp_linear(self):
        for xlog in (False, True):
            for ylog in (False, True):
                for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
                    for itest in range(self.ntests):
                        F = PolynomialFunctor(1, -10.0, 10.0)
                        Finterp = self.generateInterpolator(nx, ny, xlog, ylog, F)
                        tol = self.tol[1] / sqrt(nx*ny)
                        for x, y in xygen(self.n, self.xmin, self.xmax, self.ymin, self.ymax):
                            passing = err(Finterp(x, y), F(x, y)) < tol
                            # if not passing:
                            #     self.plotem(x, y, F, Finterp)
                            self.assertTrue(passing,
                                            "Interpolation off @ (%g,%g) (xlog,ylog)=(%s,%s) (nx,ny)=(%i,%i): %g != %g, err=%g" % (x, y, xlog, ylog, nx, ny, Finterp(x, y), F(x, y), err(Finterp(x,y), F(x,y))))

    #===========================================================================
    # Interpolate a quadratic function 
    #===========================================================================
    def test_interp_quadratic(self):
        for xlog in (False, True):
            for ylog in (False, True):
                for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
                    for itest in range(self.ntests):
                        F = PolynomialFunctor(2, -10.0, 10.0)
                        Finterp = self.generateInterpolator(nx, ny, xlog, ylog, F)
                        tol = self.tol[2] / sqrt(nx*ny)
                        for x, y in xygen(self.n, self.xmin, self.xmax, self.ymin, self.ymax):
                            passing = err(Finterp(x, y), F(x, y)) < tol
                            # if not passing:
                            #     self.plotem(x, y, F, Finterp)
                            self.assertTrue(passing,
                                            "Interpolation off @ (%g,%g) (xlog,ylog)=(%s,%s) (nx,ny)=(%i,%i): %g != %g, err=%g" % (x, y, xlog, ylog, nx, ny, Finterp(x, y), F(x, y), err(Finterp(x,y), F(x,y))))

    #===========================================================================
    # Interpolate a cubic function
    #===========================================================================
    def test_interp_cubic(self):
        for xlog in (False, True):
            for ylog in (False, True):
                for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
                    for itest in range(self.ntests):
                        F = PolynomialFunctor(3, -10.0, 10.0)
                        Finterp = self.generateInterpolator(nx, ny, xlog, ylog, F)
                        tol = self.tol[3] / sqrt(nx*ny)
                        # sys.stderr.write("(%i, %i): %s\n" % (nx, ny, Finterp))

                        for x, y in xygen(self.n, self.xmin, self.xmax, self.ymin, self.ymax):
                            passing = err(Finterp(x, y), F(x, y)) < tol
                            # if not passing:
                            #     self.plotem(x, y, F, Finterp)
                            self.assertTrue(passing,
                                            "Interpolation off @ (%g,%g) (xlog,ylog)=(%s,%s) (nx,ny)=(%i,%i): %g != %g, err=%g" % (x, y, xlog, ylog, nx, ny, Finterp(x, y), F(x, y), err(Finterp(x,y), F(x,y))))
