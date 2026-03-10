#ATS:test(SELF, label="CubicHermiteInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest
import matplotlib.pyplot as plt
import numpy as np

# Create a global random number generator.
# We force a fixed seed to cut down random failures in CI testing.
import random
random.seed(49884928910350901743)

#===========================================================================
# Measure the relative difference between two numbers
#===========================================================================
def err(a, b):
    return abs(a - b)/(1e-20 + abs(a) + abs(b))

#===========================================================================
# A generator for creating a range of x values to test
#===========================================================================
def xgen(n, xmin, xmax):
    count = 0
    if n < 3:
        yield random.uniform(xmin, xmax)
        count += 1
    else:
        while count < n:
            if count == 0:
                yield xmin
            elif count == n - 1:
                yield xmax
            else:
                yield random.uniform(xmin, xmax)
            count += 1

#===============================================================================
# Functors to generate the problem
#===============================================================================
class Fquad(ScalarScalarFunctor):
    def __init__(self, A, B, C):
        ScalarScalarFunctor.__init__(self)
        self.A = A
        self.B = B
        self.C = C
        return
    def __call__(self, x):
        return self.A + self.B*x + self.C*x*x
    def prime(self, x):
        return self.B + 2.0*self.C*x
    def prime2(self, x):
        return 2.0*self.C

class Fcubic(ScalarScalarFunctor):
    def __init__(self, A, B, C, D):
        ScalarScalarFunctor.__init__(self)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        return
    def __call__(self, x):
        return self.A + self.B*x + self.C*x*x + self.D*x*x*x
    def prime(self, x):
        return self.B + 2.0*self.C*x + 3.0*self.D*x*x
    def prime2(self, x):
        return 2.0*self.C + 6.0*self.D*x

class Fquintic(ScalarScalarFunctor):
    def __init__(self, A, B, C, D, E, F):
        ScalarScalarFunctor.__init__(self)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.F = F
        return
    def __call__(self, x):
        return self.A + self.B*x + self.C*x*x + self.D*x**3 + self.E*x**4 + self.F*x**5
    def prime(self, x):
        return self.B + 2.0*self.C*x + 3.0*self.D*x*x + 4.0*self.E*x**3 + 5.0*self.F*x**4
    def prime2(self, x):
        return 2.0*self.C + 6.0*self.D*x + 12.0*self.E*x*x + 20.0*self.F*x**3

class Fgrad(ScalarScalarFunctor):
    def __init__(self, F):
        ScalarScalarFunctor.__init__(self)
        self.F = F
        return
    def __call__(self, x):
        return self.F.prime(x)

#===============================================================================
# TestCubicHermiteInterpolator
#===============================================================================
class TestCubicHermiteInterpolator(unittest.TestCase):

    #===========================================================================
    # Set up
    #===========================================================================
    def setUp(self):
        self.nfunc = 1000 # 50000
        self.nsamples = 1000
        self.n = 100
        return

    #===========================================================================
    # Plotting fun
    #===========================================================================
    def plotem(self, x, xmin, xmax, F, Finterp):
        nx = 100
        x = np.arange(xmin, 1.001*xmax, (xmax - xmin)/nx)
        y0 = np.array([F(xi) for xi in x])
        y1 = np.array([Finterp(xi) for xi in x])
        ydiff = y1 - y0
        dydx0 = np.array([F.prime(xi) for xi in x])
        dydx1 = np.array([Finterp.prime(xi) for xi in x])
        dydiff = dydx1 - dydx0
        d2ydx0 = np.array([F.prime2(xi) for xi in x])
        d2ydx1 = np.array([Finterp.prime2(xi) for xi in x])
        d2ydiff = d2ydx1 - d2ydx0

        def _plotit(x, y, label):
            fig = plt.figure()
            plt.plot(x, y, label=label)
            plt.legend()
            return fig

        fig1 = _plotit(x, y0, label="F(x) Analytic")
        fig2 = _plotit(x, y1, label="<F(x)> Interpolation")
        fig3 = _plotit(x, ydiff, label="F(x) Error")
        fig4 = _plotit(x, dydx0, label="dF/dx(x) Analytic")
        fig5 = _plotit(x, dydx1, label="<dF/dx(x)> Interpolation")
        fig6 = _plotit(x, dydiff, label="dF/dx(x) Error")
        fig7 = _plotit(x, d2ydx0, label="d$^2$F/dx$^2$(x) Analytic")
        fig8 = _plotit(x, d2ydx1, label="<d$^2$F/dx$^2$(x)> Interpolation")
        fig9 = _plotit(x, d2ydiff, label="d$^2$F/dx$^2$(x) Error")

        plt.show()

    #===========================================================================
    # Boilerplate for checking solution
    #===========================================================================
    def checkError(self, xmin, xmax,
                   func,        # analytic function
                   F,           # interpolator
                   Ftol,        # func interpolation tolerance
                   F1tol,       # first derivative tolerance
                   F2tol,       # second derivative tolerance
                   errorLabel,
                   tolFunc = None,
                   checkMonotonicity = False):

        def tolFunc(params, x):
            tol0, dx, hx, A = params
            eta = min(x - xmin, xmax - x)
            if eta < dx:
                return 1e6
            else:
                return tol0*(1.0 + A*exp(-(eta/hx)**2))

        for x in xgen(self.nsamples, xmin, xmax):
            passing = err(F(x), func(x)) < tolFunc(Ftol, x)
            # if not passing:
            #     self.plotem(x, xmin, xmax, func, F)
            self.assertTrue(passing,
                            "Error interpolating F({}) for {}: {} != {}, err = {}".format(x, errorLabel, F(x), func(x), err(F(x), func(x))))

            # Check the first derivative
            passing = err(F.prime(x), func.prime(x)) < tolFunc(F1tol, x)
            # if not passing:
            #     self.plotem(x, xmin, xmax, func, F)
            self.assertTrue(passing,
                            "Error interpolating dF/dx({}) for {}: {} != {}, err = {}".format(x, errorLabel, F.prime(x), func.prime(x), err(F.prime(x), func.prime(x))))

            # Check the second derivative
            passing = err(F.prime2(x), func.prime2(x)) < tolFunc(F2tol, x)
            # if not passing:
            #     self.plotem(x, xmin, xmax, func, F)
            self.assertTrue(passing,
                            "Error interpolating d^2F/dx^2({}) for {}: {} != {}, err = {}".format(x, errorLabel, F.prime2(x), func.prime2(x), err(F.prime2(x), func.prime2(x))))

            # If requested, check for monotonicity in interpolation
            if checkMonotonicity:
                i0 = F.lowerBound(x)
                passing = (F(x) - F[i0])*(F(x) - F[i0 + 1]) <= 0.0
                if not passing:
                    #print(F.vals)
                    self.plotem(x, xmin, xmax, func, F)
                self.assertTrue(passing,
                                "Failing monotonicity test for {}: F({}) = {} not in [{}, {}]".format(errorLabel, x, F(x), F[i0], F[i0 + 1]))

        return

    #===========================================================================
    # Interpolate a quadratic function (without gradient info)
    #===========================================================================
    def test_quad_interp(self):
        xmin = -10.0
        xmax =  40.0
        dx = (xmax - xmin)/self.n
        hx = 0.05*(xmax - xmin)
        for ifunc in range(self.nfunc):
            A = random.uniform(-100.0, 100.0)
            B = random.uniform(-100.0, 100.0)
            C = random.uniform(-100.0, 100.0)
            func = Fquad(A, B, C)
            F = CubicHermiteInterpolator(xmin, xmax, 10*self.n, func)  # Without the analytic gradient we benefit from more fiting points
            self.checkError(xmin, xmax, func, F, 
                            Ftol  = (1.0e-4, 0.0, hx, 100.0),
                            F1tol = (1.0e-3, dx, hx, 500.0),
                            F2tol = (5.0e-3, dx, hx, 1000.0),
                            errorLabel = "quadratic function")

    #===========================================================================
    # Interpolate a quadratic function (with gradient info)
    #===========================================================================
    def test_quad_interp_with_grad(self):
        xmin = -10.0
        xmax =  40.0
        for ifunc in range(self.nfunc):
            A = random.uniform(-100.0, 100.0)
            B = random.uniform(-100.0, 100.0)
            C = random.uniform(-100.0, 100.0)
            func = Fquad(A, B, C)
            F = CubicHermiteInterpolator(xmin, xmax, self.n, func, Fgrad(func))
            self.checkError(xmin, xmax, func, F,
                            Ftol  = (5e-9, 0.0, 1.0, 0.0),
                            F1tol = (5e-8, 0.0, 1.0, 0.0),
                            F2tol = (1e-6, 0.0, 1.0, 0.0),
                            errorLabel = "quadratic function")

    #===========================================================================
    # Interpolate a quadratic function enforcing monotonicity
    #===========================================================================
    def test_quad_interp_monotonic(self):
        xmin = -10.0
        xmax =  40.0
        for ifunc in range(self.nfunc):
            A = random.uniform(-100.0, 100.0)
            B = random.uniform(-100.0, 100.0)
            C = random.uniform(-100.0, 100.0)
            func = Fquad(A, B, C)
            F = CubicHermiteInterpolator(xmin, xmax, self.n, func)
            F.makeMonotonic()
            self.checkError(xmin, xmax, func, F, 
                            Ftol  = (2.0, 0.0, 1.0, 0.0),           # Tolerance has to be way looser when using monotonicity
                            F1tol = (2.0, 0.0, 1.0, 0.0),
                            F2tol = (2.0, 0.0, 1.0, 0.0),
                            errorLabel = "quadratic function with monotonicity",
                            checkMonotonicity = True)

    #===========================================================================
    # Interpolate a cubic function (func only)
    #===========================================================================
    def test_cubic_interp(self):
        xmin = -10.0
        xmax =  40.0
        dx = (xmax - xmin)/self.n
        hx = 0.05*(xmax - xmin)
        for ifunc in range(self.nfunc):
            A = random.uniform(-100.0, 100.0)
            B = random.uniform(-100.0, 100.0)
            C = random.uniform(-100.0, 100.0)
            D = random.uniform(-100.0, 100.0)
            func = Fcubic(A, B, C, D)
            F = CubicHermiteInterpolator(xmin, xmax, 10*self.n, func)  # Without the analytic gradient we benefit from more fiting points
            self.checkError(xmin, xmax, func, F,
                            Ftol  = (0.5, 0.0, hx, 100.0),
                            F1tol = (1.0, dx, hx, 500.0),
                            F2tol = (2.0, dx, hx, 1000.0),
                            errorLabel = "cubic function")

    #===========================================================================
    # Interpolate a cubic function (func + grad)
    #===========================================================================
    def test_cubic_interp_with_grad(self):
        xmin = -10.0
        xmax =  40.0
        for ifunc in range(self.nfunc):
            A = random.uniform(-100.0, 100.0)
            B = random.uniform(-100.0, 100.0)
            C = random.uniform(-100.0, 100.0)
            D = random.uniform(-100.0, 100.0)
            func = Fcubic(A, B, C, D)
            F = CubicHermiteInterpolator(xmin, xmax, self.n, func, Fgrad(func))
            self.checkError(xmin, xmax, func, F, 
                            Ftol  = (5e-9, 0.0, 1.0, 0.0),
                            F1tol = (5e-8, 0.0, 1.0, 0.0),
                            F2tol = (1e-6, 0.0, 1.0, 0.0),
                            errorLabel = "cubic function")

    #===========================================================================
    # Interpolate a cubic function enforcing monotonicity
    #===========================================================================
    def test_cubic_interp_monotonic(self):
        xmin = -10.0
        xmax =  40.0
        for ifunc in range(self.nfunc):
            A = random.uniform(-100.0, 100.0)
            B = random.uniform(-100.0, 100.0)
            C = random.uniform(-100.0, 100.0)
            D = random.uniform(-100.0, 100.0)
            func = Fcubic(A, B, C, D)
            F = CubicHermiteInterpolator(xmin, xmax, self.n, func)
            F.makeMonotonic()
            self.checkError(xmin, xmax, func, F, 
                            Ftol  = (2.0, 0.0, 1.0, 0.0),
                            F1tol = (2.0, 0.0, 1.0, 0.0),
                            F2tol = (2.0, 0.0, 1.0, 0.0),
                            errorLabel = "cubic function with monotonicity",
                            checkMonotonicity = True)

if __name__ == "__main__":
    unittest.main()
