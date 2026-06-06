#ATS:test(SELF, label="QuadraticInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
random.seed(4599281940)

#===========================================================================
# A generator for creating a range of x values to test
#===========================================================================
def xgen(n, xmin, xmax):
    assert n > 2
    count = 0
    while count < n:
        if count == 0:
            yield xmin
        elif count == n - 1:
            yield xmax
        else:
            yield random.uniform(xmin, xmax)
        count += 1

#===============================================================================
# Functor to generate the problem
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

#===============================================================================
# TestQuadraticInterpolator
#===============================================================================
class TestQuadraticInterpolator(unittest.TestCase):

    #===========================================================================
    # Set up
    #===========================================================================
    def setUp(self):
        self.ntests = 100000
        self.n = 100
        self.xmin = -10.0
        self.xmax =  40.0
        self.dx = (self.xmax - self.xmin)/(self.n - 1)
        self.A = random.uniform(-100.0, 100.0)
        self.B = random.uniform(-100.0, 100.0)
        self.C = random.uniform(-100.0, 100.0)
        self.func = Fquad(self.A, self.B, self.C)
        self.F = QuadraticInterpolator(self.xmin, self.xmax, self.n, self.func)
        self.fuzz = 1.0e-10
        return

    #===========================================================================
    # Interpolate y
    #===========================================================================
    def test_yinterp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            self.assertTrue(fuzzyEqual(self.F(x), self.func(x), self.fuzz),
                            "Error interpolating F(x): %g != %g" % (self.F(x), self.func(x)))

    #===========================================================================
    # Interpolate dy/dx
    #===========================================================================
    def test_dyinterp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            self.assertTrue(fuzzyEqual(self.F.prime(x), self.func.prime(x), self.fuzz),
                            "Error interpolating F'(x): %g != %g" % (self.F.prime(x), self.func.prime(x)))

    #===========================================================================
    # Interpolate d^2y/dx^2
    #===========================================================================
    def test_dy2interp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            self.assertTrue(fuzzyEqual(self.F.prime2(x), self.func.prime2(x), self.fuzz),
                            "Error interpolating F''(x): %g != %g" % (self.F.prime2(x), self.func.prime2(x)))

if __name__ == "__main__":
    unittest.main()
