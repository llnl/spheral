#ATS:test(SELF, label="Test numerical integration using Simpson's rule.")

from Spheral import *
from SpheralTestUtilities import *

import unittest

# Build a random number generator.
import random
random.seed(4599281940)

#===============================================================================
# Implement a simple linear function in x.
# y = a + b*x
#===============================================================================
class LinearFunction(ScalarScalarFunctor):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        ScalarScalarFunctor.__init__(self)
        return
    def __call__(self, x):
        return self.a + self.b*x

#===============================================================================
# Implement a quadratic function in x.
# y = a + b*x * c*x^2
#===============================================================================
class QuadraticFunction(ScalarScalarFunctor):
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        ScalarScalarFunctor.__init__(self)
        return
    def __call__(self, x):
        return self.a + (self.b + self.c*x)*x

#===============================================================================
# Test the simpsonsIntegration numerical integration function.
#===============================================================================
class TestSimpsonsRuleIntegration(unittest.TestCase):

    #===========================================================================
    # Set up.
    #===========================================================================
    def setUp(self):
        self.ntests = 100
        self.xmin, self.xmax = -1e100, 1e100
        self.ymin, self.ymax = -100.0, 100.0
        self.numBins = 1000
        return

    #===========================================================================
    # Iterate over a bunch of randomly selected linear functions.
    #===========================================================================
    def testRandomLinearFunctions(self):
        for i in range(self.ntests):
            a = random.uniform(self.ymin, self.ymax)
            b = random.uniform(-100.0, 100.0)
            x0 = random.uniform(self.xmin, 0.5*self.xmax)
            x1 = random.uniform(x0, self.xmax)
            result = simpsonsIntegrationDouble(LinearFunction(a, b),
                                               x0, x1, self.numBins)
            ans = (a + 0.5*b*x1)*x1 - (a + 0.5*b*x0)*x0
            self.assertTrue(fuzzyEqual(result, ans, 1.0e-10),
                            "Failed linear integration:  %g != %g" % (result, ans))
        return

    #===========================================================================
    # Iterate over a bunch of randomly selected quadratic functions.
    #===========================================================================
    def testRandomQuadraticFunctions(self):
        for i in range(self.ntests):
            a = random.uniform(self.ymin, self.ymax)
            b = random.uniform(-100.0, 100.0)
            c = random.uniform(-100.0, 100.0)
            x0 = random.uniform(self.xmin, 0.5*self.xmax)
            x1 = random.uniform(x0, self.xmax)
            result = simpsonsIntegrationDouble(QuadraticFunction(a, b, c),
                                               x0, x1, self.numBins)
            ans = (a*x1 + 0.5*b*x1**2 + c/3.0*x1**3) - (a*x0 + 0.5*b*x0**2 + c/3.0*x0**3)
            self.assertTrue(fuzzyEqual(result, ans, 1.0e-10),
                            "Failed quadratic integration:  %g != %g" % (result, ans))
        return

if __name__ == "__main__":
    unittest.main()

