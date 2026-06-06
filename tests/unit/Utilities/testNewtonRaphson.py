from Spheral import *
from SpheralTestUtilities import *

import unittest
import Gnuplot
import CXXTests

# Build a random number generator.
import random
random.seed(4599281940)

#===============================================================================
# Test the newtonRaphson root finding function.
#===============================================================================
class TestNewtonRaphson(unittest.TestCase):

    #===========================================================================
    # Set up, create arrays of the function values.
    #===========================================================================
    def setUp(self):
        self.ntests = 1000
        self.xmin = -1e10
        self.xmax = 1e10
        self.xaccuracy = 1e-15
        self.yaccuracy = 1e-15
        self.maxIterations = 100
        return

    #===========================================================================
    # Iterate over a bunch of randomly selected cubic functions, and find the
    # known roots.
    #===========================================================================
    def testRoots(self):
        for i in range(self.ntests):

            # Randomly pick three roots.  We want to know them
            # in sorted order too.
            xlist = [random.uniform(self.xmin, self.xmax),
                     random.uniform(self.xmin, self.xmax),
                     random.uniform(self.xmin, self.xmax)]
            xlist.sort()
            x0 = xlist[0]
            x1 = xlist[1]
            x2 = xlist[2]

            # Build our test functor.
            func = CXXTests.TestFunction(x0, x1, x2)

            # Can we find each root individually?
            # x0
            x0test = CXXTests.testNewtonRaphsonRoot(func,
                                                    self.xmin,
                                                    0.5*(x0 + x1))
            self.assertTrue(fuzzyEqual(x0test, x0, self.xaccuracy),
                            "Failed to find root %g != %g of (%g, %g, %g)" % (x0test, x0, x0, x1, x2))

            # x1
            x1test = CXXTests.testNewtonRaphsonRoot(func,
                                                    0.5*(x0 + x1),
                                                    0.5*(x1 + x2))
            self.assertTrue(fuzzyEqual(x1test, x1, self.xaccuracy),
                            "Failed to find root %g != %g of (%g, %g, %g)" % (x1test, x1, x0, x1, x2))

            # x2
            x2test = CXXTests.testNewtonRaphsonRoot(func,
                                                    0.5*(x1 + x2),
                                                    self.xmax)
            self.assertTrue(fuzzyEqual(x2test, x2, self.xaccuracy),
                            "Failed to find root %g != %g of (%g, %g, %g)" % (x2test, x2, x0, x1, x2))

            # Do we find one of the roots if we span all three?
            xalltest = CXXTests.testNewtonRaphsonRoot(func,
                                                      self.xmin,
                                                      self.xmax)
            self.assertTrue(fuzzyEqual(xalltest, x0, self.xaccuracy) or
                            fuzzyEqual(xalltest, x1, self.xaccuracy) or
                            fuzzyEqual(xalltest, x2, self.xaccuracy),
                            "Failed to find root in full range: %g != (%g, %g, %g)" % (xalltest, x0, x1, x2))

        return

if __name__ == "__main__":
    unittest.main()

