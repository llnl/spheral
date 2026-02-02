#ATS:test(SELF, label="3-D line segment & plane distance unit tests")
from Spheral3d import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
random.seed(4599281940)

#===============================================================================
# Test our methods for computing distances in 3-D.
#===============================================================================
class TestDistances3d(unittest.TestCase):

    #===========================================================================
    # 
    #===========================================================================
    def setUp(self):
        self.ntests = 100
        self.multMin = 0.001
        self.multMax = 1e6
        return

    #===========================================================================
    # Randomly distort two line segments.
    #===========================================================================
    def randomDistortion(self, a0, a1, b0, b1):
        l = random.uniform(self.multMin, self.multMax)
        T = l*rotationMatrix(Vector(random.uniform(0.0, 1.0),
                                    random.uniform(0.0, 1.0),
                                    random.uniform(0.0, 1.0)).unitVector())
        return T*a0, T*a1, T*b0, T*b1, l

    #===========================================================================
    # Non-overlapping, not parallel
    #===========================================================================
    def testNonOverlappingNonParallelSegments(self):
        a0 = Vector(1.0, 1.0, 1.0)
        a1 = Vector(2.0, 2.0, 2.0)
        b0 = Vector(2.0, 2.0, 3.0)
        b1 = Vector(5.0, 5.0, 3.0)
        answer = 1.0
        result1, result2 = Vector(), Vector()
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, l = self.randomDistortion(a0, a1, b0, b1)
            result = segmentSegmentDistance(aa0, aa1, bb0, bb1)
            self.assertTrue(fuzzyEqual(result, l*answer),
                            "Distance error:  %g != %g" % (result, l*answer))

    #===========================================================================
    # Non-overlapping, parallel
    #===========================================================================
    def testNonOverlappingParallelSegments1(self):
        a0 = Vector(1.0, 1.0, 1.0)
        a1 = Vector(2.0, 2.0, 2.0)
        b0 = a0 + Vector(10.0, 10.0, 1.0)
        b1 = a1 + Vector(10.0, 10.0, 1.0)
        answer = (b0 - a1).magnitude()
        result1, result2 = Vector(), Vector()
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, l = self.randomDistortion(a0, a1, b0, b1)
            result = segmentSegmentDistance(aa0, aa1, bb0, bb1)
            self.assertTrue(fuzzyEqual(result, l*answer),
                            "Distance error:  %g != %g" % (result, l*answer))

    #===========================================================================
    # Non-overlapping, parallel, but overlapping in a projected sense.
    #===========================================================================
    def testNonOverlappingParallelSegments2(self):
        a0 = Vector(1.0, 1.0, 1.0)
        a1 = Vector(2.0, 2.0, 2.0)
        b0 = a0 + Vector(0.0, 0.0, 10.0)
        b1 = a1 + Vector(0.0, 0.0, 10.0)
        answer = (b0 - a1).magnitude()
        result1, result2 = Vector(), Vector()
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, l = self.randomDistortion(a0, a1, b0, b1)
            result = segmentSegmentDistance(aa0, aa1, bb0, bb1)
            self.assertTrue(fuzzyEqual(result, l*answer),
                            "Distance error:  %g != %g" % (result, l*answer))

    #===========================================================================
    # Not parallel, non-intersecting.
    #===========================================================================
    def testNonOverlappingNonParallelSegments1(self):
        a0 = Vector(1.0, 1.0, 1.0)
        a1 = Vector(2.0, 2.0, 2.0)
        b0 = a0 + Vector(3.0, 3.0, 3.0)
        b1 = a1 + Vector(10.0, -10.0, 14.0)
        answer = (b0 - a1).magnitude()
        result1, result2 = Vector(), Vector()
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, l = self.randomDistortion(a0, a1, b0, b1)
            result = segmentSegmentDistance(aa0, aa1, bb0, bb1)
            self.assertTrue(fuzzyEqual(result, l*answer),
                            "Distance error:  %g != %g" % (result, l*answer))

    #===========================================================================
    # Not parallel, endpoint of b on middle of a
    #===========================================================================
    def testNonOverlappingIntersectingSegments1(self):
        a0 = Vector(1.0, 1.0, 1.0)
        a1 = Vector(2.0, 2.0, 2.0)
        b0 = 0.5*(a0 + a1)
        b1 = Vector(5.0, 3.0, 2.0)
        answer = 0.0
        result1, result2 = Vector(), Vector()
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, l = self.randomDistortion(a0, a1, b0, b1)
            result = segmentSegmentDistance(aa0, aa1, bb0, bb1)
            self.assertTrue(fuzzyEqual(result, l*answer),
                            "Distance error:  %g != %g" % (result, l*answer))

    #===========================================================================
    # Not parallel, intersecting
    #===========================================================================
    def testNonOverlappingIntersectingSegments2(self):
        a0 = Vector(1.0, 1.0, 1.0)
        a1 = Vector(2.0, 2.0, 2.0)
        b0 = 0.5*(a0 + a1) - Vector(1.0, 1.0, -1.0)
        b1 = 0.5*(a0 + a1) + Vector(1.0, 1.0, -1.0)
        answer = 0.0
        result1, result2 = Vector(), Vector()
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, l = self.randomDistortion(a0, a1, b0, b1)
            result = segmentSegmentDistance(aa0, aa1, bb0, bb1)
            self.assertTrue(fuzzyEqual(result, l*answer),
                            "Distance error:  %g != %g" % (result, l*answer))

if __name__ == "__main__":
    unittest.main()

