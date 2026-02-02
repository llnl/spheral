#ATS:test(SELF, label="2-D line segment intersections unit tests")
from Spheral2d import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
random.seed(4599281940)

#===============================================================================
# Test our various segement-segment intersection scenarios.
#===============================================================================
class TestSegmentSegmentIntersection(unittest.TestCase):

    #===========================================================================
    # Set up, create arrays of the function values.
    #===========================================================================
    def setUp(self):
        self.ntests = 100
        self.multMin = 0.001
        self.multMax = 1.0e5
        return

    #===========================================================================
    # Randomly distort two line segments.
    #===========================================================================
    def randomDistortion(self, a0, a1, b0, b1):
        T = (random.uniform(self.multMin, self.multMax)*
             rotationMatrix(Vector(random.uniform(0.0, 1.0),
                                   random.uniform(0.0, 1.0)).unitVector()))
        return T*a0, T*a1, T*b0, T*b1, T

    #===========================================================================
    # Non-overlapping
    #===========================================================================
    def testNonOverlappingSegments(self):
        a0 = Vector(1.0, 1.0)
        a1 = Vector(2.0, 2.0)
        b0 = Vector(1.0, 2.0)
        b1 = Vector(2.0, 3.0)
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, T = self.randomDistortion(a0, a1, b0, b1)
            assert segmentSegmentIntersection(aa0, aa1, bb0, bb1)[0] == '0'

    #===========================================================================
    # Intersecting 
    #===========================================================================
    def testIntersectingSegments(self):
        a0 = Vector(1.0, 1.0)
        a1 = Vector(2.0, 2.0)
        b0 = Vector(1.0, 2.0)
        b1 = Vector(2.0, 1.0)
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, T = self.randomDistortion(a0, a1, b0, b1)
            code, result1, result2 = segmentSegmentIntersection(aa0, aa1, bb0, bb1)
            assert code == '1'
            assert result1 == result2
            assert fuzzyEqual((result1 - T*Vector(1.5, 1.5)).magnitude(), 0.0)

    #===========================================================================
    # Intersecting at endpoint
    #===========================================================================
    def testIntersectingSegmentsOnEndpoint(self):
        a0 = Vector(1.0, 1.0)
        a1 = Vector(2.0, 2.0)
        b0 = Vector(1.0, 2.0)
        b1 = Vector(1.5, 1.5)
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, T = self.randomDistortion(a0, a1, b0, b1)
            code, result1, result2 = segmentSegmentIntersection(aa0, aa1, bb0, bb1)
            assert code == 'v'
            assert result1 == result2
            assert fuzzyEqual((result1 - T*Vector(1.5, 1.5)).magnitude(), 0.0)

    #===========================================================================
    # Overlapping segments.
    #===========================================================================
    def testOverlappingSegments(self):
        a0 = Vector(1.0, 1.0)
        a1 = Vector(2.0, 2.0)
        b0 = Vector(3.0, 3.0)
        b1 = Vector(1.5, 1.5)
        for i in range(self.ntests):
            aa0, aa1, bb0, bb1, T = self.randomDistortion(a0, a1, b0, b1)
            code, result1, result2 = segmentSegmentIntersection(aa0, aa1, bb0, bb1)
            assert code == 'e'
            assert result1 != result2
            if result1.magnitude2() > result2.magnitude2():
                result1, result2 = result2, result1
            assert fuzzyEqual((result1 - T*Vector(1.5, 1.5)).magnitude(), 0.0)
            assert fuzzyEqual((result2 - T*Vector(2.0, 2.0)).magnitude(), 0.0)

if __name__ == "__main__":
    unittest.main()

