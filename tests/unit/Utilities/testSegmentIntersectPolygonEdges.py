#ATS:test(SELF, label="2-D line segment intersecting edges of a polygon.")
from math import *
from Spheral2d import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
random.seed(4599281940)


#===============================================================================
# Test whether a line segment intersects a polygon.
#===============================================================================
class TestLineSegmentPolygonIntersection(unittest.TestCase):

    #===========================================================================
    # Randomly distort two line segments.
    #===========================================================================
    def randomDistortion(self, a0, a1, vertices):
        l = random.uniform(self.multMin, self.multMax)
        T = l*rotationMatrix(Vector(random.uniform(0.0, 1.0),
                                    random.uniform(0.0, 1.0),
                                    random.uniform(0.0, 1.0)).unitVector())
        verts = vector_of_Vector()
        for x in vertices:
            verts.append(T*x)
        return T*a0, T*a1, Polygon(verts), T

    #===========================================================================
    # setUp
    #===========================================================================
    def setUp(self):
        self.ntests = 1000
        self.multMin = 0.001
        self.multMax = 1.0e5
        self.vertices = [Vector(1.0, 1.0), Vector(2.0, 1.0),
                         Vector(2.0, 2.0), Vector(1.0, 2.0)]
        return

    #===========================================================================
    # Segment entirely outside the polygon
    #===========================================================================
    def testNonintersectingSegment1(self):
        a0 = Vector(10.0, 10.0)
        a1 = Vector(20.0, 0.0)
        for i in range(self.ntests):
            aa0, aa1, polygon, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polygon)
            self.assertTrue(result == False,
                            "Incorrectly intersected edge %s->%s with polygon" %
                            (aa0, aa1))

    #===========================================================================
    # Segment entirely inside the polygon
    #===========================================================================
    def testNonintersectingSegment2(self):
        a0 = Vector(1.25, 1.25)
        a1 = Vector(1.75, 1.75)
        for i in range(self.ntests):
            aa0, aa1, polygon, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polygon)
            Tinv = T.Inverse()
            self.assertTrue(result == False,
                            "Incorrectly intersected edge %s->%s with polygon" %
                            (Tinv*aa0, Tinv*aa1))

    #===========================================================================
    # Segment intersecting a random side of the polygon
    #===========================================================================
    def testSegmentIntersectingRandomEdge1(self):
        a0 = Vector(1.5, 1.5)
        for i in range(self.ntests):
            theta = random.uniform(0.0, 2.0*pi)
            a1 = a0 + Vector(cos(theta), sin(theta))
            aa0, aa1, polygon, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polygon)
            Tinv = T.Inverse()
            self.assertTrue(result == True,
                            "Incorrectly missed intersection of edge %s->%s with polygon" %
                            (Tinv*aa0, Tinv*aa1))

    #===========================================================================
    # Interior segment with an endpoint on a random point of the polygon
    #===========================================================================
    def testSegmentIntersectingRandomEdge2(self):
        a0 = Vector(1.5, 1.5)
        deltas = []
        nverts = len(self.vertices)
        for i in range(nverts):
            j = (i + 1) % nverts
            deltas.append(self.vertices[j] - self.vertices[i])
        assert len(deltas) == nverts
        for i in range(self.ntests):
            j = random.randint(0, nverts - 1)
            a1 = self.vertices[j] + random.random()*deltas[j]
            aa0, aa1, polygon, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polygon)
            Tinv = T.Inverse()
            self.assertTrue(result == True,
                            "Incorrectly missed intersection of edge %s->%s with polygon" %
                            (Tinv*aa0, Tinv*aa1))

if __name__ == "__main__":
    unittest.main()

