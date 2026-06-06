#ATS:test(SELF, label="3-D line segment intersecting edges of a polyhedron.")
from math import *
from Spheral3d import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
random.seed(4599281940)

#===============================================================================
# Test whether a line segment intersects a polyhedron.
#===============================================================================
class TestLineSegmentPolyhedronIntersection(unittest.TestCase):

    #===========================================================================
    # Randomly distort two line segments.
    #===========================================================================
    def randomDistortion(self, a0, a1, vertices):
        l = random.uniform(self.multMin, self.multMax)
        theta = random.uniform(0.0, 2.0*pi)
        phi = random.uniform(0.0, pi)
        T = l*rotationMatrix(Vector(cos(theta)*sin(phi),
                                    sin(theta)*sin(phi),
                                    cos(phi)))
        verts = vector_of_Vector()
        for x in vertices:
            verts.append(T*x)
        return T*a0, T*a1, Polyhedron(verts), T

    #===========================================================================
    # setUp
    #===========================================================================
    def setUp(self):
        self.ntests = 1000
        self.multMin = 0.001
        self.multMax = 1.0e5
        self.vertices = [Vector(1.0, 1.0, 1.0), Vector(2.0, 1.0, 1.0),
                         Vector(2.0, 2.0, 1.0), Vector(1.0, 2.0, 1.0),
                         Vector(1.0, 1.0, 2.0), Vector(2.0, 1.0, 2.0),
                         Vector(2.0, 2.0, 2.0), Vector(1.0, 2.0, 2.0)]
        self.edgeVertices = [(0,1), (1,2), (2,3), (3,0),
                             (0,4), (1,5), (2,6), (3,7),
                             (4,5), (5,6), (6,7), (7,4)]
        self.edges = [(self.vertices[e[0]], self.vertices[e[1]] - self.vertices[e[0]]) for e in self.edgeVertices]
        self.faces = [(Vector(1.5, 1.5, 1.0), 0, 1),
                      (Vector(1.0, 1.5, 1.5), 1, 2),
                      (Vector(1.5, 1.0, 1.5), 0, 2),
                      (Vector(2.0, 1.5, 1.5), 1, 2),
                      (Vector(1.5, 2.0, 1.5), 0, 2),
                      (Vector(1.5, 1.5, 2.0), 0, 1)]

        return

    #===========================================================================
    # Segment entirely outside the polyhedron
    #===========================================================================
    def testNonintersectingSegment1(self):
        a0 = Vector(10.0, 10.0, 10.0)
        a1 = Vector(20.0, 20.0, 20.0)
        for i in range(self.ntests):
            aa0, aa1, polyhedron, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polyhedron)
            self.assertTrue(result == False,
                            "Incorrectly intersected edge %s->%s with polyhedron" %
                            (aa0, aa1))

    #===========================================================================
    # Segment entirely inside the polyhedron
    #===========================================================================
    def testNonintersectingSegment2(self):
        a0 = Vector(1.25, 1.25, 1.25)
        a1 = Vector(1.75, 1.75, 1.75)
        for i in range(self.ntests):
            aa0, aa1, polyhedron, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polyhedron)
            Tinv = T.Inverse()
            self.assertTrue(result == False,
                            "Incorrectly intersected edge %s->%s with polyhedron" %
                            (Tinv*aa0, Tinv*aa1))

    #===========================================================================
    # Segment passing through two polyhedron faces without hitting an edge
    #===========================================================================
    def testNonintersectingSegment3(self):
        for i in range(self.ntests):
            faces = random.sample(self.faces, 2)
            theta0, theta1 = random.uniform(0.0, 2.0*pi), random.uniform(0.0, 2.0*pi)
            l0, l1 = random.uniform(0.0, 0.49), random.uniform(0.0, 0.49)

            a0 = Vector(faces[0][0])
            a0[faces[0][1]] += l0*cos(theta0)
            a0[faces[0][2]] += l0*sin(theta0)
            
            a1 = Vector(faces[1][0])
            a1[faces[1][1]] += l1*cos(theta1)
            a1[faces[1][2]] += l1*sin(theta1)

            assert min([x >= 1.0 and x <= 2.0 for x in a0]) == True
            assert min([x >= 1.0 and x <= 2.0 for x in a1]) == True

            aa0, aa1, polyhedron, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polyhedron)
            Tinv = T.Inverse()
            self.assertTrue(result == False,
                            "Incorrectly intersected edge %s->%s with polyhedron\n%s->%s" %
                            (aa0, aa1, Tinv*aa0, Tinv*aa1))

    #===========================================================================
    # Segment intersecting a random edge of the polyhedron
    #===========================================================================
    def testSegmentIntersectingRandomEdge1(self):
        a0 = Vector(1.5, 1.5, 1.5)
        for i in range(self.ntests):
            edge = random.choice(self.edges)
            a1 = edge[0] + random.random()*edge[1]
            a1 = a0 + 2.0*(a1 - a0)
            aa0, aa1, polyhedron, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polyhedron)
            Tinv = T.Inverse()
            self.assertTrue(result == True,
                            "Incorrectly missed intersection of edge %s->%s with polyhedron" %
                            (Tinv*aa0, Tinv*aa1))

    #===========================================================================
    # Segment with an endpoint on a random edge of the polyhedron
    #===========================================================================
    def testSegmentIntersectingRandomEdge2(self):
        a0 = Vector(1.5, 1.5, 1.5)
        for i in range(self.ntests):
            edge = random.choice(self.edges)
            a1 = edge[0] + random.random()*edge[1]
            aa0, aa1, polyhedron, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polyhedron)
            Tinv = T.Inverse()
            self.assertTrue(result == True,
                            "Incorrectly missed intersection of edge %s->%s with polyhedron" %
                            (Tinv*aa0, Tinv*aa1))

    #===========================================================================
    # Segment with an endpoint on a random vertex of the polyhedron
    #===========================================================================
    def testSegmentIntersectingRandomVertex(self):
        a0 = Vector(1.5, 1.5, 1.5)
        for i in range(self.ntests):
            a1 = random.choice(self.vertices)
            aa0, aa1, polyhedron, T = self.randomDistortion(a0, a1, self.vertices)
            result = segmentIntersectEdges(aa0, aa1, polyhedron)
            Tinv = T.Inverse()
            self.assertTrue(result == True,
                            "Incorrectly missed intersection of edge %s->%s with polyhedron" %
                            (Tinv*aa0, Tinv*aa1))

if __name__ == "__main__":
    unittest.main()

