#ATS:test(SELF, label="Polyhedron unit tests")
# Unit tests for the Polyhedron class

import unittest
from math import *
import time
from SpheralTestUtilities import fuzzyEqual

from Spheral3d import *

# Create a global random number generator.
import random
random.seed(630)

#===============================================================================
# Generate random points in the give box, optionally rotating the results to
# create a non-axis aligned distribution.
#===============================================================================
def randomPoints(numPoints,
                 xmin, xmax,
                 ymin, ymax,
                 zmin, zmax,
                 theta = None,
                 phi = None):
    result = [] # vector_of_Vector()

    # Determine the rotational transform.
    if theta is None:
        theta = random.uniform(0.0, 2.0*pi)
    if phi is None:
        phi = random.uniform(0.0, pi)
    R = rotationMatrix(Vector(cos(theta)*sin(phi),
                              sin(theta)*sin(phi),
                              cos(phi)))
    
    for i in range(numPoints):
        result.append(R*Vector(random.uniform(xmin, xmax),
                               random.uniform(ymin, ymax),
                               random.uniform(zmin, zmax)))

    return R, result

#===============================================================================
# Return a random vector with at most the given magnitude.
#===============================================================================
def randomVector(rmin, rmax):
    xhat = random.uniform(0.0, 1.0)
    yhat = random.uniform(0.0, sqrt(1.0 - xhat*xhat))
    zhat = sqrt(1.0 - xhat*xhat - yhat*yhat)
    r = random.uniform(rmin, rmax)
    return Vector(r*xhat, r*yhat, r*zhat)

#===============================================================================
# Test class
#===============================================================================
class TestPolyhedron(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.ntests = 5000
        self.npoints = 100000
        self.xmin, self.xmax = -5.0, 2.0
        self.ymin, self.ymax = -1.0, 3.0
        self.zmin, self.zmax = -3.0, 10.0
        self.R, self.points = randomPoints(self.npoints,
                                           self.xmin, self.xmax,
                                           self.ymin, self.ymax,
                                           self.zmin, self.zmax)
        #t0 = time.clock()
        self.polyhedron = Polyhedron(self.points)
        #t1 = time.clock()
        #print "Required %s seconds to generate polyhedron." % (t1 - t0)
        return

    #---------------------------------------------------------------------------
    # Find the inner & outer radii of the polyhedron.
    #---------------------------------------------------------------------------
    def innerOuterRadii(self, polyhedron):
        rinner = 1e10
        router = 0.0
        centroid = polyhedron.centroid
        for f in self.polyhedron.facets:
            r = abs((f.point(0) - centroid).dot(f.normal))
            rinner = min(rinner, r)
        for v in self.polyhedron.vertices:
            r = (v - centroid).magnitude()
            router = max(router, r)
        router *= 1.0 + 1.0e-5
        return rinner, router

    #---------------------------------------------------------------------------
    # Check that a random subset the seed points are contained in the polyhedron
    # with the normal contain method.
    #---------------------------------------------------------------------------
    def testContainSeeds(self):
        for p in self.points: # random.sample(self.points, 10000):
            result = self.polyhedron.contains(p)
            if not result:
                print("Bad polyhedron:  ", [str(x) for x in self.polyhedron.vertices])
                print("Test if point on polyhedron:  ", pointOnPolyhedron(p, self.polyhedron))
            self.assertTrue(result,
                            "Polyhedron does not contain seed point: %s" % str(p))
            self.assertTrue(result == self.polyhedron.contains(p, useAxom=True),
                            "Polyhedron contains inconsistent using Axom")
        return

    #---------------------------------------------------------------------------
    # Check that a random subset the seed points are contained in the polyhedron
    # with the generic contain method.
    #---------------------------------------------------------------------------
    def testGenericContainSeeds(self):
        for p in random.sample(self.points, 5000):
            result = pointInPolyhedron(p, self.polyhedron, True)
            if not result:
                print("Bad polyhedron:  ", [str(x) for x in self.polyhedron.vertices])
                print("Distance from polyhedron: ", self.polyhedron.distance(p))
                print("Test if point on polyhedron:  ", pointOnPolyhedron(p, self.polyhedron))
                print("Min zray distance:  ", min([(Vector2d(p.x, p.y) - Vector2d(vi.x, vi.y)).magnitude() for vi in self.polyhedron.vertices]))
            self.assertTrue(result,
                            "Polyhedron does not contain seed point: %s" % str(p))
            self.assertTrue(result == self.polyhedron.contains(p, useAxom=True),
                            "Polyhedron contains inconsistent using Axom")
        return

    #---------------------------------------------------------------------------
    # Check that all the seed points are contained in the polyhedron using
    # the convexContain method.
    #---------------------------------------------------------------------------
    def testConvexContainSeeds(self):
        for p in self.points:
            result = self.polyhedron.convexContains(p)
            self.assertTrue(result,
                            "Polyhedron does not contain seed point: %s" % str(p))
            self.assertTrue(result == self.polyhedron.contains(p, useAxom=True),
                            "Polyhedron contains inconsistent using Axom")
        return

    #---------------------------------------------------------------------------
    # Generate random points in the polyhedron and test they are contained.
    #---------------------------------------------------------------------------
    def testRandomInnerPoints(self):
        rinner, router = self.innerOuterRadii(self.polyhedron)
        centroid = self.polyhedron.centroid
        for i in range(self.ntests):
            p = centroid + randomVector(0.0, rinner)
            result = self.polyhedron.contains(p)
            self.assertTrue(result,
                            "Polyhedron should contain %s but reports it does not." % str(p))
            self.assertTrue(result == self.polyhedron.contains(p, useAxom=True),
                            "Polyhedron contains inconsistent using Axom")
        return

    #---------------------------------------------------------------------------
    # Generate random points outside the polyhedron and test they are not
    # contained.
    #---------------------------------------------------------------------------
    def testRandomOuterPoints(self):
        rinner, router = self.innerOuterRadii(self.polyhedron)
        centroid = self.polyhedron.centroid
        for i in range(self.ntests):
            p = centroid + randomVector(router, 2.0*router)
            result = self.polyhedron.contains(p)
            self.assertTrue(not result,
                            "%s should be outside polyhedron but polyhedron reports it is contained." % str(p))
            self.assertTrue(result == self.polyhedron.contains(p, useAxom=True),
                            "Polyhedron contains inconsistent using Axom")
        return

    #---------------------------------------------------------------------------
    # Test vertex containment.
    #---------------------------------------------------------------------------
    def testVertexPointContainment(self):
        vertices = self.polyhedron.vertices
        for v in vertices:
            self.assertTrue(self.polyhedron.contains(v),
                            "%s vertex position should be contained." % str(v))
            self.assertTrue(not self.polyhedron.contains(v, False),
                            "%s vertex position should not be contained." % str(v))
        return

    #---------------------------------------------------------------------------
    # Test that nested polyhedrons intersect.
    #---------------------------------------------------------------------------
    def testIntersectInterior(self):
        centroid = self.polyhedron.centroid
        vertices = vector_of_Vector(self.polyhedron.vertices)
        numVerts = len(vertices)
        for i in range(numVerts):
            vertices[i] = 0.5*(centroid + vertices[i])
            assert self.polyhedron.contains(vertices[i])
        polyhedron2 = Polyhedron(vertices)
        self.assertTrue(self.polyhedron.intersect(polyhedron2),
                        "Failed to intersect with a contained polyhedron.")
        self.assertTrue(polyhedron2.intersect(self.polyhedron),
                        "Failed to intersect when entirely contained within a polyhedron.")
        return

    #---------------------------------------------------------------------------
    # Test that two polyhedrons just touching intersect.
    #---------------------------------------------------------------------------
    def testIntersectTouchingPolyhedrons(self):
        xmin = self.polyhedron.xmin.x
        vertices = vector_of_Vector(self.polyhedron.vertices)
        numVerts = len(vertices)
        for i in range(numVerts):
            xv = vertices[i].x
            vertices[i].x -= 2.0*(xv - xmin)
        polyhedron2 = Polyhedron(vertices)
        self.assertTrue(self.polyhedron.intersect(polyhedron2),
                        "Failed to intersect with polyhedron touching at a vertex")
        self.assertTrue(polyhedron2.intersect(self.polyhedron),
                        "Failed to intersect with polyhedron touching at a vertex")
        return

    #---------------------------------------------------------------------------
    # Test that two separate polyhedrons do not intersect.
    #---------------------------------------------------------------------------
    def testNotIntersectPolyhedrons(self):
        xmin = self.polyhedron.xmin.x
        xlength = self.polyhedron.xmax.x - self.polyhedron.xmin.x
        vertices = vector_of_Vector(self.polyhedron.vertices)
        numVerts = len(vertices)
        for i in range(numVerts):
            xv = vertices[i].x
            vertices[i].x -= 2.0*(xv - xmin) + xlength
        polyhedron2 = Polyhedron(vertices)
        self.assertTrue(not self.polyhedron.intersect(polyhedron2),
                        "Erroneously claiming polyhedrons intersect : [[%g,%g], [%g,%g], [[%g,%g], [%g,%g]]" % (self.polyhedron.xmin.x, 
                                                                                                             self.polyhedron.xmax.x,
                                                                                                             self.polyhedron.xmin.y,
                                                                                                             self.polyhedron.xmax.y,
                                                                                                             polyhedron2.xmin.x, 
                                                                                                             polyhedron2.xmax.x,
                                                                                                             polyhedron2.xmin.y,
                                                                                                             polyhedron2.xmax.y))
        self.assertTrue(not polyhedron2.intersect(self.polyhedron),
                        "Erroneously claiming polyhedrons intersect : [[%g,%g], [%g,%g], [[%g,%g], [%g,%g]]" % (self.polyhedron.xmin.x, 
                                                                                                             self.polyhedron.xmax.x,
                                                                                                             self.polyhedron.xmin.y,
                                                                                                             self.polyhedron.xmax.y,
                                                                                                             polyhedron2.xmin.x, 
                                                                                                             polyhedron2.xmax.x,
                                                                                                             polyhedron2.xmin.y,
                                                                                                             polyhedron2.xmax.y))
        return

    #---------------------------------------------------------------------------
    # Test that a nested box intersects.
    #---------------------------------------------------------------------------
    def testIntersectBoxInterior(self):
        rinner, router = self.innerOuterRadii(self.polyhedron)
        centroid = self.polyhedron.centroid
        box = (centroid - 0.95*rinner*Vector.one,
               centroid + 0.95*rinner*Vector.one)
        self.assertTrue(self.polyhedron.intersect(box),
                        "Failed to intersect with a contained box.")
        return

    #---------------------------------------------------------------------------
    # Test that a circumscribing box intersects.
    #---------------------------------------------------------------------------
    def testIntersectBoxExterior(self):
        rinner, router = self.innerOuterRadii(self.polyhedron)
        centroid = self.polyhedron.centroid
        box = (centroid - 2.0*router*Vector.one,
               centroid + 2.0*router*Vector.one)
        self.assertTrue(self.polyhedron.intersect(box),
                        "Failed to intersect with a box we are in.")
        return
   
    #---------------------------------------------------------------------------
    # Test that a box just touching a polyhedron intersects.
    #---------------------------------------------------------------------------
    def testIntersectTouchingBox(self):
        xmin = self.polyhedron.xmin.x
        vertex = None
        for v in self.polyhedron.vertices:
            if fuzzyEqual(v.x, xmin, 1.0e-10):
                vertex = v
        assert not vertex is None
        box = (Vector(v.x - 2.0, v.y - 2.0, v.z - 2.0), v)
        self.assertTrue(self.polyhedron.intersect(box),
                        "Failed to intersect with a box touching on one side.")

    #---------------------------------------------------------------------------
    # Test that that we don't intersect a box outside the polyhedron.
    #---------------------------------------------------------------------------
    def testNotIntersectBox(self):
        xmin, xmax = self.polyhedron.xmin, self.polyhedron.xmax
        delta = xmax - xmin
        box = (xmin - 2.0*delta, xmax - 2.0*delta)
        self.assertTrue(not self.polyhedron.intersect(box),
                        "Erroneously intersect external box.")

    #---------------------------------------------------------------------------
    # Test reconstructing.
    #---------------------------------------------------------------------------
    def testReconstruct(self):
        polyhedron2 = Polyhedron()
        polyhedron2.reconstruct(self.polyhedron.vertices,
                                self.polyhedron.facetVertices)
        self.assertTrue(polyhedron2 == self.polyhedron,
                        "Failed to properly reconstruct polyhedron from vertices and facets.")
        return

    #---------------------------------------------------------------------------
    # Test volume
    #---------------------------------------------------------------------------
    def testVolume(self):
        facets = self.polyhedron.facets
        c = self.polyhedron.centroid
        answer = 0.0
        for f in facets:
            answer += f.area * (f.normal.dot(f.point(0) - c))
        answer /= 3.0
        self.assertTrue(fuzzyEqual(self.polyhedron.volume, answer, 1.0e-10),
                        "Failed volume computation: %g != %g" % (self.polyhedron.volume,
                                                                 answer))

    #---------------------------------------------------------------------------
    # Closest point to vertices.
    #---------------------------------------------------------------------------
    def testClosestPointToVertices(self):
        verts = self.polyhedron.vertices
        for p in verts:
            cp = self.polyhedron.closestPoint(p)
            self.assertTrue(fuzzyEqual((cp - p).magnitude(), 0.0, 1.0e-10),
                            "Closest point to vertex %s : %s" % (p, cp))

    #---------------------------------------------------------------------------
    # Closest point on facets.
    #---------------------------------------------------------------------------
    def testClosestPointOnFacets(self):
        facets = self.polyhedron.facets
        for f in facets:
            p = f.position
            cp = self.polyhedron.closestPoint(p)
            self.assertTrue(fuzzyEqual((cp - p).magnitude(), 0.0, 1.0e-10),
                            "Closest point to facet position %s : %s" % (p, cp))

    #---------------------------------------------------------------------------
    # Closest point above facets.
    #---------------------------------------------------------------------------
    def testClosestPointAboveFacets(self):
        facets = self.polyhedron.facets
        verts = self.polyhedron.vertices
        for f in facets:
            iverts = f.ipoints
            n = len(iverts)
            assert n >= 3
            minedge = 1e100
            for k in range(len(iverts)):
                i0, i1 = iverts[k], iverts[(k + 1) % n]
                minedge = min(minedge, (verts[i1] - verts[i0]).magnitude())
            #chi = random.uniform(0.1, 10.0)
            cp0 = f.position
            p = cp0 + 0.5*minedge*f.normal
            cp = self.polyhedron.closestPoint(p)
            self.assertTrue(fuzzyEqual((cp0 - cp).magnitude(), 0.0, 1.0e-10),
                            "Closest point to position off of facet position %s : %s" % (cp0, cp))

    #---------------------------------------------------------------------------
    # Test ==
    #---------------------------------------------------------------------------
    def testEqual(self):
        self.assertTrue(self.polyhedron == self.polyhedron,
                        "Failed self equivalence.")

    #---------------------------------------------------------------------------
    # Test !=
    #---------------------------------------------------------------------------
    def testNotEqual(self):
        polyhedron2 = Polyhedron()
        self.assertTrue(polyhedron2 != self.polyhedron,
                        "Failed not equal.")

    #---------------------------------------------------------------------------
    # Test copy constructor
    #---------------------------------------------------------------------------
    def testCopy(self):
        polyhedron2 = Polyhedron(self.polyhedron)
        self.assertTrue(polyhedron2 == self.polyhedron,
                        "Failed to copy construct.")

    #---------------------------------------------------------------------------
    # Test shift in-place
    #---------------------------------------------------------------------------
    def testShiftInPlace(self):
        shift = Vector(random.uniform(-10.0, -10.0),
                       random.uniform(-10.0, -10.0))
        polyhedron2 = Polyhedron(self.polyhedron)
        polyhedron2 += shift
        for p0, p1 in zip([self.polyhedron.xmin, self.polyhedron.xmax] + list(self.polyhedron.vertices),
                          [polyhedron2.xmin, polyhedron2.xmax] + list(polyhedron2.vertices)):
            pshift = p0 + shift
            self.assertTrue(pshift == p1, "In-place shift point comparison failed: %s != %s" % (pshift, p1))

if __name__ == "__main__":
    unittest.main()
