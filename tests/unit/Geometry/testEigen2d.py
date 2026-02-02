#ATS:test(SELF, label="GeomTensor eigen values/vectors unit tests")
# Unit tests for the eigen value/vector methods of the GeomTensor (2D)

import unittest
from math import *
from SpheralTestUtilities import fuzzyEqual

from Spheral import *

# Create a global random number generator.
import random
random.seed(37549927891710)
ranrange = 1.0e8

#===============================================================================
# Compute an accuracy criterion based on how degenerate the eigen values are.
#===============================================================================
def degenerateFuzz(eigenvalues):
    assert eigenvalues[0] <= eigenvalues[1]
    normalization = max([abs(x) for x in eigenvalues] + [1.0e-10])
    dx = abs(eigenvalues[1] - eigenvalues[0])/normalization
    assert dx >= 0.0
    return max(1.0e-5, 1.0/(1.0 + 50.0*dx))

#===============================================================================
# Generate a random 2x2 symmetric tensor with known eigen values and eigen
# vectors.
#===============================================================================
def randomSymTensor2d(lam1 = None,
                      lam2 = None):

    if lam1 is None:
        lam1 = random.uniform(-ranrange, ranrange)
    if lam2 is None:
        lam2 = random.uniform(-ranrange, ranrange)

    # Pick random Euler angles.
    theta = random.uniform(0.0, 2.0*pi)

    # Build the rotation matrix of eigen vectors to rotate from the principle to
    # the lab frame (so transpose of what we usually mean)
    R = Tensor2d( cos(theta), sin(theta),
                 -sin(theta), cos(theta))
    assert fuzzyEqual(R.Determinant(), 1.0)
    check = R*R.Transpose()
    for i in range(2):
        for j in range(2):
            if i == j:
                assert fuzzyEqual(check(i,j), 1.0)
            else:
                assert fuzzyEqual(check(i,j), 0.0)

    # Check the eigen vectors.
    vec1 = R.getColumn(0)
    vec2 = R.getColumn(1)
    assert fuzzyEqual(vec1.magnitude(), 1.0)
    assert fuzzyEqual(vec2.magnitude(), 1.0)
    assert fuzzyEqual(vec1.dot(vec2), 0.0)

    # Now put it all together into our final symmetric matrix.
    A = SymTensor2d(lam1, 0.0,
                    0.0, lam2)
    A.rotationalTransform(R)

    # Return the tensor, it's eigen values, and the tensor of eigenvectors.
    return A, Vector2d(lam1, lam2), R

#===============================================================================
# Test class for Tensor2d.eigenValues and Tensor2d.eigenVectors
#===============================================================================
class TestEigenVectors(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.ntests = 10000
        return

    #---------------------------------------------------------------------------
    # eigenValues (random input)
    #---------------------------------------------------------------------------
    def testRandomEigenValues(self):
        for i in range(self.ntests):
            A, vlam0, vectors0 = randomSymTensor2d()
            lam0 = [x for x in vlam0]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-10),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenValues (two equal eigenvalues)
    #---------------------------------------------------------------------------
    def testDoublyDegenerateEigenValues(self):
        for i in range(self.ntests):
            lam12 = random.uniform(-ranrange, ranrange)
            A, vlam0, vectors0 = randomSymTensor2d(lam1 = lam12,
                                                   lam2 = lam12)
            lam0 = [x for x in vlam0]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-8),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenValues (diagonal matrix input)
    #---------------------------------------------------------------------------
    def testDiagonalEigenValues(self):
        for i in range(self.ntests):
            lam1 = random.uniform(-ranrange, ranrange)
            lam2 = random.uniform(-ranrange, ranrange)
            A = SymTensor2d(lam1, 0.0,
                            0.0, lam2)
            lam0 = [lam1, lam2]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-10),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenVectors (random input)
    #---------------------------------------------------------------------------
    def testRandomEigenVectors(self):
        for i in range(self.ntests):
            A, vlam0, vectors0 = randomSymTensor2d()
            lam0 = [(vlam0(i), vectors0.getColumn(i)) for i in range(2)]
            lam0.sort()
            eigenVecs0 = [x[1] for x in lam0]
            eigenStruct = A.eigenVectors()
            lam = [(eigenStruct.eigenValues(i), eigenStruct.eigenVectors.getColumn(i)) for i in range(2)]
            lam.sort()
            eigenVecs = [x[1] for x in lam]
            for i in range(len(lam0)):
                lami = lam0[i]
                veci = eigenVecs[i]
                vec0 = eigenVecs0[i]
                self.assertTrue(fuzzyEqual(veci.magnitude(), 1.0),
                                "Eigen vector %s does not have unit magnitude" % str(veci))
                self.assertTrue(fuzzyEqual(abs(veci.dot(vec0)), 1.0, degenerateFuzz([x[0] for x in lam0])),
                                "Eigen vector %s does not equal expected value %s for eigen values %s" % (str(veci), str(vec0), str(vlam0)))
        return

    #---------------------------------------------------------------------------
    # eigenVectors (two equal eigenvalues)
    #---------------------------------------------------------------------------
    def testDoublyDegenerateEigenVectors(self):
        for i in range(self.ntests):
            lam12 = random.uniform(-ranrange, ranrange)
            A, vlam0, vectors0 = randomSymTensor2d(lam1 = lam12,
                                                   lam2 = lam12)
            lam0 = [(vlam0(i), vectors0.getColumn(i)) for i in range(2)]
            lam0.sort()
            eigenVecs0 = [x[1] for x in lam0]
            eigenStruct = A.eigenVectors()
            lam = [(eigenStruct.eigenValues(i), eigenStruct.eigenVectors.getColumn(i)) for i in range(2)]
            lam.sort()
            eigenVecs = [x[1] for x in lam]

            for x in eigenVecs:
                self.assertTrue(fuzzyEqual(x.magnitude(), 1.0),
                                "Eigen vector %s does not have unit magnitude %s" % (str(x), str(eigenStruct.eigenVectors)))

            # The eigen vectors need only be perpendicular to each other
            self.assertTrue(fuzzyEqual(eigenVecs[0].dot(eigenVecs[1]), 0.0),
                            "Eigen vectors (%s, %s) are not orthogonal\n%s" % (str(eigenVecs[0]),
                                                                               str(eigenVecs[1]),
                                                                               str(eigenStruct.eigenValues)))


        return

    #---------------------------------------------------------------------------
    # eigenVectors (diagonal matrix input)
    #---------------------------------------------------------------------------
    def testDiagonalEigenVectors(self):
        for i in range(self.ntests):
            lam1 = random.uniform(-ranrange, ranrange)
            lam2 = random.uniform(-ranrange, ranrange)
            A = SymTensor2d(lam1, 0.0,
                            0.0, lam2)
            lam0 = [(lam1, Vector2d(1, 0)),
                    (lam2, Vector2d(0, 1))]
            lam0.sort()
            eigenVecs0 = [x[1] for x in lam0]
            eigenStruct = A.eigenVectors()
            lam = [(eigenStruct.eigenValues(i), eigenStruct.eigenVectors.getColumn(i)) for i in range(2)]
            lam.sort()
            eigenVecs = [x[1] for x in lam]
            for (x, x0) in zip(eigenVecs, eigenVecs0):
                self.assertTrue(fuzzyEqual(x.magnitude(), 1.0),
                                "Eigen vector %s does not equal expected value %s" % (str(x), str(x0)))
                self.assertTrue(fuzzyEqual(abs(x.dot(x0)), 1.0),
                                "Eigen vector %s does not equal expected value %s" % (str(x), str(x0)))
        return

if __name__ == "__main__":
    unittest.main()
