#ATS:test(SELF, label="GeomTensor eigen values/vectors unit tests")
# Unit tests for the eigen value/vector methods of the GeomTensor.

import unittest
from math import *
from SpheralTestUtilities import fuzzyEqual

from Spheral import *

# Create a global random number generator.
import random
random.seed(3754589209)
ranrange = 1.0e8

#===============================================================================
# Compute an accuracy criterion based on how degenerate the eigen values are.
#===============================================================================
def degenerateFuzz(i, eigenvalues):
    assert eigenvalues[0] <= eigenvalues[1] and eigenvalues[1] <= eigenvalues[2]
    normalization = max([abs(x) for x in eigenvalues] + [1.0e-10])
    if i == 0:
        dx = abs(eigenvalues[1] - eigenvalues[0])/normalization
    else:
        dx = abs(eigenvalues[2] - eigenvalues[1])/normalization
    assert dx >= 0.0
    return max(1.0e-5, 1.0/(1.0 + 50.0*dx))

#===============================================================================
# Generate a random 3x3 symmetric tensor with known eigen values and eigen
# vectors.
#===============================================================================
def randomSymTensor3d(lam1 = None,
                      lam2 = None,
                      lam3 = None):

    if lam1 is None:
        lam1 = random.uniform(-ranrange, ranrange)
    if lam2 is None:
        lam2 = random.uniform(-ranrange, ranrange)
    if lam3 is None:
        lam3 = random.uniform(-ranrange, ranrange)

    # Pick random Euler angles.
    theta = random.uniform(0.0, 2.0*pi)
    phi = random.uniform(0.0, pi)
    psi = random.uniform(0.0, pi)

    # Build the rotation matrix of eigen vectors.
    R = Tensor3d(cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                 -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                 sin(theta)*sin(phi),
                 cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                 -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                 -sin(theta)*cos(phi),
                 sin(theta)*sin(psi),
                 sin(theta)*cos(psi),
                 cos(theta))
    assert fuzzyEqual(R.Determinant(), 1.0)
    check = R*R.Transpose()
    for i in range(3):
        for j in range(3):
            if i == j:
                assert fuzzyEqual(check(i,j), 1.0)
            else:
                assert fuzzyEqual(check(i,j), 0.0)

    # Check the eigen vectors.
    vec1 = R.getColumn(0)
    vec2 = R.getColumn(1)
    vec3 = R.getColumn(2)
    assert fuzzyEqual(vec1.magnitude(), 1.0)
    assert fuzzyEqual(vec2.magnitude(), 1.0)
    assert fuzzyEqual(vec3.magnitude(), 1.0)
    assert fuzzyEqual(vec1.dot(vec2), 0.0)
    assert fuzzyEqual(vec3.dot(vec1), 0.0)
    assert fuzzyEqual(vec3.dot(vec2), 0.0)

    # Now put it all together into our final symmetric matrix.
    A = SymTensor3d(lam1, 0.0, 0.0,
                    0.0, lam2, 0.0,
                    0.0, 0.0, lam3)
    A.rotationalTransform(R)

    # Return the tensor, it's eigen values, and the tensor of eigenvectors.
    return A, Vector3d(lam1, lam2, lam3), R

#===============================================================================
# Test class for Tensor3d.eigenValues and Tensor3d.eigenVectors
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
            A, vlam0, vectors0 = randomSymTensor3d()
            lam0 = [x for x in vlam0]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-5),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenValues (two equal eigenvalues)
    #---------------------------------------------------------------------------
    def testDoublyDegenerateEigenValues(self):
        for i in range(self.ntests):
            lam12 = random.uniform(-ranrange, ranrange)
            A, vlam0, vectors0 = randomSymTensor3d(lam1 = lam12,
                                                   lam2 = lam12)
            lam0 = [x for x in vlam0]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-3),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenValues (three equal eigenvalues)
    #---------------------------------------------------------------------------
    def testTriplyDegenerateEigenValues(self):
        for i in range(self.ntests):
            lam123 = random.uniform(-ranrange, ranrange)
            A, vlam0, vectors0 = randomSymTensor3d(lam1 = lam123,
                                                   lam2 = lam123,
                                                   lam3 = lam123)
            lam0 = [x for x in vlam0]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-5),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenValues (diagonal matrix input)
    #---------------------------------------------------------------------------
    def testDiagonalEigenValues(self):
        for i in range(self.ntests):
            lam1 = random.uniform(-ranrange, ranrange)
            lam2 = random.uniform(-ranrange, ranrange)
            lam3 = random.uniform(-ranrange, ranrange)
            A = SymTensor3d(lam1, 0.0, 0.0,
                            0.0, lam2, 0.0,
                            0.0, 0.0, lam3)
            lam0 = [lam1, lam2, lam3]
            lam0.sort()
            vlam = A.eigenValues()
            lam = [x for x in vlam]
            lam.sort()
            for (x, x0) in zip(lam, lam0):
                self.assertTrue(fuzzyEqual(x, x0, 1e-5),
                                "Eigen values %s do not equal expected values %s" % (str(lam), str(lam0)))
        return

    #---------------------------------------------------------------------------
    # eigenVectors (random input)
    #---------------------------------------------------------------------------
    def testRandomEigenVectors(self):
        for i in range(self.ntests):
            A, vlam0, vectors0 = randomSymTensor3d()
            lam0 = [(vlam0(i), vectors0.getColumn(i)) for i in range(3)]
            lam0.sort()
            eigenVecs0 = [x[1] for x in lam0]
            eigenStruct = A.eigenVectors()
            lam = [(eigenStruct.eigenValues(i), eigenStruct.eigenVectors.getColumn(i)) for i in range(3)]
            lam.sort()
            eigenVecs = [x[1] for x in lam]
            for i in range(len(lam0)):
                lami = lam0[i]
                veci = eigenVecs[i]
                vec0 = eigenVecs0[i]
                self.assertTrue(fuzzyEqual(veci.magnitude(), 1.0),
                                "Eigen vector %s does not have unit magnitude" % str(veci))
                self.assertTrue(fuzzyEqual(abs(veci.dot(vec0)), 1.0, degenerateFuzz(i, [x[0] for x in lam0])),
                                "Eigen vector %s does not equal expected value %s for eigen values %s" % (str(veci), str(vec0), str(vlam0)))
        return

    #---------------------------------------------------------------------------
    # eigenVectors (two equal eigenvalues)
    #---------------------------------------------------------------------------
    def testDoublyDegenerateEigenVectors(self):
        for i in range(self.ntests):
            lam12 = random.uniform(-ranrange, ranrange)
            A, vlam0, vectors0 = randomSymTensor3d(lam1 = lam12,
                                                   lam2 = lam12)
            lam0 = [(vlam0(i), vectors0.getColumn(i)) for i in range(3)]
            lam0.sort()
            eigenVecs0 = [x[1] for x in lam0]
            eigenStruct = A.eigenVectors()
            lam = [(eigenStruct.eigenValues(i), eigenStruct.eigenVectors.getColumn(i)) for i in range(3)]
            lam.sort()
            eigenVecs = [x[1] for x in lam]

            for x in eigenVecs:
                self.assertTrue(fuzzyEqual(x.magnitude(), 1.0),
                                "Eigen vector %s does not have unit magnitude %s" % (str(x), str(eigenStruct.eigenVectors)))

            # Identify the unique eigen value.
            unique = -1
            thpt = 0.0
            degenerate = []
            if (lam0[0][0] == lam0[1][0]):
                unique = 2
                degenerate = [0, 1]
                thpt = abs(vlam0(2))/(abs(vlam0(0)) + 1.0e-10)
            else:
                unique = 0
                degenerate = [1, 2]
                thpt = abs(vlam0(0))/(abs(vlam0(1)) + 1.0e-10)
            assert thpt > 0.0
            if thpt > 1.0:
                thpt = 1.0/thpt

            # Does the eigenvector for the unique eigen value match?
            self.assertTrue(fuzzyEqual(abs(eigenVecs[unique].dot(eigenVecs0[unique])), 1.0, degenerateFuzz(unique, [x[0] for x in lam0])),
                            "Eigen vector %s does not equal expected value %s for eigen values %s, %s" % (str(eigenVecs[unique]),
                                                                                                          str(eigenVecs0[unique]),
                                                                                                          str(vlam0),
                                                                                                          str(eigenStruct.eigenValues)))

            # The remaining eigen values need only be perpendicular to each other and the unique
            # value.
            self.assertTrue(fuzzyEqual(eigenVecs[0].dot(eigenVecs[1]), 0.0) and
                            fuzzyEqual(eigenVecs[0].dot(eigenVecs[2]), 0.0) and
                            fuzzyEqual(eigenVecs[1].dot(eigenVecs[2]), 0.0),
                            "Eigen vectors (%s, %s, %s) are not orthogonal\n%s" % (str(eigenVecs[0]),
                                                                                   str(eigenVecs[1]),
                                                                                   str(eigenVecs[2]),
                                                                                   str(eigenStruct.eigenValues)))


        return

    #---------------------------------------------------------------------------
    # eigenVectors (three equal eigenvalues)
    #---------------------------------------------------------------------------
    def testTriplyDegenerateEigenVectors(self):
        for i in range(self.ntests):
            lam123 = random.uniform(-ranrange, ranrange)
            A = SymTensor3d(lam123, 0.0, 0.0,
                            0.0, lam123, 0.0,
                            0.0, 0.0, lam123)
            vlam0 = Vector3d(lam123, lam123, lam123)
            vectors0 = [Vector3d(1, 0, 0),
                        Vector3d(0, 1, 0),
                        Vector3d(0, 0, 1)]
            eigenStruct = A.eigenVectors()
            for i in range(3):
                vec = eigenStruct.eigenVectors.getColumn(i)
                match = [fuzzyEqual(abs(vec.dot(vectors0[j])), 1.0) for j in range(len(vectors0))]
                assert len(match) == len(vectors0)
                assert sum(match) == 1
                del vectors0[match.index(True)]
            self.assertTrue(len(vectors0) == 0,
                            "Failed triply degenerate eigen vector decomposition: %s %s." %
                            (str(eigenStruct.eigenVectors), str(lam123)))
        return

    #---------------------------------------------------------------------------
    # eigenVectors (diagonal matrix input)
    #---------------------------------------------------------------------------
    def testDiagonalEigenVectors(self):
        for i in range(self.ntests):
            lam1 = random.uniform(-ranrange, ranrange)
            lam2 = random.uniform(-ranrange, ranrange)
            lam3 = random.uniform(-ranrange, ranrange)
            A = SymTensor3d(lam1, 0.0, 0.0,
                            0.0, lam2, 0.0,
                            0.0, 0.0, lam3)
            lam0 = [(lam1, Vector3d(1, 0, 0)),
                    (lam2, Vector3d(0, 1, 0)),
                    (lam3, Vector3d(0, 0, 1))]
            lam0.sort()
            eigenVecs0 = [x[1] for x in lam0]
            eigenStruct = A.eigenVectors()
            lam = [(eigenStruct.eigenValues(i), eigenStruct.eigenVectors.getColumn(i)) for i in range(3)]
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
