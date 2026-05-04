from math import *
import unittest

import Generate3dTestSetup

from Spheral import *
import Gnuplot
from SpheralGnuPlotUtilities import *

#===============================================================================
# Main testing class.
#===============================================================================
class TestASPHIdealH3d(unittest.TestCase):

    #===========================================================================
    # Setup, create a uniform node distribution.
    #===========================================================================
    def setUp(self):
        self.ratiofuzz = 0.15
        self.volfuzz = 0.1
        self.aligncutoff = 0.8
        self.alignfuzz = 0.1
        self.iterations = 10
        self.ntests = 1

        self.nx = 11
        halfnx = self.nx/2
        self.inode = halfnx*self.nx*self.nx + halfnx*self.nx + halfnx

        self.nPerh = 2.01
        self.testdata = Generate3dTestSetup.Generate3dTestSetup(asph = True,
                                                                seed = 'lattice',
                                                                nodesPerh = self.nPerh,
                                                                nx1 = self.nx,
                                                                nx2 = 0,
                                                                nx3 = 0)
        self.q = MonaghanGingoldViscosity3d()
        self.WT = TableKernel3d(BSplineKernel3d())
        self.hydro = Hydro3d(self.WT, self.WT, self.q)

        return

    #===========================================================================
    # Given a vector with desired x, y multiplication factors, apply the
    # desired mapping to the node positions.
    #===========================================================================
    def distortNodeDistribution(self, T):
        for i in range(self.testdata.nodes1.numNodes):
            self.testdata.nodes1.positions()[i] = T*self.testdata.nodes1.positions()[i]
        self.testdata.nodes1.neighbor().updateNodes()
        return

    #===========================================================================
    # The actual test.
    #===========================================================================
    def testIdealH(self):

        # The state of the node we're testing.
        ri = self.testdata.nodes1.positions()[self.inode]
        Hi = SymTensor3d(self.testdata.nodes1.Hfield()[self.inode])
        print(ri, Hi)
        assert ri == Vector3d(0.5, 0.5, 0.5)

        for i in range(self.ntests):

            # Pick a random multiplicative distortion to apply to the node
            # position field.
            mx = self.testdata.g.uniform(0.5, 2.0)
            my = self.testdata.g.uniform(0.5, 2.0)
            mz = self.testdata.g.uniform(0.5, 2.0)
            theta = self.testdata.g.uniform(0.0, 2.0*pi)
            phi = self.testdata.g.uniform(0.0, pi)
            mx, my, mz, theta, phi = 0.5, 2.0, 1.0, 0.0*pi, 0.0
            T0 = SymTensor3d(mx, 0.0, 0.0,
                             0.0, my, 0.0,
                             0.0, 0.0, mz)
            Rtheta = Tensor3d(cos(theta), sin(theta), 0.0,
                              -sin(theta), cos(theta), 0.0,
                              0.0, 0.0, 1.0)
            Rphi = Tensor3d(cos(phi), 0.0, sin(phi),
                            0.0,      1.0, 0.0,
                            -sin(phi), 0.0, cos(phi))
            R = Rphi*Rtheta
            T = R*T0
            Rinverse = R.Inverse()
            Tinverse = T.Inverse()
            print("mx, my, mz, theta/pi, phi/pi: ", mx, my, mz, theta/pi, phi/pi)
            print("T, Ti: ", T, Tinverse)

            # Compute the expected answer for the ideal H.
            answer = (R*(T0.Inverse()*Hi).Symmetric()*R.Inverse()).Symmetric()

            # Distort the node positions.
            self.distortNodeDistribution(T)

            # Run the ideal H algorithm by evaluating the fluid derivatives.
            for iter in range(self.iterations):
                print("Iteration ", iter)
                derivs = StateDerivatives3d(self.testdata.dataBase)
                derivs.Zero()
                self.hydro.evaluateDerivatives(1.0, 1.0, self.testdata.dataBase, derivs)
                print("Hooo boy")
                Hip = SymTensor3d(derivs.Hideal()[0][self.inode])
                print("Original, Answer, Hip: ", Hi, answer, Hip)
                self.testdata.nodes1.Hfield()[self.inode] = Hip
                self.testdata.nodes1.neighbor().updateNodes()

##            print "Plotting..."
##            rt = T*ri
##            plot = plotKernel3d(self.testdata.nodes1,
##                                xmin = ri.x - 0.1,
##                                xmax = ri.x + 0.1,
##                                ymin = ri.y - 0.1,
##                                ymax = ri.y + 0.1,
##                                h = 0.2,
##                                plot = plot)

##            from SpheralVisitDump import SpheralVisitDump
##            dumper = SpheralVisitDump(self.testdata.dataBase,
##                                      "Htest",
##                                      listOfFields = [self.testdata.nodes1.Hfield()]
##                                      )
##            dumper.dump(0.0, 0)


##            print "Second moments:"
##            for mom in derivs.massSecondMoment()[0].internalValues():
##                print mom

##            # Test the measured ideal H against our answer.
##            eigenanswer = answer.eigenVectors()
##            if eigenanswer.eigenValues.x < eigenanswer.eigenValues.y:
##                h1, v1 = eigenanswer.eigenValues.x, eigenanswer.eigenVectors.getColumn(0)
##                h2, v2 = eigenanswer.eigenValues.y, eigenanswer.eigenVectors.getColumn(1)
##            else:
##                h2, v2 = eigenanswer.eigenValues.x, eigenanswer.eigenVectors.getColumn(0)
##                h1, v1 = eigenanswer.eigenValues.y, eigenanswer.eigenVectors.getColumn(1)

##            eigenp = Hip.eigenVectors()
##            if eigenp.eigenValues.x < eigenp.eigenValues.y:
##                h1p, v1p = eigenp.eigenValues.x, eigenp.eigenVectors.getColumn(0)
##                h2p, v2p = eigenp.eigenValues.y, eigenp.eigenVectors.getColumn(1)
##            else:
##                h2p, v2p = eigenp.eigenValues.x, eigenp.eigenVectors.getColumn(0)
##                h1p, v1p = eigenp.eigenValues.y, eigenp.eigenVectors.getColumn(1)

##            self.assertTrue(fuzzyEqual(h1p/h2p, h1/h2, self.ratiofuzz),
##                            "H Shape does not match expectation: (%s/%s = %s) != (%s/%s = %s)" %
##                            (h1p, h2p, h1p/h2p, h1, h2, h1/h2))
##            self.assertTrue(fuzzyEqual(1.0/sqrt(Hip.Determinant()), 1.0/sqrt(answer.Determinant()), self.volfuzz),
##                            "H Determinant does not match expectation: %s != %s" %
##                            (Hip.Determinant(), answer.Determinant()))
##            self.assertTrue(h1/h2 > self.aligncutoff or
##                            (fuzzyEqual(abs(v1.dot(v1p)), 1.0, self.alignfuzz) and
##                             fuzzyEqual(abs(v2.dot(v2p)), 1.0, self.alignfuzz)),
##                            "H not alinged with expected directions: %s %s %s %s != %s %s %s %s" %
##                            (h1p, v1p, h2p, v2p, h1, v1, h2, v2))

            # Restore the node positions.
            self.distortNodeDistribution(Tinverse)
            assert self.testdata.nodes1.positions()[self.inode] == ri

        return

if __name__ == "__main__":
    unittest.main()

