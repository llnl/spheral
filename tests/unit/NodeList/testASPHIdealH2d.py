from math import *
import unittest

import Generate2dTestSetup

from Spheral import *
import Gnuplot
from PlotKernel2d import *
from SpheralGnuPlotUtilities import *

def W(eta):
    if abs(eta) < 2.0:
        return 2.0 - abs(eta)
    else:
        return 0.0

def findMoment(nodes, inode):
    #W = HatKernel2d(2.0, 2.0) # HKernel2d()
    nodes.neighbor().setMasterList(inode)
    nodes.neighbor().setRefineNeighborList(inode)
    r = nodes.positions()
    H = nodes.Hfield()
    ri = r[inode]
    Hi = H[inode]

    psi = SymTensor2d()
    for j in nodes.neighbor().refineNeighborList():
        rij = r[j] - ri
        eta = Hi*rij
        psi += eta.selfdyad()/(eta.magnitude()**4 + 1.0e-10)*W(eta.magnitude())

    psi /= sqrt(psi.Determinant())
    psii = psi.Inverse()

    eigen = psii.eigenVectors()
    result = SymTensor2d(eigen.eigenValues.x, 0.0,
                         0.0, eigen.eigenValues.y)
    result.rotationalTransform(eigen.eigenVectors)

    return result

#===============================================================================
# Main testing class.
#===============================================================================
class TestASPHIdealH(unittest.TestCase):

    #===========================================================================
    # Setup, create a uniform node distribution.
    #===========================================================================
    def setUp(self):
        self.ratiofuzz = 0.15
        self.volfuzz = 0.1
        self.aligncutoff = 0.8
        self.alignfuzz = 0.1
        self.iterations = 50
        self.ntests = 1

##        self.nx = 101
##        self.inode = 50*101 + 50

        self.nx = 21
        halfnx = self.nx // 2
        self.inode = halfnx*self.nx + halfnx
        self.inode0 = 0
        self.inode1 = halfnx*self.nx

        self.nPerh = 2.01
        self.testdata = Generate2dTestSetup.Generate2dTestSetup(asph = True,
                                                                seed = 'lattice',
                                                                nodesPerh = self.nPerh,
                                                                n1 = self.nx*self.nx,
                                                                n2 = 0,
                                                                n3 = 0)
        self.testdata.nodes1.HsmoothFraction = 0.01
        self.testdata.nodes1.nodesPerSmoothingScale = self.nPerh
        self.q = MonaghanGingoldViscosity2d(1.0, 1.0)
        self.WT = TableKernel2d(BSplineKernel2d())
        self.hydro = Hydro2d(self.WT, self.WT, self.q)
        self.hydro.HsmoothMin = 1.0e-5
        self.hydro.HsmoothMax = 1.0e5
        self.hydro.HEvolution = Hydro2d.HEvolutionType.IdealH

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
        Hi = SymTensor2d(self.testdata.nodes1.Hfield()[self.inode])
        assert ri == Vector2d(0.5, 0.5)

        self.testdata.nodes1.Hfraction = 0.5

        plot = Gnuplot.Gnuplot(persist = True)

        for i in range(self.ntests):

            # Pick a random multiplicative distortion to apply to the node
            # position field.
            mx = self.testdata.g.uniform(0.5, 2.0)
            my = self.testdata.g.uniform(0.5, 2.0)
            theta = self.testdata.g.uniform(0.0, 2.0*pi)
            mx, my, theta = 0.5, 2.0, 0.*pi
            T0 = SymTensor2d(mx, 0.0,
                             0.0, my)
            R = Tensor2d(cos(theta), sin(theta),
                         -sin(theta), cos(theta))
            T = R*T0
            Rinverse = R.Inverse()
            Tinverse = T.Inverse()
            print("mx, my, theta/pi: ", mx, my, theta/pi)
            print("T, Ti: ", T, Tinverse)

            # Compute the expected answer for the ideal H.
            answer = (R*(T0.Inverse()*Hi).Symmetric()*R.Inverse()).Symmetric()

            # Distort the node positions.
            self.distortNodeDistribution(T)

            # Run the ideal H algorithm by evaluating the fluid derivatives.
            H = self.testdata.dataBase.globalHfield
            for iter in range(self.iterations):
                self.testdata.dataBase.updateConnectivityMap()
                packages = vector_of_Physics2d_ptr()
                packages.append(self.hydro)
                state = State2d(self.testdata.dataBase, packages)
                derivs = StateDerivatives2d(self.testdata.dataBase, packages)
                self.hydro.evaluateDerivatives(1.0, 1.0, self.testdata.dataBase, state, derivs)
                Hnew = self.hydro.Hideal()
                print("(%s %s) : (%s %s %s)" % (str(Hi), str(answer),
                                                str(Hnew[0][self.inode]),
                                                str(Hnew[0][self.inode0]),
                                                str(Hnew[0][self.inode1])))
##                self.testdata.nodes1.Hfield()[self.inode] = Hip
                H.assignFields(Hnew)
                self.testdata.nodes1.neighbor().updateNodes()

            rt = T*ri
##            if plot:
##                del plot
##            plot = plotKernel2d(self.testdata.nodes1,
##                                xmin = ri.x - 0.1,
##                                xmax = ri.x + 0.1,
##                                ymin = ri.y - 0.1,
##                                ymax = ri.y + 0.1,
##                                h = 0.2,
##                                np = 10)

            Hinverse = self.testdata.dataBase.fluidHinverse
            from SpheralVisitDump import SpheralVisitDump
            dumper = SpheralVisitDump("Htest",
                                      listOfFieldLists = [Hinverse]
                                      )
            dumper.dump(0.0, 0)

##            print "Second moments:"
##            for mom in derivs.massSecondMoment()[0].internalValues():
##                print mom

            # Test the measured ideal H against our answer.
            eigenanswer = answer.eigenVectors()
            if eigenanswer.eigenValues.x < eigenanswer.eigenValues.y:
                h1, v1 = eigenanswer.eigenValues.x, eigenanswer.eigenVectors.getColumn(0)
                h2, v2 = eigenanswer.eigenValues.y, eigenanswer.eigenVectors.getColumn(1)
            else:
                h2, v2 = eigenanswer.eigenValues.x, eigenanswer.eigenVectors.getColumn(0)
                h1, v1 = eigenanswer.eigenValues.y, eigenanswer.eigenVectors.getColumn(1)

            Hip = SymTensor2d(Hnew[0][self.inode])
            eigenp = Hip.eigenVectors()
            if eigenp.eigenValues.x < eigenp.eigenValues.y:
                h1p, v1p = eigenp.eigenValues.x, eigenp.eigenVectors.getColumn(0)
                h2p, v2p = eigenp.eigenValues.y, eigenp.eigenVectors.getColumn(1)
            else:
                h2p, v2p = eigenp.eigenValues.x, eigenp.eigenVectors.getColumn(0)
                h1p, v1p = eigenp.eigenValues.y, eigenp.eigenVectors.getColumn(1)

            self.assertTrue(fuzzyEqual(h1p/h2p, h1/h2, self.ratiofuzz),
                            "H Shape does not match expectation: (%s/%s = %s) != (%s/%s = %s)" %
                            (h1p, h2p, h1p/h2p, h1, h2, h1/h2))
            self.assertTrue(fuzzyEqual(1.0/sqrt(Hip.Determinant()), 1.0/sqrt(answer.Determinant()), self.volfuzz),
                            "H Determinant does not match expectation: %s != %s" %
                            (Hip.Determinant(), answer.Determinant()))
            self.assertTrue(h1/h2 > self.aligncutoff or
                            (fuzzyEqual(abs(v1.dot(v1p)), 1.0, self.alignfuzz) and
                             fuzzyEqual(abs(v2.dot(v2p)), 1.0, self.alignfuzz)),
                            "H not alinged with expected directions: %s %s %s %s != %s %s %s %s" %
                            (h1p, v1p, h2p, v2p, h1, v1, h2, v2))

            # Restore the node positions.
            self.distortNodeDistribution(Tinverse)
            assert self.testdata.nodes1.positions()[self.inode] == ri

        return

##    def testBiteme(self):
##        theta = pi/4.0 # self.testdata.g.uniform(0.0, 2.0*pi)
##        mx = 1.0 # self.testdata.g.uniform(0.1, 2.0)
##        my = 0.1 # self.testdata.g.uniform(0.1, 2.0)
##        T = SymTensor2d(mx, 0.0,
##                        0.0, my)
##        R = Tensor2d(cos(theta), -sin(theta),
##                     sin(theta), cos(theta))
##        T.rotationalTransform(R)
##        self.distortNodeDistribution(T)
##        p = plotNodePositions2d(self.testdata.dataBase,
##                                persist = True)

##        result = findMoment(self.testdata.nodes1, self.inode)
##        print "Moment @ %s: " % self.testdata.nodes1.positions()[self.inode], result, result.eigenVectors().eigenValues

if __name__ == "__main__":
    unittest.main()

