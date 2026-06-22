#ATS:test(SELF,       label="sampleMultipleFields2Lattice 1D unit tests (serial)")
#ATS:test(SELF, np=4, label="sampleMultipleFields2Lattice 1D unit tests (4 proc)")

import unittest

from Spheral import *
import mpi

from testSampleMultipleFields2Lattice import TestSampleMultipleFields2Lattice

#===============================================================================
# 1-D tests.
#===============================================================================
class TestSampleMultipleFields2Lattice1d(TestSampleMultipleFields2Lattice,
                                         unittest.TestCase):

    #---------------------------------------------------------------------------
    # Initialize the problem.
    #---------------------------------------------------------------------------
    def setUp(self):

        from Spheral1d import (vector_of_int, Vector, Tensor, GammaLawGasMKS,
                               TableKernel, BSplineKernel, makeFluidNodeList,
                               ScalarField, VectorField, DataBase, Plane,
                               PeriodicBoundary)

        self.ndim = 1
        self.xmin = Vector(0.0)
        self.xmax = Vector(1.0)
        self.nsample = vector_of_int()
        self.nsample.append(100)

        # Tolerances for the test
        self.scalarTol = 1.0e-5
        self.vectorTol = 1.0e-3
        self.tensorTol = 1.0e-4

        n = 100
        self.rho0 = 10.0
        self.v0 = Vector(1.0)
        self.eps0 = -1.0
        self.gradv0 = Tensor(8.0)
        x0, x1 = 0.0, 1.0
        
        # Create the nodes and such.
        self.eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.WT = TableKernel(BSplineKernel())
        self.nodes = makeFluidNodeList("nodes", self.eos)

        # Distribute the nodes.
        from DistributeNodes import distributeNodesInRange1d
        distributeNodesInRange1d([(self.nodes, n, self.rho0, (x0, x1))])

        # Set the velocities and energies.
        self.nodes.velocity(VectorField("tmp", self.nodes, self.v0))
        self.nodes.specificThermalEnergy(ScalarField("tmp", self.nodes, self.eps0))

        self.db = DataBase()
        self.db.appendNodeList(self.nodes)

        # Create the boundary conditions.
        p0 = Plane(Vector(0.0), Vector(1.0))
        p1 = Plane(Vector(1.0), Vector(-1.0))
        xbc = PeriodicBoundary(p0, p1)
        self.bcs = [xbc]
        try:
            self.bcs.append(TreeDistributedBoundary1d.instance())
        except:
            if mpi.procs > 1:
                raise RuntimeError("Unable to get parallel boundary condition")
            else:
                pass

        # Enforce boundaries.
        db = DataBase()
        db.appendNodeList(self.nodes)
        for bc in self.bcs:
            bc.setAllGhostNodes(db)
            bc.finalizeGhostBoundary()
            self.nodes.neighbor().updateNodes()
        for bc in self.bcs:
            bc.applyGhostBoundary(self.nodes.mass())
            bc.applyGhostBoundary(self.nodes.massDensity())
            bc.applyGhostBoundary(self.nodes.specificThermalEnergy())
            bc.applyGhostBoundary(self.nodes.velocity())
        for bc in self.bcs:
            bc.finalizeGhostBoundary()

        self.H0 = self.nodes.Hfield()[0]
        return

    def tearDown(self):
        del self.nodes

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
