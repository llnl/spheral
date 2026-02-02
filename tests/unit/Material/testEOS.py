#ATS:test(SELF, label="EquationOfState unit tests")

from SpheralTestUtilities import fuzzyEqual
import numpy as np
from math import *
import unittest

from Spheral1d import *

# Create a global random number generator.
import random
random.seed(5479029084)

# We'll just use a gamma-law gas to base our tests on.
gamma = 5.0/3.0
mu = 5.0
units = CGS()
eos = GammaLawGas(gamma = gamma,
                  mu = mu,
                  constants = units)
gamma1 = gamma - 1.0

#-------------------------------------------------------------------------------
# Test harness.
#-------------------------------------------------------------------------------
class EOSTest(unittest.TestCase):

    def testSpecificThermalEnergyLookup(self):
        n = 1000000                      # How many tests?
        rhoMin, rhoMax = 1e-5, 1000.0    # range of densities
        Pmin, Pmax = 0.0, 1e5            # range of target pressures
        epsTol, Ptol = 1e-10, 1e-10      # tolerances
        Perrcheck = 1e-4                 # error check
        maxIterations = 100
        rho0 = np.random.uniform(rhoMin, rhoMax, n)
        P0 = np.random.uniform(Pmin, Pmax, n)

        # for (rho0i, P0i) in zip(rho0, P0):
        #     epsi = eos.specificThermalEnergyForPressure(P0i, rho0i, 0.0, 1e10, epsTol, Ptol, maxIterations, verbose=False)
        #     Pi = eos.pressure(rho0i, epsi)
        #     Perri = abs(Pi - P0i)/max(1e-10, P0i)
        #     self.assertTrue(Perri < Perrcheck,
        #                     f"Pressure error out of tolerance: {Perri} > {Perrcheck}")

        # Now some cuteness to call a method across elements of numpy arrays.
        epsLookup = np.vectorize(eos.specificThermalEnergyForPressure)
        Plookup = np.vectorize(eos.pressure)

        # Find the estimated eps(P, rho)
        eps = epsLookup(P0, rho0, 0.0, 1e10, epsTol, Ptol, maxIterations)
        
        # The corresponding pressures.
        P = Plookup(rho0, eps)

        Perr = np.minimum(P0, np.abs((P - P0)/np.maximum(1.0e-10, P0)))
        self.assertTrue((Perr < Perrcheck).all(),
                        "Pressure error out of tolerance: %s > %s" % (Perr.max(), Perrcheck))
        return

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
