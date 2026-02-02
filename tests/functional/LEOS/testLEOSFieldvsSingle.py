#ATS:test(SELF, label="LEOS Field vs single value tests")
#-------------------------------------------------------------------------------
# LEOS unit test: compare looking up values in the Field and single value
# interfaces
#-------------------------------------------------------------------------------
from Spheral1d import *
import numpy as np
import unittest

import random
random.seed(45577959828927)

# LEOS is not always built, so if we don't have it just pass
import Spheral
if not hasattr(Spheral, "LEOS1d"):
    exit(0)

# We'll work in CGuS units.
units = PhysicalConstants(0.01,     # Unit length in meters
                          0.001,    # Unit mass in kg
                          1.0e-6)   # Unit time in sec

#-------------------------------------------------------------------------------
# Build a LEOS equation of state for SiO2
#-------------------------------------------------------------------------------
eosLEOS = LEOS(materialNumber = 2210,
               constants = units)
#print(eosLEOS.descriptor)
rho0 = eosLEOS.referenceDensity
T0 = eosLEOS.referenceTemperature
eps0 = eosLEOS.specificThermalEnergy(rho0, T0)

# The ranges for random values
rhoMin, rhoMax = 1e-3*rho0, 1e2*rho0
Tmin, Tmax = 1.0, 1.0e8
epsMin, epsMax = eosLEOS.specificThermalEnergy(rho0, Tmin), eosLEOS.specificThermalEnergy(rho0, Tmax)

nodes = makeFluidNodeList(name = "test nodes",
                          eos = eosLEOS,
                          numInternal = 1000)

# Generate a random set of (rho, eps, T) test values
rhoField = ScalarField("test mass density", nodes)
epsField = ScalarField("test specific thermal energy", nodes)
Tfield = ScalarField("test temperature", nodes)
for i in range(nodes.numInternalNodes):
    rhoField[i] = random.uniform(rhoMin, rhoMax)
    epsField[i] = random.uniform(epsMin, epsMax)
    Tfield[i] = random.uniform(Tmin, Tmax)

#-------------------------------------------------------------------------------
# The test class
#
# For now we just check that the Field methods return the same value as the
# individual value lookup.
#-------------------------------------------------------------------------------
class TestLEOSFieldVsSingleValue(unittest.TestCase):

    #---------------------------------------------------------------------------
    # Generic test F(x,y) without derivatives
    #---------------------------------------------------------------------------
    def genericLookupTest(self,
                          label,
                          fieldMethod,
                          valueMethod,
                          x,
                          y):
        F = ScalarField(label, nodes)
        fieldMethod(F, x, y)
        for i in range(nodes.numInternalNodes):
            self.assertTrue(np.isclose(F[i], valueMethod(x[i], y[i])),
                            "Failed {} equality @ ({}, {}): {} != {}".format(label,
                                                                             x[i],
                                                                             y[i],
                                                                             F[i],
                                                                             valueMethod(x[i], y[i])))
        return

    #---------------------------------------------------------------------------
    # Generic test F(x,y) with derivatives
    #---------------------------------------------------------------------------
    def genericLookupWithDerivsTest(self,
                                    label,
                                    fieldMethod,
                                    valueMethod,
                                    x,
                                    y):
        F = ScalarField(label, nodes)
        dFde = ScalarField(label + " dFde", nodes)
        dFdr = ScalarField(label + " dPdr", nodes)
        fieldMethod(F, dFde, dFdr, x, y)
        for i in range(nodes.numInternalNodes):
            F_a, dFde_a, dFdr_a = valueMethod(x[i], y[i])
            self.assertTrue(np.isclose(F[i], F_a),
                            "Failed {} equality @ ({}, {}): {} != {}".format(label,
                                                                             x[i],
                                                                             y[i],
                                                                             F[i],
                                                                             F_a))
            self.assertTrue(np.isclose(dFde[i], dFde_a),
                            "Failed {} dFde equality @ ({}, {}): {} != {}".format(label,
                                                                                  x[i],
                                                                                  y[i],
                                                                                  dFde[i],
                                                                                  dFde_a))
            self.assertTrue(np.isclose(dFdr[i], dFdr_a),
                            "Failed {} dPdr equality @ ({}, {}): {} != {}".format(label,
                                                                                  x[i],
                                                                                  y[i],
                                                                                  dFdr[i],
                                                                                  dFdr_a))
        return

    #---------------------------------------------------------------------------
    # Pressure
    #---------------------------------------------------------------------------
    def testPressureField(self):
        self.genericLookupTest("pressure", eosLEOS.setPressure, eosLEOS.pressure, rhoField, epsField)
        return

    #---------------------------------------------------------------------------
    # Pressure and derivs
    #---------------------------------------------------------------------------
    def testPressureAndDerivsField(self):
        self.genericLookupWithDerivsTest("pressure", eosLEOS.setPressureAndDerivs, eosLEOS.pressureAndDerivs, rhoField, epsField)
        return

    #---------------------------------------------------------------------------
    # Temperature
    #---------------------------------------------------------------------------
    def testTemperatureField(self):
        self.genericLookupTest("temperature", eosLEOS.setTemperature, eosLEOS.temperature, rhoField, Tfield)
        return

    #---------------------------------------------------------------------------
    # Specific thermal energy
    #---------------------------------------------------------------------------
    def testSpecificThermalEnergyField(self):
        self.genericLookupTest("specific thermal energy", eosLEOS.setSpecificThermalEnergy, eosLEOS.specificThermalEnergy, rhoField, Tfield)
        return

    #---------------------------------------------------------------------------
    # Specific heat
    #---------------------------------------------------------------------------
    def testSpecificHeatField(self):
        self.genericLookupTest("specific heat", eosLEOS.setSpecificHeat, eosLEOS.specificHeat, rhoField, Tfield)
        return

    #---------------------------------------------------------------------------
    # Sound speed
    #---------------------------------------------------------------------------
    def testSoundSpeedField(self):
        self.genericLookupTest("sound speed", eosLEOS.setSoundSpeed, eosLEOS.soundSpeed, rhoField, epsField)
        return

    #---------------------------------------------------------------------------
    # Gamma
    #---------------------------------------------------------------------------
    def testGammaField(self):
        self.genericLookupTest("gamma", eosLEOS.setGammaField, eosLEOS.gamma, rhoField, epsField)
        return

    #---------------------------------------------------------------------------
    # bulk modulus
    #---------------------------------------------------------------------------
    def testBulkModulusField(self):
        self.genericLookupTest("bulk modulus", eosLEOS.setBulkModulus, eosLEOS.bulkModulus, rhoField, epsField)
        return

    #---------------------------------------------------------------------------
    # Entropy
    #---------------------------------------------------------------------------
    def testEntropyField(self):
        self.genericLookupTest("entropy", eosLEOS.setEntropy, eosLEOS.entropy, rhoField, epsField)
        return

    #---------------------------------------------------------------------------
    # Melt temperature
    #---------------------------------------------------------------------------
    def testMeltTemperatureField(self):
        self.genericLookupTest("melt temperature", eosLEOS.setMeltTemperature, eosLEOS.meltTemperature, rhoField, epsField)
        return

#-------------------------------------------------------------------------------
# Run all tests
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
