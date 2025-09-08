#ATS:test(SELF, label="LEOS analytic EOS comparison (Tillotson and gamma law gas)", nt=10)
#-------------------------------------------------------------------------------
# LEOS unit test: generate an LEOS input from an analytic EOS, and compare
# looking up values between the two
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import fuzzyEqual
import numpy as np
import unittest
from writeLEOSfile import *

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
# The test class (generic harness)
#
# Test calls to the LEOS interpolated version of a reference eos match the
# reference values
#-------------------------------------------------------------------------------
class TestLEOSvsReference:

    #---------------------------------------------------------------------------
    # Test an individual array
    #---------------------------------------------------------------------------
    def fieldTest(self,
                  label,
                  x,
                  y,
                  refField,
                  testField,
                  rtol,
                  atol):
        for i in range(self.nodes.numInternalNodes):
            self.assertTrue(np.allclose(testField[i], refField[i], rtol, atol),
                            "Failed {} equality @ ({}, {}): {} != {}".format(label,
                                                                             x[i],
                                                                             y[i],
                                                                             testField[i],
                                                                             refField[i]))
        return

    #---------------------------------------------------------------------------
    # Generic test F(x,y) without derivatives
    #---------------------------------------------------------------------------
    def genericLookupTest(self,
                          label,
                          refMethod,
                          leosMethod,
                          x,
                          y,
                          rtol = 1.0e-3,
                          atol = 1.0e-8):
        Fr = ScalarField(label + " ref", self.nodes)
        Fl = ScalarField(label + " leos", self.nodes)
        refMethod(Fr, x, y)
        leosMethod(Fl, x, y)
        self.fieldTest(label, x, y, Fr, Fl, rtol, atol)
        return

    #---------------------------------------------------------------------------
    # Generic test F(x,y) with derivatives
    #---------------------------------------------------------------------------
    def genericLookupWithDerivsTest(self,
                                    label,
                                    refMethod,
                                    leosMethod,
                                    x,
                                    y,
                                    rtol = 1.0e-3,
                                    atol = 1.0e-3,
                                    dxrtol = 1.0e-3,
                                    dxatol = 1.0e-3,
                                    dyrtol = 1.0e-3,
                                    dyatol = 1.0e-3):
        Fr = ScalarField(label + " ref", self.nodes)
        dFdx_r = ScalarField(label + " ref dFdx", self.nodes)
        dFdy_r = ScalarField(label + " ref dFdy", self.nodes)
        Fl = ScalarField(label + " leos", self.nodes)
        dFdx_l = ScalarField(label + " leos dFdx", self.nodes)
        dFdy_l = ScalarField(label + " leos dFdy", self.nodes)
        refMethod(Fr, dFdy_r, dFdx_r, x, y)
        leosMethod(Fl, dFdy_l, dFdx_l, x, y)
        self.fieldTest(label, x, y, Fr, Fl, rtol, atol)
        self.fieldTest("dFdx " + label, x, y, dFdx_r, dFdx_l, dxrtol, dxatol)
        self.fieldTest("dFdy " + label, x, y, dFdy_r, dFdy_l, dyrtol, dyatol)
        return

    #---------------------------------------------------------------------------
    # Pressure
    #---------------------------------------------------------------------------
    def testPressureField(self):
        self.genericLookupTest("pressure", self.eosRef.setPressure, self.eosLEOS.setPressure, self.rhoField, self.epsField, self.Prtol, self.Patol)
        return

    #---------------------------------------------------------------------------
    # Pressure and derivs
    #---------------------------------------------------------------------------
    def testPressureAndDerivsField(self):
        self.genericLookupWithDerivsTest("pressure", self.eosRef.setPressureAndDerivs, self.eosLEOS.setPressureAndDerivs, self.rhoField, self.epsField,
                                         self.Prtol, self.Patol, self.Pdxrtol, self.Pdxatol, self.Pdyrtol, self.Pdyatol)
        return

    #---------------------------------------------------------------------------
    # Temperature
    #---------------------------------------------------------------------------
    def testTemperatureField(self):
        self.genericLookupTest("temperature", self.eosRef.setTemperature, self.eosLEOS.setTemperature, self.rhoField, self.Tfield, self.Trtol, self.Tatol)
        return

    #---------------------------------------------------------------------------
    # Specific thermal energy
    #---------------------------------------------------------------------------
    def testSpecificThermalEnergyField(self):
        self.genericLookupTest("specific thermal energy", self.eosRef.setSpecificThermalEnergy, self.eosLEOS.setSpecificThermalEnergy, self.rhoField, self.Tfield,
                               self.epsrtol, self.epsatol)
        return

    #---------------------------------------------------------------------------
    # Specific heat
    #---------------------------------------------------------------------------
    def testSpecificHeatField(self):
        self.genericLookupTest("specific heat", self.eosRef.setSpecificHeat, self.eosLEOS.setSpecificHeat, self.rhoField, self.Tfield, self.CVrtol, self.CVatol)
        return

    #---------------------------------------------------------------------------
    # Sound speed
    #---------------------------------------------------------------------------
    def testSoundSpeedField(self):
        self.genericLookupTest("sound speed", self.eosRef.setSoundSpeed, self.eosLEOS.setSoundSpeed, self.rhoField, self.epsField, self.CSrtol, self.CSatol)
        return

    #---------------------------------------------------------------------------
    # Gamma
    #---------------------------------------------------------------------------
    def testGammaField(self):
        self.genericLookupTest("gamma", self.eosRef.setGammaField, self.eosLEOS.setGammaField, self.rhoField, self.epsField, self.gamrtol, self.gamatol)
        return

    #---------------------------------------------------------------------------
    # bulk modulus
    #---------------------------------------------------------------------------
    def testBulkModulusField(self):
        self.genericLookupTest("bulk modulus", self.eosRef.setBulkModulus, self.eosLEOS.setBulkModulus, self.rhoField, self.epsField, self.Krtol, self.Katol)
        return

    #---------------------------------------------------------------------------
    # Entropy
    #---------------------------------------------------------------------------
    def testEntropyField(self):
        self.genericLookupTest("entropy", self.eosRef.setEntropy, self.eosLEOS.setEntropy, self.rhoField, self.epsField, self.Srtol, self.Satol)
        return

#     #---------------------------------------------------------------------------
#     # Melt temperature
#     #---------------------------------------------------------------------------
#     def testMeltTemperatureField(self):
#         self.genericLookupTest("melt temperature", self.eosRef.setMeltTemperature, self.eosLEOS.meltTemperature, self.rhoField, self.epsField)
#         return

#-------------------------------------------------------------------------------
# Tillotson reference (basalt)
#-------------------------------------------------------------------------------
class TestLEOSvsTillotson(TestLEOSvsReference,
                          unittest.TestCase):

    #---------------------------------------------------------------------------
    # Make the EOS static members
    #---------------------------------------------------------------------------
    # Build a Tillotson version of basalt
    etaMin, etaMax = 0.2, 5.0
    etaMinTest, etaMaxTest = 0.5, 2.0
    Tmin, Tmax = 1.0e0, 1.0e8
    TminTest, TmaxTest = 50.0, 1.0e6
    eosRef = TillotsonEquationOfState(materialName = "basalt",
                                      etamin = etaMin,
                                      etamax = etaMax,
                                      units = units)

    # Write an LEOS input file sampling the Tillotson
    rho0 = eosRef.referenceDensity
    rhoMin = etaMin*rho0
    rhoMax = etaMax*rho0
    rhoMinTest = etaMinTest*rho0
    rhoMaxTest = etaMaxTest*rho0
    writeLEOSfile(eos = eosRef,
                  rhoMin = rhoMin,
                  rhoMax = rhoMax,
                  Tmin = Tmin,
                  Tmax = Tmax,
                  nrho = 1000,
                  nT = 1000,
                  basename = "eosTillotsonBasalt_testLEOS",
                  name = "Basalt",
                  eosNumber = 999,
                  T0 = 300.0,
                  atomicWeight = [28.086,  # Si,
                                  15.99,   # O
                                  26.982,  # Al
                                  55.854,  # Fe
                                  40.078,  # Ca
                                  24.305], # Mg
                  atomicNumber = [14.0,    # Si
                                  8.0,     # O
                                  13.0,    # Al
                                  26.0,    # Fe
                                  20.0,    # Ca
                                  12.0],   # Mg
                  atomicFrac =   [0.3333333333333,
                                  0.6666666666666,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0])

    # Read the generated file into an LEOS equation of state
    eosLEOS = LEOS(materialNumber = 999,
                   constants = units,
                   dbname = "eosTillotsonBasalt_testLEOS.ascii",
                   leosFileFormat = "ascii")

    # Make some nodes to use making test lookups
    nodes = makeFluidNodeList(name = "test nodes Tillotson",
                                   eos = eosLEOS,
                                   numInternal = 1000)

    # Generate a random set of (rho, eps, T) test values
    epsMinTest, epsMaxTest = eosLEOS.specificThermalEnergy(rho0, TminTest), eosLEOS.specificThermalEnergy(rho0, TmaxTest)
    rhoField = ScalarField("test mass density", nodes)
    epsField = ScalarField("test specific thermal energy", nodes)
    Tfield = ScalarField("test temperature", nodes)
    for i in range(nodes.numInternalNodes):
        rhoField[i] = random.uniform(rhoMinTest, rhoMaxTest)
        epsField[i] = random.uniform(epsMinTest, epsMaxTest)
        Tfield[i] = random.uniform(TminTest, TmaxTest)
    
    # Define our tolerances
    Prtol,   Patol   = 1.0e-2, 5.0e-1
    Pdxrtol, Pdxatol = 1.0e-0, 1.0e-1
    Pdyrtol, Pdyatol = 1.0e-0, 2.0e-0
    Trtol,   Tatol   = 1.0e-0, 1.0
    epsrtol, epsatol = 1.0e-3, 1.0e-2
    CVrtol,  CVatol  = 1.0e-3, 1.0e-2
    CSrtol,  CSatol  = 8.0e-1, 1.0e-0
    gamrtol, gamatol = 1.0e-0, 1.0e0
    Krtol,   Katol   = 1.0e-0, 1.0e-2
    Srtol,   Satol   = 1.0e-1, 1.0e-0

#-------------------------------------------------------------------------------
# Gamma-law gas reference
#
# We use a Tillotson with coefficients to give us gamma-law gas behavior
#-------------------------------------------------------------------------------
class TestLEOSvsGammaLaw(TestLEOSvsReference,
                         unittest.TestCase):

    #---------------------------------------------------------------------------
    # Make the EOS static members
    #---------------------------------------------------------------------------
    # Build a gamma-law gas analogous to Air
    etaMin, etaMax = 1e-3, 1e3
    etaMinTest, etaMaxTest = 0.1, 100.0
    Tmin, Tmax = 1.0e0, 1.0e8
    TminTest, TmaxTest = 50.0, 1.0e6
    eosRef = TillotsonEquationOfState(referenceDensity = 1.293e-3,  # g/cc
                                      etamin = etaMin,
                                      etamax = etaMax,
                                      etamin_solid = etaMin,
                                      etamax_solid = etaMax,
                                      a = 1.4,
                                      b = 0.0,
                                      A = 0.0,
                                      B = 0.0,
                                      alpha = 0.0,
                                      beta = 0.0,
                                      eps0 = 0.0,
                                      epsLiquid = -1e100,
                                      epsVapor = -1e100,
                                      atomicWeight = 28.96,
                                      constants = units)

    # Write an LEOS input file sampling the Tillotson
    rho0 = eosRef.referenceDensity
    rhoMin = etaMin*rho0
    rhoMax = etaMax*rho0
    rhoMinTest = etaMinTest*rho0
    rhoMaxTest = etaMaxTest*rho0
    writeLEOSfile(eos = eosRef,
                  rhoMin = rhoMin,
                  rhoMax = rhoMax,
                  Tmin = Tmin,
                  Tmax = Tmax,
                  nrho = 1000,
                  nT = 1000,
                  basename = "eosGammaLawAir_testLEOS",
                  name = "Air",
                  eosNumber = 999,
                  T0 = 300.0,
                  atomicWeight = [14.007,  # N
                                  15.99,   # O
                                  39.948], # Ar
                  atomicNumber = [7.0,     # N
                                  8.0,     # O
                                  18.0],   # Ar
                  atomicFrac =   [0.78,    # N
                                  0.21,    # O
                                  0.01])   # Ar

    # Read the generated file into an LEOS equation of state
    eosLEOS = LEOS(materialNumber = 999,
                   constants = units,
                   dbname = "eosGammaLawAir_testLEOS.ascii",
                   leosFileFormat = "ascii")

    # Make some nodes to use making test lookups
    nodes = makeFluidNodeList(name = "test nodes Air",
                                   eos = eosLEOS,
                                   numInternal = 1000)

    # Generate a random set of (rho, eps, T) test values
    epsMinTest, epsMaxTest = eosLEOS.specificThermalEnergy(rho0, TminTest), eosLEOS.specificThermalEnergy(rho0, TmaxTest)
    rhoField = ScalarField("test mass density", nodes)
    epsField = ScalarField("test specific thermal energy", nodes)
    Tfield = ScalarField("test temperature", nodes)
    for i in range(nodes.numInternalNodes):
        rhoField[i] = random.uniform(rhoMinTest, rhoMaxTest)
        epsField[i] = random.uniform(epsMinTest, epsMaxTest)
        Tfield[i] = random.uniform(TminTest, TmaxTest)
    
    # Define our tolerances
    Prtol,   Patol   = 1.0e-3, 1.0e-2
    Pdxrtol, Pdxatol = 1.0e-3, 1.0e-2
    Pdyrtol, Pdyatol = 1.0e-3, 1.0e-2
    Trtol,   Tatol   = 1.0e2, 1.0e3
    epsrtol, epsatol = 1.0e-3, 1.0e-2
    CVrtol,  CVatol  = 1.0e-3, 1.0e-2
    CSrtol,  CSatol  = 8.0e-1, 1.0e-1
    gamrtol, gamatol = 1.0e-0, 1.0e-2
    Krtol,   Katol   = 1.0e-0, 1.0e-2
    Srtol,   Satol   = 1.0e2, 1.0e2

# #-------------------------------------------------------------------------------
# # ANEOS SiO2 reference
# #-------------------------------------------------------------------------------
# class TestLEOSvsANEOS(TestLEOSvsReference,
#                       unittest.TestCase):

#     #---------------------------------------------------------------------------
#     # Make the EOS static members
#     #---------------------------------------------------------------------------
#     # Build an ANEOS version of SiO2
#     etaMin, etaMax = 0.2, 5.0
#     Tmin, Tmax = 1.0e0, 1.0e8
#     eosRef = ANEOS("quartz",
#                    numRhoVals = 800,
#                    numTvals = 800,
#                    Tmin = Tmin,
#                    Tmax = Tmax,
#                    constants = units)

#     # Write an LEOS input file sampling the Tillotson
#     rho0 = eosRef.referenceDensity
#     rhoMin = etaMin*rho0
#     rhoMax = etaMax*rho0
#     writeLEOSfile(eos = eosRef,
#                   rhoMin = rhoMin,
#                   rhoMax = rhoMax,
#                   Tmin = Tmin,
#                   Tmax = Tmax,
#                   nrho = 1000,
#                   nT = 1000,
#                   basename = "eosANEOSSiO2_testLEOS",
#                   name = "SiO2",
#                   eosNumber = 999,
#                   T0 = 300.0,
#                   atomicWeight = [28.086,  # Si,
#                                   15.99],  # O
#                   atomicNumber = [14.0,    # Si
#                                   8.0],    # O
#                   atomicFrac =   [0.3333333333333,
#                                   0.6666666666666])

#     # Read the generated file into an LEOS equation of state
#     eosLEOS = LEOS(materialNumber = 999,
#                    constants = units,
#                    dbname = "eosANEOSSiO2_testLEOS.ascii",
#                    leosFileFormat = "ascii")

#     # Make some nodes to use making test lookups
#     nodes = makeFluidNodeList(name = "test nodes ANEOS",
#                                    eos = eosLEOS,
#                                    numInternal = 1000)

#     # Generate a random set of (rho, eps, T) test values
#     epsMin, epsMax = eosLEOS.specificThermalEnergy(rho0, Tmin), eosLEOS.specificThermalEnergy(rho0, Tmax)
#     rhoField = ScalarField("test mass density", nodes)
#     epsField = ScalarField("test specific thermal energy", nodes)
#     Tfield = ScalarField("test temperature", nodes)
#     for i in range(nodes.numInternalNodes):
#         rhoField[i] = random.uniform(rhoMin, rhoMax)
#         epsField[i] = random.uniform(epsMin, epsMax)
#         Tfield[i] = random.uniform(Tmin, Tmax)
    
#-------------------------------------------------------------------------------
# Run all tests
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
