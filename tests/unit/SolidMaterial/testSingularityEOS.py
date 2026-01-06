#ATS:test(SELF, "", label="Test singularity EOS")

#-------------------------------------------------------------------------------
# Test the singularity EOS
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
import numpy as np
title("Singularity EOS test")

commandLine(
    # Ideal gas
    gam = 5./3.,
    mu = 2.0,

    # Initial conditions
    rho = 1.5e3,
    eps = 2.5e6,
    temp = 405.0,

    # Testing
    tolerance = 1.0e-4,
)

units = MKS()
output("units.c")

def check(a, b, desc, tol = tolerance):
    err = abs(a - b)
    relErr = 2.0 * err / (abs(a) + abs(b))
    message = "    {}: ref={:.4e}, new={:.4e}, err={:.4e}, rerr={:.4e}".format(desc, a, b, err, relErr)
    if relErr > tol:
        raise ValueError(message)
    else:
        print(message)

#-------------------------------------------------------------------------------
# Test the gamma law gas
#-------------------------------------------------------------------------------
def testGamma():
    print("checking gamma law gas")

    weirdUnits = PhysicalConstants(0.03, 2.3, 40.5, 0.7, 7.1)
    eos = GammaLawGas(gam, mu, weirdUnits)
    seos = SingularityGammaLawGas(gam, mu, weirdUnits)
    output("eos")
    output("seos")

    check(eos.pressure(rho, eps), seos.pressure(rho, eps), "pressure")
    check(eos.temperature(rho, eps), seos.temperature(rho, eps), "temperature")
    check(eos.specificThermalEnergy(rho, temp), seos.specificThermalEnergy(rho, temp), "specificThermalEnergy")
    check(eos.specificHeat(rho, temp), seos.specificHeat(rho, temp), "specificHeat")
    check(eos.soundSpeed(rho, eps), seos.soundSpeed(rho, eps), "soundSpeed")
    check(eos.gamma, seos.gamma(rho, eps), "gamma")
    check(eos.bulkModulus(rho, eps), seos.bulkModulus(rho, eps), "bulkModulus")

testGamma()

