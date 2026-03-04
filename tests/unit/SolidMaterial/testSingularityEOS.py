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
    rhoWeird = 1.5e-2,
    eps = 2.5e6,
    temp = 405.0,

    # Spiner
    filename = None,
    matid = None,
    printAttr = False,
    spinerMult = 1e3,

    # Testing
    tolerance = 1.0e-3,
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

    check(eos.pressure(rhoWeird, eps), seos.pressure(rhoWeird, eps), "pressure")
    check(eos.temperature(rhoWeird, eps), seos.temperature(rhoWeird, eps), "temperature")
    check(eos.specificThermalEnergy(rhoWeird, temp), seos.specificThermalEnergy(rhoWeird, temp), "specificThermalEnergy")
    check(eos.specificHeat(rhoWeird, temp), seos.specificHeat(rhoWeird, temp), "specificHeat")
    check(eos.soundSpeed(rhoWeird, eps), seos.soundSpeed(rhoWeird, eps), "soundSpeed")
    check(eos.gamma, seos.gamma(rhoWeird, eps), "gamma")
    check(eos.bulkModulus(rhoWeird, eps), seos.bulkModulus(rhoWeird, eps), "bulkModulus")

testGamma()

#-------------------------------------------------------------------------------
# Test spiner EOS
#-------------------------------------------------------------------------------

if printAttr:
    import h5py
    def print_hdf5(name, obj, depth=0):
        indent = "  " * depth
        kind = "Group" if isinstance(obj, h5py.Group) else "Dataset"
        print(f"{indent}{kind}: {name}")
        for k, v in obj.attrs.items():
            print(f"{indent}  ATTR {k} = {v}")

    with h5py.File(filename, "r") as f:
        def visitor(name, obj):
            depth = name.count("/")  # nesting level
            print_hdf5(name, obj, depth)

        f.visititems(visitor)
    quit()

def testSpiner():
    # Check inversion
    print("testing spiner EOS inversion")

    eos = SingularitySpiner(filename, matid, units)
    output("eos")

    highTemp = spinerMult * temp
    check(highTemp, eos.temperature(rho, eos.specificThermalEnergy(rho, highTemp)),
          "inversion of temperature")

    highEps = spinerMult * eps
    check(highEps, eos.specificThermalEnergy(rho, eos.temperature(rho, highEps)),
          "inversion of energy")

    # Check against itself
    print("testing spiner EOS against T version")
    
    eosT = SingularitySpinerT(filename, matid, units)
    output("eosT")

    check(eos.pressure(rho, highEps), eosT.pressure(rho, highEps), "pressure")
    check(eos.temperature(rho, highEps), eosT.temperature(rho, highEps), "temperature")
    check(eos.specificThermalEnergy(rho, highTemp), eosT.specificThermalEnergy(rho, highTemp), "specificThermalEnergy")
    check(eos.specificHeat(rho, highTemp), eosT.specificHeat(rho, highTemp), "specificHeat")
    check(eos.gamma(rho, highEps), eosT.gamma(rho, highEps), "gamma")

    # These are bugged in the current singularity for rhoT
    check(eos.soundSpeed(rho, highEps), eosT.soundSpeed(rho, highEps), "soundSpeed", 1e10)
    check(eos.bulkModulus(rho, highEps), eosT.bulkModulus(rho, highEps), "bulkModulus", 1e10)

if filename is not None and matid is not None:
    testSpiner()
