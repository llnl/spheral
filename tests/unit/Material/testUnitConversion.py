#ATS:t1 = test(SELF, label="Unit conversion")

from Spheral import *
from SpheralTestUtilities import *

title("Unit conversion")
commandLine(
    tol = 1.0e-12,
)

def check(a, b, desc):
    err = abs(a - b)
    relErr = 2.0 * err / (abs(a) + abs(b))
    message = "    {}: calc={:.4e}, ref={:.4e}, err={:.4e}, rerr={:.4e}".format(desc, a, b, err, relErr)
    if relErr > tol:
        raise ValueError(message)
    else:
        print(message)

unitsFrom = CGuS()
unitsTo = MKS()

conv = UnitConversion(unitsFrom, unitsTo)

check(conv.length, 1e-2, "length")
check(conv.mass, 1e-3, "mass")
check(conv.time, 1e-6, "time")
check(conv.temperature, 1, "temperature")
check(conv.charge, 1, "charge")
check(conv.massDensity, 1e3, "mass density")
check(conv.specificEnergy, 1e8, "specific energy")
check(conv.energyDensity, 1e11, "energy density")
check(conv.pressure, 1e11, "pressure")
check(conv.bulkModulus, 1e11, "bulk modulus")
check(conv.specificHeat, 1e8, "specific heat")
check(conv.entropy, 1e8, "entropy")
check(conv.speed, 1e4, "speed")
