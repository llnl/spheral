#ATS:t1 = test(SELF, "--unitType cgs", label="CGS")
#ATS:t2 = test(SELF, "--unitType mixed", label="Mixed")

from Spheral import *
from SpheralTestUtilities import *


commandLine(
    unitType = "cgs", # "mixed"
)

checksum = 0

#-------------------------------------------------------------------------------
# Check units against reference data
#-------------------------------------------------------------------------------
if unitType == "cgs":
    units = CGS()
    vals = [[units.protonMass, 1.67262192369e-24, "proton mass"],
            [units.electronMass, 9.1093837015e-28, "electron mass"],
            [units.electronCharge, 1.602176634e-19, "electron charge"],
            [units.G, 6.67430e-8, "gravitational constant"],
            [units.c, 2.99792458e10, "speed of light"],
            [units.kB, 1.380649e-16, "Boltzmann constant"],
            [units.molarGasConstant, 8.314462618e7, "molar gas constant"],
            [units.unitMassDensity, 1000.0, "unit mass density"],
            [units.stefanBoltzmannConstant, 5.670374419e-5, "Stefan Boltzmann constant"],
            [units.blackBodyConstant, 7.56573325003e-15, "Radiation constant"],
            [units.planckConstant, 6.62607015e-27, "Planck constant"],
            [units.unitEnergyJ, 1.0e-7, "unit energy"]]
elif unitType == "mixed":
    units = PhysicalConstants(1.e-3, # millimeter
                              1.e-18, # femtograms
                              1.e6, # megasecond
                              1.e9, # gigakelvin
                              1.e-21) # zeptocoulomb
    vals = [[units.protonMass, 1.67262192369e-9, "proton mass"],
            [units.electronMass, 9.1093837015e-13, "electron mass"],
            [units.electronCharge, 1.602176634e2, "electron charge"],
            [units.G, 6.67430e-8, "gravitational constant"],
            [units.c, 2.99792458e17, "speed of light"],
            [units.kB, 1.380649e22, "Boltzmann constant"],
            [units.molarGasConstant, 8.314462618e45, "molar gas constant"],
            [units.unitMassDensity, 1.0e-9, "unit mass density"],
            [units.stefanBoltzmannConstant, 5.670374419e64, "Stefan Boltzmann constant"],
            [units.blackBodyConstant, 7.56573325003e47, "Radiation constant"],
            [units.planckConstant, 6.62607015e-4, "Planck constant"],
            [units.unitEnergyJ, 1.0e-36, "unit energy"]]
    
for val, ref, desc in vals:
    if abs(val - ref) / ref > 1.e-10:
        print(("{} test for {} failed\n\tcalculated: {} \t expected: {}".format(unitType, desc, val, ref)))
        checksum += 1
            
if checksum > 0:
    raise ValueError("number of tests failed: {}".format(checksum))
else:
    print("all tests passed")
